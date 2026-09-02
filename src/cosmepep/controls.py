"""The nulls, without which none of the scores above mean anything.

Four control families, ordered from the weakest test to the strongest:

``uniprot_background``
    Random sequence at UniProt-wide amino-acid frequencies. The naive null, and the wrong
    one here -- the skin matrix is 15% glycine and 10% proline, so a matrix-derived peptide
    is separable from this on composition alone and any scoring function clears the bar
    without knowing anything about peptides.
``ecm_background``
    Random sequence at the pooled composition of the panel itself. Removes the compositional
    freebie that ``uniprot_background`` hands over.
``natural_fragment``
    Real matrix sequence, same length, sampled at random positions from the panel. Real
    residue order, real local grammar, selected for nothing. This is the control that asks
    whether the pipeline found anything, or merely re-described the proteome it read.
``scrambled``
    The design's own residues in a different order. Composition is held exactly fixed, so
    separation here can only come from order.

The last one is decisive and its expected result is known in advance from the term flags in
``objectives``: any weight sitting on a permutation-invariant term contributes exactly
nothing to it. Measuring it anyway is the point -- it converts a design's headline score into
a statement about how much of that score is a claim about the sequence rather than about its
amino-acid counts.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

from .composition import BACKGROUND_FREQUENCIES
from .ecm import fetch_panel, residue_frequencies

CONTROL_FAMILIES = ("scrambled", "natural_fragment", "ecm_background", "uniprot_background")


@dataclass(frozen=True)
class ControlSet:
    """One control family generated against a set of designs."""

    kind: str
    sequences: tuple[str, ...]
    note: str


def _weighted(frequencies: dict[str, float], length: int, rng: random.Random) -> str:
    residues = list(frequencies)
    weights = [frequencies[residue] for residue in residues]
    return "".join(rng.choices(residues, weights=weights, k=length))


def build_controls(
    designs: list[str],
    *,
    seed: int = 0,
    replicates: int = 5,
) -> dict[str, ControlSet]:
    """Generate every control family, ``replicates`` decoys per design.

    Length is matched design by design rather than in aggregate, because several terms in
    every programme are explicit length windows and an unmatched length control would be
    separated by those alone.
    """
    rng = random.Random(seed)
    panel = fetch_panel()
    ecm_frequencies = residue_frequencies(panel)
    joined = [protein.sequence for protein in panel]

    scrambled: list[str] = []
    natural: list[str] = []
    ecm_random: list[str] = []
    uniprot_random: list[str] = []

    for design in designs:
        length = len(design)
        for _ in range(replicates):
            residues = list(design)
            rng.shuffle(residues)
            scrambled.append("".join(residues))

            protein = rng.choice(joined)
            start = rng.randrange(0, max(1, len(protein) - length))
            natural.append(protein[start : start + length])

            ecm_random.append(_weighted(ecm_frequencies, length, rng))
            uniprot_random.append(_weighted(BACKGROUND_FREQUENCIES, length, rng))

    return {
        "scrambled": ControlSet(
            "scrambled",
            tuple(scrambled),
            "same residues, different order; separation here is separation on sequence",
        ),
        "natural_fragment": ControlSet(
            "natural_fragment",
            tuple(natural),
            "real panel sequence at matched length, selected for nothing",
        ),
        "ecm_background": ControlSet(
            "ecm_background",
            tuple(ecm_random),
            "random draw at the panel's own composition",
        ),
        "uniprot_background": ControlSet(
            "uniprot_background",
            tuple(uniprot_random),
            "random draw at UniProt-wide composition; the naive null",
        ),
    }


def composition_distance(left: list[str], right: list[str]) -> float:
    """Total-variation distance between two sequence sets' pooled compositions.

    Reported so that a control's claim to be composition-matched is a measurement rather
    than an assertion. The scrambled family must come out at exactly zero; if it does not,
    the shuffling is broken.
    """

    def pooled(sequences: list[str]) -> dict[str, float]:
        counts = dict.fromkeys(BACKGROUND_FREQUENCIES, 0)
        for sequence in sequences:
            for residue in sequence:
                if residue in counts:
                    counts[residue] += 1
        total = max(sum(counts.values()), 1)
        return {residue: count / total for residue, count in counts.items()}

    left_pooled, right_pooled = pooled(left), pooled(right)
    return 0.5 * sum(
        abs(left_pooled[residue] - right_pooled[residue]) for residue in BACKGROUND_FREQUENCIES
    )
