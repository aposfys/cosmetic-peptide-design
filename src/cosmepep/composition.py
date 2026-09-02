"""How matrix-like a peptide's residues are, measured against the panel rather than guessed.

An earlier version of this scored "collagen character" as the fraction of Gly, Pro, Ala and
Lys. That hand-picked set assigns zero to LSVD, which is one of the two tetrapeptide
matrikines that went through a split-face clinical study -- so the term was rejecting a
validated compound for not looking like the author's idea of collagen. The replacement is
derived: each residue is scored by how enriched it is in the skin ECM panel relative to the
UniProt-wide background, and the peptide takes the mean over its residues.

It remains a composition statistic, and therefore permutation-invariant. It is not a claim
about activity. It says only that the residues would not be out of place in a matrix
fragment, which is the most a composition can say.
"""

from __future__ import annotations

from functools import lru_cache
from math import log2

from .ecm import fetch_panel, residue_frequencies

#: UniProt-wide amino-acid frequencies, the denominator for enrichment.
BACKGROUND_FREQUENCIES: dict[str, float] = {
    "A": 0.0825,
    "R": 0.0553,
    "N": 0.0406,
    "D": 0.0545,
    "C": 0.0137,
    "Q": 0.0393,
    "E": 0.0675,
    "G": 0.0707,
    "H": 0.0227,
    "I": 0.0596,
    "L": 0.0966,
    "K": 0.0584,
    "M": 0.0242,
    "F": 0.0386,
    "P": 0.0470,
    "S": 0.0656,
    "T": 0.0534,
    "W": 0.0108,
    "Y": 0.0292,
    "V": 0.0687,
}


@lru_cache(maxsize=1)
def ecm_log_enrichment() -> dict[str, float]:
    """log2(f_ECM / f_background) per residue, from the cached panel."""
    frequencies = residue_frequencies(fetch_panel())
    return {
        residue: log2(max(frequencies[residue], 1e-6) / BACKGROUND_FREQUENCIES[residue])
        for residue in BACKGROUND_FREQUENCIES
    }


@lru_cache(maxsize=200_000)
def matrix_enrichment(sequence: str) -> float:
    """Mean per-residue log2 enrichment in the skin ECM panel."""
    table = ecm_log_enrichment()
    return sum(table[residue] for residue in sequence) / len(sequence)
