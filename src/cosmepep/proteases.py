"""Skin proteases, as cleavage rules, used twice and for opposite purposes.

Forwards, they enumerate the design space: a peptide that skin can actually liberate from
its own matrix is a plausible endogenous signal, and the set of such fragments is the
"encrypted" matrikine space. Backwards, the same rules are a liability model: a designed
peptide carrying many of these sites will be chewed up in the stratum corneum before it
reaches a fibroblast.

The rules are P1/P1' specificity approximations, not the neural cleavage predictors of the
Manchester group (PROSPER, DeepCleave, MPSC). They are deliberately coarse. A regular
expression over four residues cannot know about tertiary structure or exosite binding, and
treating its output as a prediction of real proteolysis would be a mistake -- what it gives
is an enumeration with the right *flavour*, and a site count that is comparable between two
candidate sequences.
"""

from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class Protease:
    """One protease as a regex whose match end is the scissile bond."""

    name: str
    source: str
    #: Matches up to and including P1; the cleavage falls after the match.
    pattern: str
    note: str

    def sites(self, sequence: str) -> list[int]:
        """Indices *after* which this protease cuts."""
        return [match.end() for match in re.finditer(self.pattern, sequence)]


#: Skin-resident and skin-adjacent proteases. Host enzymes first, then the commensal
#: glutamyl endopeptidases -- V8 and its Staphylococcus epidermidis relatives are abundant
#: on healthy skin and cut after Glu, so they are as much a part of the skin degradome as
#: the MMPs are.
PROTEASES: tuple[Protease, ...] = (
    Protease(
        "MMP-1",
        "fibroblast/keratinocyte",
        r"P.[GA](?=[LIMV])",
        "interstitial collagenase; Pro at P3, Gly/Ala at P1, hydrophobic P1'",
    ),
    Protease(
        "MMP-2/9",
        "fibroblast/neutrophil",
        r"P.[GAN](?=[LIMFV])",
        "gelatinases; denatured-collagen preference, same Pro-x-x-Hy frame",
    ),
    Protease(
        "MMP-12",
        "macrophage",
        r"[PA].[AG](?=[LIMFV])",
        "macrophage elastase; the principal elastin-degrading MMP in photoaged dermis",
    ),
    Protease(
        "ELANE",
        "neutrophil",
        r"[AVI](?![P])",
        "neutrophil elastase; small aliphatic P1, blocked by P1' proline",
    ),
    Protease(
        "CTSG",
        "neutrophil",
        r"[FYW](?![P])",
        "cathepsin G; chymotrypsin-like",
    ),
    Protease(
        "CTSK",
        "fibroblast/osteoclast-type",
        r"[LPV].(?![P])",
        "cathepsin K; the one collagenase that cuts the triple helix at several sites",
    ),
    Protease(
        "KLK5",
        "stratum corneum",
        r"[KR](?![P])",
        "kallikrein-5; trypsin-like desquamation protease, active at skin surface pH",
    ),
    Protease(
        "KLK7",
        "stratum corneum",
        r"[YFML](?![P])",
        "kallikrein-7; chymotryptic desquamation protease",
    ),
    Protease(
        "CMA1",
        "mast cell",
        r"[FYW](?![P])",
        "chymase; released on degranulation, degrades dermal matrix",
    ),
    Protease(
        "GluC/Esp",
        "S. aureus / S. epidermidis",
        r"E(?![P])",
        "commensal glutamyl endopeptidases; cut after Glu, ubiquitous on skin",
    ),
)

#: Fragment lengths worth enumerating. Every cosmetic peptide in commerce is 3-8 residues:
#: shorter and there is no recognition surface, longer and neither synthesis cost nor
#: stratum corneum permeation is tolerable.
MIN_FRAGMENT, MAX_FRAGMENT = 3, 8


@dataclass(frozen=True)
class EncryptedFragment:
    """A peptide that the skin degradome could liberate from a matrix protein."""

    sequence: str
    gene: str
    start: int
    #: Proteases whose rules place a cut at both ends. A fragment released by several
    #: enzymes is more likely to actually appear in tissue than one needing a single
    #: specific enzyme to fire twice.
    n_terminal_enzymes: tuple[str, ...]
    c_terminal_enzymes: tuple[str, ...]

    @property
    def support(self) -> int:
        return len(set(self.n_terminal_enzymes) | set(self.c_terminal_enzymes))


def cut_map(sequence: str) -> dict[int, set[str]]:
    """Map every cut position in a sequence to the proteases that produce it."""
    cuts: dict[int, set[str]] = {}
    for protease in PROTEASES:
        for position in protease.sites(sequence):
            cuts.setdefault(position, set()).add(protease.name)
    return cuts


def enumerate_fragments(
    sequence: str,
    gene: str,
    min_support: int = 2,
) -> list[EncryptedFragment]:
    """All double-cut fragments of a protein within the cosmetic length window.

    ``min_support`` filters on how many distinct proteases bracket the fragment. At the
    default of 2 the enumeration keeps fragments that more than one enzyme could release,
    which is the weakest defensible stand-in for "this actually occurs in tissue".
    """
    cuts = cut_map(sequence)
    # Position 0 and len(sequence) are the protein's own termini, not proteolysis, so they
    # are excluded: a fragment must be cut out, not merely trimmed off one end.
    positions = sorted(cuts)
    fragments: list[EncryptedFragment] = []
    for index, start in enumerate(positions):
        for stop in positions[index + 1 :]:
            span = stop - start
            if span < MIN_FRAGMENT:
                continue
            if span > MAX_FRAGMENT:
                break
            fragment = EncryptedFragment(
                sequence=sequence[start:stop],
                gene=gene,
                start=start,
                n_terminal_enzymes=tuple(sorted(cuts[start])),
                c_terminal_enzymes=tuple(sorted(cuts[stop])),
            )
            if fragment.support >= min_support:
                fragments.append(fragment)
    return fragments


def site_density(sequence: str) -> float:
    """Cut sites per internal peptide bond -- the liability reading of the same rules.

    Two corrections that matter more for short peptides than they would for a protein.
    Distinct positions are counted rather than protease-position pairs, because two enzymes
    cutting the same bond destroy the peptide once and not twice. And only the ``len - 1``
    internal bonds count: a rule matching at the C-terminal position describes a bond that
    is not in this molecule, and including it inflates the density of a tetrapeptide by 25%.
    """
    if len(sequence) < 2:
        return 0.0
    internal = [position for position in cut_map(sequence) if 0 < position < len(sequence)]
    return len(internal) / (len(sequence) - 1)
