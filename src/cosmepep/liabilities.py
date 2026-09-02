"""What will go wrong with the molecule before it ever reaches a fibroblast.

Three failure modes, kept apart because they are fixed by different people. Chemical
liabilities are decided at the sequence and survive into the jar as degradants. Synthesis
liabilities decide whether a contract manufacturer will quote for it at all. Proteolytic
liability decides how long it lasts on skin.

Unlike the descriptors in ``physchem``, most of these terms are *order-sensitive*: Asn-Gly
is a deamidation hotspot and Gly-Asn is not, a Pro at position 2 forms a diketopiperazine
and a Pro at position 5 does not. That matters for what the scoring function is capable of
seeing, so the split is tracked explicitly.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

from .proteases import site_density

#: (name, regex, penalty, rationale). Penalties are on a common 0-1 scale where 1.0 is a
#: liability that would stop development on its own.
CHEMICAL_LIABILITIES: tuple[tuple[str, str, float, str], ...] = (
    ("deamidation_NG", r"N[GS]", 0.60, "Asn-Gly/Ser succinimide route; fastest deamidation motif"),
    ("deamidation_QG", r"QG", 0.25, "slower Gln deamidation"),
    ("acid_cleavage_DP", r"DP", 0.55, "Asp-Pro is the classic acid-labile bond"),
    ("isomerisation_DG", r"D[GS]", 0.40, "Asp-Gly/Ser isoaspartate formation"),
    ("oxidation_M", r"M", 0.35, "methionine sulfoxide; a real shelf-life problem in an aerated emulsion"),
    ("oxidation_W", r"W", 0.30, "tryptophan photo-oxidation, and it discolours the product"),
    ("free_cysteine", r"C", 0.70, "free thiol: disulfide scrambling and metal-catalysed oxidation"),
    ("n_terminal_Q", r"^Q", 0.45, "N-terminal Gln cyclises to pyroglutamate and loses the free amine"),
    ("diketopiperazine", r"^.P", 0.50, "Pro at position 2 drives N-terminal diketopiperazine excision"),
)

#: Solid-phase synthesis pain. Not safety, cost -- but cost is what kills a cosmetic
#: ingredient, since the whole category competes against retinol at a few cents a dose.
SYNTHESIS_LIABILITIES: tuple[tuple[str, str, float, str], ...] = (
    ("beta_sheet_run", r"[VILFY]{3,}", 0.55, "hydrophobic run: on-resin aggregation, deletion sequences"),
    ("polyarginine", r"R{2,}", 0.30, "consecutive Arg; slow couplings and expensive deprotection"),
    ("polyproline", r"P{3,}", 0.35, "polyproline stretch couples badly"),
    ("bulky_pair", r"[WFY][WFY]", 0.25, "adjacent bulky aromatics; steric coupling failure"),
    ("cysteine_pair", r"C.*C", 0.45, "two cysteines: an oxidation state to control at every step"),
)


@dataclass(frozen=True)
class LiabilityReport:
    """Everything wrong with a candidate, itemised."""

    sequence: str
    chemical_hits: list[tuple[str, str]] = field(default_factory=list)
    synthesis_hits: list[tuple[str, str]] = field(default_factory=list)
    chemical_penalty: float = 0.0
    synthesis_penalty: float = 0.0
    protease_site_density: float = 0.0

    @property
    def clean(self) -> bool:
        return not self.chemical_hits and not self.synthesis_hits


def _scan(
    sequence: str, rules: tuple[tuple[str, str, float, str], ...]
) -> tuple[list[tuple[str, str]], float]:
    hits: list[tuple[str, str]] = []
    penalty = 0.0
    for name, pattern, weight, rationale in rules:
        for match in re.finditer(pattern, sequence):
            hits.append((name, match.group(0)))
            penalty += weight
            _ = rationale
    return hits, penalty


def assess(sequence: str) -> LiabilityReport:
    """Full liability report for one sequence."""
    chemical_hits, chemical_penalty = _scan(sequence, CHEMICAL_LIABILITIES)
    synthesis_hits, synthesis_penalty = _scan(sequence, SYNTHESIS_LIABILITIES)
    return LiabilityReport(
        sequence=sequence,
        chemical_hits=chemical_hits,
        synthesis_hits=synthesis_hits,
        chemical_penalty=chemical_penalty,
        synthesis_penalty=synthesis_penalty,
        protease_site_density=site_density(sequence),
    )
