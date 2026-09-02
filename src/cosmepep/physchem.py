"""Sequence-level descriptors, computed at the pH the skin actually has.

Almost every peptide property table in the literature is quoted at pH 7. The skin surface
sits at pH 4.7-5.5 -- the acid mantle -- and only becomes neutral in the viable epidermis.
A peptide's net charge, and therefore its solubility, its partitioning into the stratum
corneum lipids and its electrostatic fit to a target, are different numbers in those two
places. Both are computed here and the difference between them is informative: a candidate
that is cationic on the surface and neutral at its target has a delivery problem that no
amount of lipidation fixes.

A warning that applies to most of this module: molecular weight, both charges, the
isoelectric point, GRAVY, the Boman index, the aliphatic index and the two composition
fractions are functions of amino-acid *composition* only. They are permutation-invariant --
identical for a sequence and for any shuffle of it -- so no threshold on them can separate a
design from its own scramble, however well it separates designs from random draws. Only the
hydrophobic moment (a vector sum over an assumed helical geometry) and the instability index
(a sum over dipeptide weights) read sequence order at all. ``tests/test_invariance.py``
asserts this residue by residue rather than trusting the claim; ``evaluate`` measures what
it costs.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import peptides

#: Skin surface, under the acid mantle.
SURFACE_PH = 5.5
#: Viable epidermis and dermis, where a signal peptide meets its receptor.
DERMAL_PH = 7.4


@dataclass(frozen=True)
class SequenceProfile:
    """Physicochemical descriptors for one sequence."""

    sequence: str
    length: int
    molecular_weight: float
    charge_surface: float
    charge_dermal: float
    isoelectric_point: float
    gravy: float
    hydrophobic_moment: float
    boman: float
    instability_index: float
    aliphatic_index: float
    #: Fraction of residues that are Gly, Pro, Hyp-capable or Ala -- the collagen signature.
    collagen_character: float
    #: Fraction aromatic (F, W, Y) plus His: the radical-scavenging and metal-handling set.
    redox_character: float

    @property
    def amphipathic_index(self) -> float:
        """Hydrophobic moment scaled by net surface charge.

        The membrane-lytic combination, and therefore the irritation proxy: cationic
        amphipathic peptides are antimicrobial because they permeabilise bilayers, and a
        cosmetic that permeabilises bilayers is a stinging complaint.
        """
        return self.hydrophobic_moment * max(self.charge_surface, 0.0)


@lru_cache(maxsize=200_000)
def describe(sequence: str) -> SequenceProfile:
    """Descriptor set for one sequence."""
    peptide = peptides.Peptide(sequence)
    length = len(sequence)
    collagen_residues = sum(sequence.count(residue) for residue in "GPAK")
    redox_residues = sum(sequence.count(residue) for residue in "FWYH")
    return SequenceProfile(
        sequence=sequence,
        length=length,
        molecular_weight=peptide.molecular_weight(),
        charge_surface=peptide.charge(pH=SURFACE_PH),
        charge_dermal=peptide.charge(pH=DERMAL_PH),
        isoelectric_point=peptide.isoelectric_point(),
        gravy=peptide.hydrophobicity(scale="KyteDoolittle"),
        hydrophobic_moment=peptide.hydrophobic_moment(angle=100),
        boman=peptide.boman(),
        instability_index=peptide.instability_index(),
        aliphatic_index=peptide.aliphatic_index(),
        collagen_character=collagen_residues / length,
        redox_character=redox_residues / length,
    )
