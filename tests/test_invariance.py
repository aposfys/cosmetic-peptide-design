"""The permutation-invariance claim, asserted rather than trusted.

``physchem`` states that most of its descriptors cannot see sequence order. That statement
governs how every AUC against the scrambled control should be read, so it is checked here
against actual shuffles instead of being taken on the author's word.
"""

from __future__ import annotations

import random

import pytest

from cosmepep.liabilities import assess
from cosmepep.physchem import describe

SEQUENCES = ("KTTKSGE", "EEMQRRA", "VGVAPGL", "GHKPTND", "WFYKRIV", "NGDPMCA")

INVARIANT = (
    "molecular_weight",
    "charge_surface",
    "charge_dermal",
    "isoelectric_point",
    "gravy",
    "boman",
    "aliphatic_index",
    "collagen_character",
    "redox_character",
)
ORDER_SENSITIVE = ("hydrophobic_moment", "instability_index")


def shuffles(sequence: str, count: int = 80) -> list[str]:
    rng = random.Random(hash(sequence) & 0xFFFF)
    return ["".join(rng.sample(sequence, len(sequence))) for _ in range(count)]


@pytest.mark.parametrize("field", INVARIANT)
def test_descriptor_is_permutation_invariant(field: str) -> None:
    for sequence in SEQUENCES:
        reference = getattr(describe(sequence), field)
        for shuffled in shuffles(sequence):
            assert getattr(describe(shuffled), field) == pytest.approx(reference, abs=1e-9)


@pytest.mark.parametrize("field", ORDER_SENSITIVE)
def test_descriptor_reads_order(field: str) -> None:
    changed = False
    for sequence in SEQUENCES:
        reference = getattr(describe(sequence), field)
        if any(getattr(describe(s), field) != pytest.approx(reference, abs=1e-9) for s in shuffles(sequence)):
            changed = True
    assert changed, f"{field} was declared order-sensitive but never changed under shuffling"


def test_liability_penalties_read_order() -> None:
    """Asn-Gly is a deamidation hotspot and Gly-Asn is not; the scanner must agree."""
    assert assess("NG").chemical_penalty > assess("GN").chemical_penalty
