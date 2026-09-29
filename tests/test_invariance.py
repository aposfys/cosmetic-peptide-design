"""The permutation-invariance claim, asserted rather than trusted.

``physchem`` states that most of its descriptors cannot see sequence order. That statement
governs how every AUC against the scrambled control should be read, so it is checked here
against actual shuffles instead of being taken on the author's word.
"""

from __future__ import annotations

import itertools
import random

import pytest

from cosmepep.chemistry import profile
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
        if any(
            getattr(describe(s), field) != pytest.approx(reference, abs=1e-9)
            for s in shuffles(sequence)
        ):
            changed = True
    assert changed, f"{field} was declared order-sensitive but never changed under shuffling"


def test_invariant_descriptors_differ_between_anagrams_only_at_rounding() -> None:
    """Why the declared-invariant per-term AUCs are 0.500 and not exactly 0.500.

    Molecular weight, GRAVY and the Boman index are sums taken in sequence order, so an
    anagram accumulates identical addends in a different order and lands a unit or two in
    the last place away. ``roc_auc`` gives ties half credit by comparing floats for equality,
    and these objectives clamp, so most of a pool sits on one large tie. A difference of
    1e-16 splits that tie and moves an AUC that should be exactly 0.500. The effect is real
    and it is bounded, which is the point of pinning it here.
    """
    sequence = "KKIKPLL"
    anagrams = sorted({"".join(p) for p in itertools.permutations(sequence)})
    reference = describe(sequence)
    inexact = set()
    for anagram in anagrams:
        described = describe(anagram)
        for field in INVARIANT:
            difference = abs(getattr(described, field) - getattr(reference, field))
            assert difference < 1e-12, f"{field} moved by {difference} on {anagram}"
            if difference:
                inexact.add(field)
    assert inexact, "no descriptor drifted, so the rounding explanation no longer applies"
    assert inexact <= {"molecular_weight", "gravy", "boman"}


def test_free_amine_permeation_reads_an_n_terminal_proline() -> None:
    """The one declared-invariant term with a real order dependence.

    Crippen log P types an atom from its neighbours, so a free N-terminal proline is a ring
    secondary amine where an internal proline is an amide. That makes ``log Kp`` of the
    free-amine form depend on proline position, which is why ``cu_carrier`` is the one
    programme whose ``permeation`` cell stays off 0.500 after rounding. Capping the
    N-terminus removes the dependence entirely.
    """
    anagrams = sorted({"".join(p) for p in itertools.permutations("GHKPA")})

    def log_kp(sequence: str, form: str) -> float:
        molecule = profile(sequence, form)
        assert molecule is not None
        return molecule.log_kp

    free = {round(log_kp(anagram, "free"), 9) for anagram in anagrams}
    assert len(free) == 2, "the free-amine form no longer splits on proline position"
    assert max(free) - min(free) == pytest.approx(0.0579, abs=5e-4)
    leading = {round(log_kp(a, "free"), 9) for a in anagrams if a.startswith("P")}
    assert len(leading) == 1 and min(free) in leading

    for form in ("palmitoyl", "acetyl_amide"):
        capped = {round(log_kp(anagram, form), 9) for anagram in anagrams}
        assert len(capped) == 1, f"{form} should cap the N-terminus and remove the dependence"

    without_proline = {round(log_kp(a, "free"), 9) for a in {"GHKVA", "HKVAG", "AGVKH"}}
    assert len(without_proline) == 1


def test_liability_penalties_read_order() -> None:
    """Asn-Gly is a deamidation hotspot and Gly-Asn is not; the scanner must agree."""
    assert assess("NG").chemical_penalty > assess("GN").chemical_penalty
