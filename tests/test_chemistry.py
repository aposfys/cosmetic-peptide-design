"""Palmitoylation has to happen on the alpha-amine, and nowhere else."""

from __future__ import annotations

import pytest

from cosmepep.chemistry import LIPIDATION_LOG_KP_GAIN, delivery_gain, potts_guy, profile


def test_palmitoyl_adds_one_chain_not_three() -> None:
    """KTTKS has two lysine epsilon-amines. Acylating one of those is a different molecule."""
    free = profile("KTTKS", "free")
    lipidated = profile("KTTKS", "palmitoyl")
    assert free is not None and lipidated is not None
    assert lipidated.molecular_weight - free.molecular_weight == pytest.approx(238.41, abs=0.1)


def test_arginine_guanidine_is_not_acylated() -> None:
    free = profile("EEMQRR", "free")
    lipidated = profile("EEMQRR", "palmitoyl")
    assert free is not None and lipidated is not None
    assert lipidated.molecular_weight - free.molecular_weight == pytest.approx(238.41, abs=0.1)


def test_amidation_spares_the_glutamate_side_chain() -> None:
    """Ac-EEMQRR-NH2 must lose one oxygen, not three."""
    free = profile("EEMQRR", "free")
    capped = profile("EEMQRR", "acetyl_amide")
    assert free is not None and capped is not None
    # +42.04 for acetyl, -0.98 for OH -> NH2.
    assert capped.molecular_weight - free.molecular_weight == pytest.approx(41.06, abs=0.1)


@pytest.mark.parametrize("sequence", ["KTTKS", "GHK", "EEMQRR", "GEKG", "VGVAPGLA"])
def test_lipidation_gain_is_a_constant(sequence: str) -> None:
    """The documented no-op. If this ever stops holding, the docstring is wrong."""
    assert delivery_gain(sequence) == pytest.approx(LIPIDATION_LOG_KP_GAIN, abs=0.01)


def test_potts_guy_matches_the_published_form() -> None:
    assert potts_guy(0.0, 0.0) == pytest.approx(-2.7)
    assert potts_guy(1.0, 100.0) == pytest.approx(-2.7 + 0.71 - 0.61)


@pytest.mark.parametrize("sequence", ["KTTKS", "GHK", "EEMQRR", "GEKG", "KKIKPLL"])
def test_acetyl_amidation_lowers_permeability_by_a_constant(sequence: str) -> None:
    """Capping is a stability modification, not a delivery one, and it costs 0.55 log units."""
    free = profile(sequence, "free")
    capped = profile(sequence, "acetyl_amide")
    assert free is not None and capped is not None
    assert capped.log_kp - free.log_kp == pytest.approx(-0.55, abs=0.01)
