"""Controls have to be what they claim to be."""

from __future__ import annotations

import pytest

from cosmepep.controls import build_controls, composition_distance
from cosmepep.evaluate import roc_auc

DESIGNS = ["GLPGPP", "THPG", "GKPG", "HHGP", "KTTKSG"]


def test_scrambled_control_is_exactly_composition_matched() -> None:
    controls = build_controls(DESIGNS, seed=0, replicates=20)
    distance = composition_distance(DESIGNS, list(controls["scrambled"].sequences))
    assert distance == pytest.approx(0.0, abs=1e-12)


def test_scrambled_control_matches_length_design_by_design() -> None:
    controls = build_controls(DESIGNS, seed=0, replicates=3)
    lengths = sorted(len(s) for s in controls["scrambled"].sequences)
    expected = sorted(len(d) for d in DESIGNS for _ in range(3))
    assert lengths == expected


def test_naive_null_is_further_away_than_the_matched_one() -> None:
    controls = build_controls(DESIGNS, seed=0, replicates=40)
    naive = composition_distance(DESIGNS, list(controls["uniprot_background"].sequences))
    matched = composition_distance(DESIGNS, list(controls["ecm_background"].sequences))
    assert naive > matched


def test_auc_gives_ties_half_credit() -> None:
    assert roc_auc([1.0, 1.0, 1.0], [1.0, 1.0, 1.0]) == pytest.approx(0.5)
    assert roc_auc([1.0, 1.0], [0.0, 0.0]) == pytest.approx(1.0)
