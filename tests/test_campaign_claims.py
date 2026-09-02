"""Tests for the claims the README makes about a completed campaign.

These read `results/campaign.json` rather than re-running the search, so they
pin the reported numbers against silent drift without a thirty-generation run.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

CAMPAIGN = Path(__file__).resolve().parents[1] / "results" / "campaign.json"


def _cells(campaign):
    """Every programme-track cell, with its scrambled separation."""
    for programme in campaign["programmes"]:
        for track_name, track in programme["tracks"].items():
            separations = {s["control_kind"]: s for s in track["separations"]}
            if "scrambled" in separations:
                yield programme, track_name, track, separations


@pytest.fixture(scope="module")
def campaign():
    if not CAMPAIGN.exists():
        pytest.skip(f"{CAMPAIGN} missing; run `cosmepep design` first")
    return json.loads(CAMPAIGN.read_text(encoding="utf-8"))


def test_scrambled_controls_hold_composition_exactly_fixed(campaign):
    """The scrambled arm is only meaningful if it really is composition-matched.

    A non-zero distance would mean the control differs from the designs in
    something other than residue order, and the whole comparison would be
    measuring that difference instead.
    """
    for _, _, _, separations in _cells(campaign):
        assert separations["scrambled"]["composition_distance"] == 0.0


def test_naive_controls_flatter_every_programme(campaign):
    """Against a UniProt background every programme looks near-perfect.

    This is the number a pipeline would report if it never built a
    composition-preserving control, and it is nearly uninformative.
    """
    for _, _, _, separations in _cells(campaign):
        naive = separations.get("uniprot_background")
        if naive is None:
            continue
        assert naive["auc"] > 0.98
        assert naive["auc"] > separations["scrambled"]["auc"]


def test_order_sensitive_weight_predicts_surviving_separation(campaign):
    """The README's central quantitative claim.

    Terms that are functions of composition alone return the same value for a
    sequence and its anagram, so they cannot separate designs from scrambles
    even in principle. The share of objective weight that *is* order-sensitive
    should therefore predict the scrambled AUC -- and it does, strongly.
    """
    pytest.importorskip("scipy")
    from scipy import stats

    weights, aucs = [], []
    for programme, _, _, separations in _cells(campaign):
        weights.append(programme["order_sensitive_weight"])
        aucs.append(separations["scrambled"]["auc"])

    assert len(weights) >= 10
    r, p_value = stats.pearsonr(weights, aucs)
    assert r > 0.85
    assert p_value < 0.01


def test_the_weakest_programme_is_reported_not_hidden(campaign):
    """The tyrosinase programme cannot tell a design from its own anagram.

    Pinned because it is the result most tempting to drop, and dropping it
    would leave a set of programmes that all appear to work.
    """
    weakest = min(
        (separations["scrambled"]["auc"] for _, _, _, separations in _cells(campaign))
    )
    assert weakest < 0.65


def test_reference_calibration_is_reported_even_though_it_is_circular(campaign):
    """Marketed peptides sit at the bottom of a pool optimised against them.

    The designs were produced by maximising this composite and the references
    were not, so a percentile near zero restates the objective rather than
    validating it. The field has to stay in the output for that to be visible.
    """
    percentiles = [
        entry["percentile_in_pool"]
        for _, _, track, _ in _cells(campaign)
        for entry in track.get("reference_calibration", [])
    ]
    assert percentiles, "reference calibration was dropped from the report"
    at_floor = sum(1 for value in percentiles if value == 0.0)
    assert at_floor / len(percentiles) > 0.8
