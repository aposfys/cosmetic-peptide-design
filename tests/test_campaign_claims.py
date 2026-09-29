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

    # Pinned on the encrypted track alone. That track is a ranked enumeration over every
    # fragment with no stochastic search in it, so its correlation is a property of the
    # objective and reproduces exactly. The pooled figure is not safe to assert: the de novo
    # half depends on search budget, and a legitimate `cosmepep design --generations 6
    # --population 60` run drops it from 0.929 to 0.740 -- which would fail this test for a
    # reason that has nothing to do with the claim it is guarding.
    weights, aucs = [], []
    for programme, track_name, _, separations in _cells(campaign):
        if track_name != "encrypted":
            continue
        weights.append(programme["order_sensitive_weight"])
        aucs.append(separations["scrambled"]["auc"])

    assert len(weights) >= 5
    r, p_value = stats.pearsonr(weights, aucs)
    assert r > 0.85
    assert p_value < 0.05


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


def test_tracks_of_a_programme_share_one_order_weight(campaign):
    """Why the pooled correlation must not be given a p-value.

    Order-sensitive weight is a property of the objective, not of a run, so both
    tracks of a programme carry an identical x. The ten pooled points therefore
    sit over five distinct weights, and any test that treats them as ten
    independent observations understates p by about two orders of magnitude.
    """
    weights: dict[str, set[float]] = {}
    for programme, _, _, _ in _cells(campaign):
        weights.setdefault(programme["programme"], set()).add(
            programme["order_sensitive_weight"]
        )
    assert all(len(values) == 1 for values in weights.values())
    assert len(weights) < sum(1 for _ in _cells(campaign))


def test_the_inflated_p_value_is_not_quoted_anywhere(campaign):
    """The pooled p = 0.0001 is the repository's own house error if it survives."""
    root = CAMPAIGN.parents[1]
    for document in (root / "README.md", root / "docs" / "METHODS.md"):
        text = document.read_text(encoding="utf-8")
        assert "p = 0.0001, R" not in text, f"{document.name} quotes the pooled p-value"


def test_every_reference_row_inside_its_pool_is_named_in_the_readme(campaign):
    """The README used to say the informative case happens once. It happens twice.

    A reference landing inside the design pool is the only cell where the calibration
    carries information, so the count and the identities both have to match the campaign
    rather than being carried over from an earlier run.
    """
    inside = [
        (programme["programme"], track_name, reference)
        for programme, track_name, track, _ in _cells(campaign)
        for reference in track["reference_calibration"]
        if reference["percentile_in_pool"] > 0.0
    ]
    assert len(inside) == 2
    readme = (CAMPAIGN.parents[1] / "README.md").read_text(encoding="utf-8")
    for _, _, reference in inside:
        assert reference["sequence"] in readme
        assert f"percentile {reference['percentile_in_pool']:.0f}" in readme
    results = (CAMPAIGN.parent / "RESULTS.md").read_text(encoding="utf-8")
    assert "which happens once" not in results


def test_no_document_claims_the_invariant_terms_are_exactly_half(campaign):
    """Ten of the forty-four declared-invariant cells are not 0.500, so nothing may say
    they are there *by construction*. The terms are invariant; the reported AUC is invariant
    to within floating-point rounding, which is a different sentence.
    """
    flat = [
        term
        for _, _, track, _ in _cells(campaign)
        for term in track["per_term_vs_scrambled"]
        if not term["declared_order_sensitive"]
    ]
    off = [term for term in flat if term["auc_vs_scrambled"] != 0.5]
    assert off, "if every cell is exactly 0.500 this guard can be dropped"
    assert all(abs(term["auc_vs_scrambled"] - 0.5) < 0.01 for term in flat)
    root = CAMPAIGN.parents[1]
    for document in (root / "results" / "RESULTS.md", root / "README.md"):
        assert "0.500 by construction" not in document.read_text(encoding="utf-8")


def test_the_reduced_budget_figures_come_from_a_committed_artifact():
    """The two numbers that used to live only in prose.

    `results/budget_sensitivity.json` is written by `cosmepep sensitivity`, and the README,
    METHODS and `cosmepep audit` all have to quote it rather than a remembered value.
    """
    path = CAMPAIGN.parent / "budget_sensitivity.json"
    if not path.exists():
        pytest.skip(f"{path} missing; run `cosmepep sensitivity` first")
    sensitivity = json.loads(path.read_text(encoding="utf-8"))
    assert (sensitivity["generations"], sensitivity["population_size"]) == (6, 60)
    root = CAMPAIGN.parents[1]
    for key in ("de_novo", "pooled"):
        quoted = f"{sensitivity['correlations'][key]:.3f}"
        assert quoted in (root / "README.md").read_text(encoding="utf-8")
        assert quoted in (root / "docs" / "METHODS.md").read_text(encoding="utf-8")
