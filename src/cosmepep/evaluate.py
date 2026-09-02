"""Separation between designs and their controls, with bootstrap intervals.

The reportable output of a design run is not a top-ten table. It is how far the objective
separates what it selected from each null, and what that separation is worth once the
interval is drawn. A composite that ranks designs far above UniProt-random draws and no
better than chance above their own scrambles has told you about amino-acid frequencies.

The dataclass and the rank-based AUC follow ``pepdesign.evaluate`` so numbers from the two
projects can be read side by side.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from dataclasses import dataclass
from statistics import mean, pstdev


@dataclass(frozen=True)
class Separation:
    """How well a score separates designs from one control family."""

    control_kind: str
    n_design: int
    n_control: int
    mean_design: float
    mean_control: float
    auc: float
    auc_ci_low: float
    auc_ci_high: float
    effect_size: float
    threshold_at_5pct_control: float
    design_recall_at_threshold: float
    #: Total-variation distance between the two sets' pooled compositions. Zero for a
    #: scrambled control by construction; anything else is a bug in the shuffling.
    composition_distance: float


def roc_auc(positive: Sequence[float], negative: Sequence[float]) -> float:
    """AUC as the probability a random positive outscores a random negative.

    Rank-based so ties take half credit exactly. These scores are heavily tied -- a bounded
    desirability objective puts many candidates on a plateau at 1.0 -- and a trapezoidal
    implementation would round every one of those ties in a single direction.
    """
    if not positive or not negative:
        raise ValueError("need both populations to compute an AUC")
    combined = sorted([(value, 1) for value in positive] + [(value, 0) for value in negative])
    index = 0
    rank_sum_positive = 0.0
    while index < len(combined):
        stop = index
        while stop + 1 < len(combined) and combined[stop + 1][0] == combined[index][0]:
            stop += 1
        average_rank = (index + stop) / 2.0 + 1.0
        for position in range(index, stop + 1):
            if combined[position][1] == 1:
                rank_sum_positive += average_rank
        index = stop + 1
    count_positive, count_negative = len(positive), len(negative)
    return (rank_sum_positive - count_positive * (count_positive + 1) / 2) / (
        count_positive * count_negative
    )


def bootstrap_auc(
    positive: Sequence[float],
    negative: Sequence[float],
    *,
    draws: int = 2000,
    seed: int = 0,
) -> tuple[float, float]:
    """Percentile bootstrap interval on the AUC, resampling both populations."""
    rng = random.Random(seed)
    values = []
    for _ in range(draws):
        resampled_positive = [rng.choice(positive) for _ in positive]
        resampled_negative = [rng.choice(negative) for _ in negative]
        values.append(roc_auc(resampled_positive, resampled_negative))
    values.sort()
    return values[int(0.025 * draws)], values[int(0.975 * draws)]


def cohens_d(positive: Sequence[float], negative: Sequence[float]) -> float:
    """Pooled-standard-deviation effect size, reported next to the AUC.

    A large AUC on a negligible effect is possible when distributions are tight, and reads
    misleadingly on its own.
    """
    spread_positive, spread_negative = pstdev(positive), pstdev(negative)
    pooled = ((spread_positive**2 + spread_negative**2) / 2) ** 0.5
    if pooled == 0:
        return 0.0
    return (mean(positive) - mean(negative)) / pooled


def separate(
    control_kind: str,
    design_scores: Sequence[float],
    control_scores: Sequence[float],
    composition_distance: float,
    *,
    seed: int = 0,
) -> Separation:
    """Full separation report for one control family."""
    ranked_controls = sorted(control_scores, reverse=True)
    cutoff_index = max(0, int(0.05 * len(ranked_controls)) - 1)
    threshold = ranked_controls[cutoff_index] if ranked_controls else 0.0
    recall = sum(1 for score in design_scores if score >= threshold) / len(design_scores)
    low, high = bootstrap_auc(design_scores, control_scores, seed=seed)
    return Separation(
        control_kind=control_kind,
        n_design=len(design_scores),
        n_control=len(control_scores),
        mean_design=mean(design_scores),
        mean_control=mean(control_scores),
        auc=roc_auc(design_scores, control_scores),
        auc_ci_low=low,
        auc_ci_high=high,
        effect_size=cohens_d(design_scores, control_scores),
        threshold_at_5pct_control=threshold,
        design_recall_at_threshold=recall,
        composition_distance=composition_distance,
    )
