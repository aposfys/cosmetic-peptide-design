"""Turn a campaign into Markdown, leading with the null rather than the shortlist."""

from __future__ import annotations

import json
from pathlib import Path

RESULTS = Path(__file__).resolve().parents[2] / "results"


def _separation_table(separations: list[dict]) -> list[str]:
    order = {
        "scrambled": 0,
        "natural_fragment": 1,
        "ecm_background": 2,
        "uniprot_background": 3,
    }
    rows = [
        "| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |",
        "| --- | ---: | --- | ---: | ---: |",
    ]
    for separation in sorted(
        separations, key=lambda s: order.get(s["control_kind"], 9), reverse=True
    ):
        marker = "**" if separation["control_kind"] == "scrambled" else ""
        rows.append(
            f"| {marker}{separation['control_kind']}{marker} "
            f"| {marker}{separation['auc']:.3f}{marker} "
            f"| [{separation['auc_ci_low']:.3f}, {separation['auc_ci_high']:.3f}] "
            f"| {separation['effect_size']:.2f} "
            f"| {separation['composition_distance']:.3f} |"
        )
    return rows


def write(campaign: dict, out: Path = RESULTS) -> Path:
    lines: list[str] = [
        "# Results",
        "",
        f"Seed {campaign['seed']}, {campaign['generations']} generations, "
        f"population {campaign['population_size']}, "
        f"{campaign['evaluation_pool']} candidates per track into the null comparison.",
        "",
        "Every AUC below is the same designs scored by the same function. Only the control",
        "changes. The spread between the first row and the last of each table is the size of",
        "the choice a pipeline makes when it reports a single number.",
        "",
    ]
    for programme in campaign["programmes"]:
        lines += [
            f"## {programme['title']} (`{programme['programme']}`)",
            "",
            f"**Target.** {programme['target']}  ",
            f"**Reference compound.** {programme['reference']}  ",
            f"**Deliverable form.** {programme['delivery_form']} — permeation is scored on "
            "this form, not on a lipidated one the mechanism could not survive.  ",
            f"**Weight on order-sensitive terms.** {programme['order_sensitive_weight']:.0%} — "
            "the rest cannot distinguish a design from its own scramble even in principle.",
            "",
        ]
        for track_name, track in programme["tracks"].items():
            lines += [
                f"### {track_name} track",
                "",
                f"Plateau: {track['plateau_size']} candidates tied at "
                f"{track['plateau_score']}; the shortlist is a max-min diversity selection "
                "from that plateau, not a ranking.",
                "",
            ]
            lines += _separation_table(track["separations"])
            lines.append("")
            if track["reference_calibration"]:
                inside = [
                    reference
                    for reference in track["reference_calibration"]
                    if reference["percentile_in_pool"] > 0.0
                ]
                if inside:
                    verdict = (
                        "The informative case is a reference landing *inside* the pool. "
                        + (
                            "; ".join(
                                f"{reference['inci']} reaches percentile "
                                f"{reference['percentile_in_pool']:.0f} here"
                                for reference in inside
                            )
                            + "."
                        )
                    )
                else:
                    verdict = (
                        "Every reference in this pool sits at the floor, so this table carries "
                        "no information about the designs at all."
                    )
                lines += [
                    "Marketed peptides of this class, placed in the same pool. A percentile of 0 "
                    f"is **not** a finding: the pool is the top {track['n_evaluated']} candidates "
                    "by this composite and the reference compounds were not selected on it, so "
                    "the comparison restates the objective rather than validating it. "
                    + verdict,
                    "",
                ]
                lines += [
                    "| Peptide | INCI | Evidence | Composite | Percentile in pool |",
                    "| --- | --- | --- | ---: | ---: |",
                ]
                for reference in track["reference_calibration"]:
                    lines.append(
                        f"| {reference['sequence']} | {reference['inci']} | {reference['evidence']} "
                        f"| {reference['composite']:.3f} | {reference['percentile_in_pool']:.0f} |"
                    )
                lines.append("")
            form = programme["delivery_form"]
            lines += [
                f"| Sequence | Provenance | Source | Composite | MW ({form}) | log Kp ({form}) | q(5.5) | Liabilities |",
                "| --- | --- | --- | ---: | ---: | ---: | ---: | --- |",
            ]
            for card in track["shortlist"]:
                liabilities = (
                    ", ".join(card["chemical_liabilities"] + card["synthesis_liabilities"])
                    or "—"
                )
                lines.append(
                    f"| `{card['sequence']}` | {card['provenance']} | {card['source_gene'] or '—'} "
                    f"| {card['composite']:.3f} | {card['delivery_mw']:.0f} "
                    f"| {card['delivery_log_kp']:.2f} | {card['charge_ph55']:+.2f} | {liabilities} |"
                )
            lines.append("")
            ordered = [
                t for t in track["per_term_vs_scrambled"] if t["declared_order_sensitive"]
            ]
            flat = [
                t for t in track["per_term_vs_scrambled"] if not t["declared_order_sensitive"]
            ]
            worst = max((abs(t["auc_vs_scrambled"] - 0.5) for t in flat), default=0.0)
            residual = (
                "sit at exactly 0.500 here"
                if worst == 0.0
                else (
                    f"sit at 0.500 to within {worst:.4f} here. That residual is not a term "
                    "reading sequence order; its two sources are set out in "
                    "[METHODS](../docs/METHODS.md) section 5"
                )
            )
            lines += [
                "Per-term AUC against the scrambled control. Terms declared composition-only "
                f"({', '.join(t['term'] for t in flat)}) return the same value for a sequence "
                f"and its anagram, so they {residual}. The separation such as it is comes from "
                f"{', '.join(t['term'] for t in ordered)}.",
                "",
                "| Term | Weight | Order-sensitive | AUC vs scrambled |",
                "| --- | ---: | --- | ---: |",
            ]
            for term in track["per_term_vs_scrambled"]:
                lines.append(
                    f"| {term['term']} | {term['weight']:.1f} | "
                    f"{'yes' if term['declared_order_sensitive'] else 'no'} | "
                    f"{term['auc_vs_scrambled']:.3f} |"
                )
            lines.append("")

    out.mkdir(parents=True, exist_ok=True)
    path = out / "RESULTS.md"
    path.write_text("\n".join(lines))
    return path


def load(out: Path = RESULTS) -> dict:
    return json.loads((out / "campaign.json").read_text())
