"""Command line for the campaign."""

from __future__ import annotations

import argparse
import json
import sys

from .ecm import fetch_panel, residue_frequencies
from .objectives import PROGRAMMES
from .proteases import PROTEASES, enumerate_fragments
from .report import load, write
from .run import run_all


def _panel() -> None:
    panel = fetch_panel()
    print(f"{'gene':10s} {'accession':10s} {'len':>6s}  compartment")
    for protein in panel:
        print(f"{protein.gene:10s} {protein.accession:10s} {protein.length:6d}  {protein.compartment}")
    print(f"\n{len(panel)} proteins, {sum(p.length for p in panel):,} residues")
    frequencies = residue_frequencies(panel)
    top = sorted(frequencies.items(), key=lambda pair: -pair[1])[:6]
    print("most abundant: " + ", ".join(f"{residue} {value:.1%}" for residue, value in top))


def _degradome() -> None:
    panel = fetch_panel()
    fragments: dict[str, str] = {}
    for protein in panel:
        for fragment in enumerate_fragments(protein.sequence, protein.gene):
            fragments.setdefault(fragment.sequence, fragment.gene)
    print(f"{len(PROTEASES)} protease rules over {len(panel)} proteins")
    print(f"{len(fragments):,} unique 3-8mer fragments with >=2 enzymes bracketing them")
    lengths: dict[int, int] = {}
    for sequence in fragments:
        lengths[len(sequence)] = lengths.get(len(sequence), 0) + 1
    for length in sorted(lengths):
        print(f"  {length}-mer  {lengths[length]:6,d}")


def _programmes() -> None:
    for programme in PROGRAMMES:
        print(f"{programme.key:22s} {programme.order_weight:5.0%} order-sensitive weight  {programme.title}")
        for term in programme.all_terms:
            flag = "order" if term.order_sensitive else "comp "
            print(f"    {flag}  w={term.weight:4.1f}  {term.name:22s} {term.rationale}")
        print()


def _audit() -> None:
    """The relationship between what a programme's terms can see and what it achieves."""
    campaign = load()
    rows = []
    for programme in campaign["programmes"]:
        for track_name, track in programme["tracks"].items():
            by_kind = {s["control_kind"]: s["auc"] for s in track["separations"]}
            rows.append(
                (programme["programme"], track_name, programme["order_sensitive_weight"],
                 by_kind["scrambled"], by_kind["uniprot_background"])
            )
    print(f"{'programme':22s} {'track':10s} {'order-w':>8s} {'scrambled':>10s} {'naive null':>11s} {'gap':>7s}")
    for name, track_name, weight, scrambled, naive in rows:
        print(f"{name:22s} {track_name:10s} {weight:8.2f} {scrambled:10.3f} {naive:11.3f} {naive - scrambled:7.3f}")
    weights = [row[2] for row in rows]
    scrambles = [row[3] for row in rows]
    mean_weight = sum(weights) / len(weights)
    mean_scramble = sum(scrambles) / len(scrambles)
    covariance = sum((w - mean_weight) * (s - mean_scramble) for w, s in zip(weights, scrambles))
    spread = (
        sum((w - mean_weight) ** 2 for w in weights) * sum((s - mean_scramble) ** 2 for s in scrambles)
    ) ** 0.5
    print(
        f"\nPearson r(order-sensitive weight, AUC vs scrambled) = {covariance / spread:.3f}  "
        f"over n={len(rows)} programme-tracks.\nHow well an objective separates designs from "
        "their own scrambles is predicted by how much\nof its weight sits on terms that can read "
        "order at all -- which is knowable by inspection,\nbefore anything is run."
    )


def _evaluate() -> None:
    campaign = load()
    print(f"{'programme':22s} {'track':10s} {'scrambled':>10s} {'natural':>10s} {'ecm-bg':>10s} {'uniprot':>10s}")
    for programme in campaign["programmes"]:
        for track_name, track in programme["tracks"].items():
            by_kind = {s["control_kind"]: s["auc"] for s in track["separations"]}
            print(
                f"{programme['programme']:22s} {track_name:10s} "
                f"{by_kind['scrambled']:10.3f} {by_kind['natural_fragment']:10.3f} "
                f"{by_kind['ecm_background']:10.3f} {by_kind['uniprot_background']:10.3f}"
            )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="cosmepep", description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("panel", help="the skin ECM proteins the design space is drawn from")
    subparsers.add_parser("degradome", help="enumerate protease-liberated fragments")
    subparsers.add_parser("programmes", help="design programmes and their scoring terms")
    design = subparsers.add_parser("design", help="run the full campaign")
    design.add_argument("--generations", type=int, default=30)
    design.add_argument("--population", type=int, default=250)
    design.add_argument("--seed", type=int, default=0)
    subparsers.add_parser("evaluate", help="separation table from an existing campaign")
    subparsers.add_parser("audit", help="order-sensitive weight against achieved separation")
    subparsers.add_parser("report", help="write results/RESULTS.md from campaign.json")

    args = parser.parse_args(argv)
    if args.command == "panel":
        _panel()
    elif args.command == "degradome":
        _degradome()
    elif args.command == "programmes":
        _programmes()
    elif args.command == "design":
        campaign = run_all(
            generations=args.generations, population_size=args.population, seed=args.seed
        )
        path = write(campaign)
        print(f"wrote results/campaign.json and {path}")
    elif args.command == "evaluate":
        _evaluate()
    elif args.command == "audit":
        _audit()
    elif args.command == "report":
        print(f"wrote {write(load())}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
