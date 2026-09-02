"""The campaign: enumerate, evolve, then try hard to show that it did not work.

Each programme is run on two tracks. The *encrypted* track enumerates what the skin
degradome could cut out of the panel proteins and ranks it. The *de novo* track runs the same
objective with panel substrings forbidden, so its output is sequence that does not exist in
the human matrix. Both are then put through the four control families and reported with
intervals.

The ordering is deliberate. Designs are produced first and evaluated second, against nulls
that were fixed before any of the numbers were seen, and every separation is reported
including the ones that come out at chance.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from .chemistry import profile, smiles
from .composition import matrix_enrichment
from .controls import build_controls, composition_distance
from .evaluate import Separation, separate
from .generate import Candidate, diverse_subset, evolve, levenshtein, plateau, seed_library
from .liabilities import assess
from .objectives import PROGRAMMES, Programme
from .physchem import describe
from .reference import KNOWN_PEPTIDES

RESULTS = Path(__file__).resolve().parents[2] / "results"

#: How many top-scoring candidates form the population that goes into the null comparison.
#: Large enough for the AUC interval to mean something, small enough to still be "the
#: designs" rather than the whole enumeration.
EVALUATION_POOL = 120
#: How many candidates are carried into the reported table per track.
SHORTLIST = 8


@dataclass
class DesignCard:
    """Everything needed to decide whether to make a candidate."""

    sequence: str
    programme: str
    track: str
    provenance: str
    source_gene: str | None
    composite: float
    terms: dict[str, float]
    length: int
    free_mw: float
    free_clogp: float
    free_log_kp: float
    palmitoyl_mw: float
    palmitoyl_clogp: float
    palmitoyl_log_kp: float
    #: The form this programme can actually ship, and its properties in that form. For
    #: cu_carrier these are the free-acid numbers, because a palmitoylated copper carrier is
    #: not a copper carrier.
    delivery_form: str
    delivery_mw: float
    delivery_log_kp: float
    charge_ph55: float
    charge_ph74: float
    isoelectric_point: float
    gravy: float
    matrix_enrichment: float
    protease_site_density: float
    chemical_liabilities: list[str]
    synthesis_liabilities: list[str]
    nearest_known: str
    nearest_known_distance: int
    recommended_form: str
    palmitoyl_smiles: str | None


def _recommended_form(sequence: str, programme_key: str) -> str:
    """Which chemical form to make, from the programme and the liabilities.

    Not cosmetic advice dressed up as chemistry: N-acylation removes the free alpha-amine,
    which is exactly the group the copper programmes need, so the two are mutually exclusive
    and the choice is forced rather than preferred.
    """
    if programme_key == "cu_carrier":
        return (
            "free N-terminus, C-amidation only. Acylation caps the alpha-amine that holds "
            "Cu(II), so a palmitoylated version of this peptide is a different mechanism."
        )
    if sequence.startswith("Q"):
        return "acetyl / C-amide (blocks pyroglutamate formation at the N-terminus)"
    if programme_key == "snare_competitor":
        return "acetyl / C-amide, as acetyl hexapeptide-8 is"
    return "N-palmitoyl / C-amide"


def build_card(
    candidate: Candidate, track: str, delivery_form: str = "palmitoyl"
) -> DesignCard:
    sequence = candidate.sequence
    free = profile(sequence, "free")
    lipidated = profile(sequence, "palmitoyl")
    delivered = profile(sequence, delivery_form)
    sequence_profile = describe(sequence)
    liabilities = assess(sequence)
    nearest = min(
        ((known.sequence, levenshtein(sequence, known.sequence)) for known in KNOWN_PEPTIDES),
        key=lambda pair: pair[1],
    )
    assert free is not None and lipidated is not None and delivered is not None
    return DesignCard(
        sequence=sequence,
        programme=candidate.programme,
        track=track,
        provenance=candidate.provenance,
        source_gene=candidate.source_gene,
        composite=candidate.composite,
        terms={name: round(value, 4) for name, value in candidate.terms.items()},
        length=len(sequence),
        free_mw=round(free.molecular_weight, 2),
        free_clogp=round(free.clogp, 2),
        free_log_kp=round(free.log_kp, 2),
        palmitoyl_mw=round(lipidated.molecular_weight, 2),
        palmitoyl_clogp=round(lipidated.clogp, 2),
        palmitoyl_log_kp=round(lipidated.log_kp, 2),
        delivery_form=delivery_form,
        delivery_mw=round(delivered.molecular_weight, 2),
        delivery_log_kp=round(delivered.log_kp, 2),
        charge_ph55=round(sequence_profile.charge_surface, 2),
        charge_ph74=round(sequence_profile.charge_dermal, 2),
        isoelectric_point=round(sequence_profile.isoelectric_point, 2),
        gravy=round(sequence_profile.gravy, 2),
        matrix_enrichment=round(matrix_enrichment(sequence), 3),
        protease_site_density=round(liabilities.protease_site_density, 3),
        chemical_liabilities=[f"{name}:{hit}" for name, hit in liabilities.chemical_hits],
        synthesis_liabilities=[f"{name}:{hit}" for name, hit in liabilities.synthesis_hits],
        nearest_known=nearest[0],
        nearest_known_distance=nearest[1],
        recommended_form=_recommended_form(sequence, candidate.programme),
        palmitoyl_smiles=smiles(sequence, "palmitoyl"),
    )


def run_programme(
    programme: Programme,
    *,
    generations: int = 30,
    population_size: int = 250,
    seed: int = 0,
) -> dict:
    """One programme, both tracks, with controls."""
    avoid = tuple(known.sequence for known in KNOWN_PEPTIDES)

    encrypted = seed_library(programme, limit=EVALUATION_POOL * 3)
    encrypted = [
        candidate
        for candidate in encrypted
        if all(levenshtein(candidate.sequence, other) >= 2 for other in avoid)
    ]
    de_novo = evolve(
        programme,
        generations=generations,
        population_size=population_size,
        seed=seed,
        avoid=avoid,
        require_de_novo=True,
    )

    tracks = {"encrypted": encrypted, "de_novo": de_novo}
    report: dict = {
        "programme": programme.key,
        "title": programme.title,
        "target": programme.target,
        "reference": programme.reference,
        "delivery_form": programme.delivery_form,
        "order_sensitive_weight": round(programme.order_weight, 4),
        "terms": [
            {
                "name": term.name,
                "weight": term.weight,
                "order_sensitive": term.order_sensitive,
                "rationale": term.rationale,
            }
            for term in programme.all_terms
        ],
        "tracks": {},
    }

    for track_name, candidates in tracks.items():
        pool = candidates[:EVALUATION_POOL]
        if not pool:
            continue
        top = plateau(candidates)
        # Selection band, not a ranking. Where the objective has a plateau, that plateau is
        # the band; where it does not, the band is everything within 5% of the best score,
        # floored at 40 candidates. Taking the top eight by composite instead returns a
        # single mutational family -- the barrier_amp track returns eight variants of KKIKP --
        # which is a property of the optimiser's last generation and not a set of designs.
        if len(top) > SHORTLIST:
            band = top
        else:
            best = candidates[0].composite
            band = [c for c in candidates if c.composite >= best * 0.95][: max(40, SHORTLIST)]
            band = band or candidates[:SHORTLIST]
        shortlist = diverse_subset(band, SHORTLIST)

        design_sequences = [candidate.sequence for candidate in pool]
        design_scores = [candidate.composite for candidate in pool]
        controls = build_controls(design_sequences, seed=seed, replicates=5)

        separations: list[Separation] = []
        for kind, control_set in controls.items():
            control_scores = [
                programme.score(sequence)["composite"] for sequence in control_set.sequences
            ]
            separations.append(
                separate(
                    kind,
                    design_scores,
                    control_scores,
                    composition_distance(design_sequences, list(control_set.sequences)),
                    seed=seed,
                )
            )

        # Per-term separation against the scrambled control only. This is where the
        # order-sensitivity flags get checked against behaviour instead of being trusted.
        scrambled = controls["scrambled"].sequences
        per_term: list[dict] = []
        for term in programme.all_terms:
            design_term = [term.score(sequence) for sequence in design_sequences]
            control_term = [term.score(sequence) for sequence in scrambled]
            from .evaluate import roc_auc

            per_term.append(
                {
                    "term": term.name,
                    "weight": term.weight,
                    "declared_order_sensitive": term.order_sensitive,
                    "auc_vs_scrambled": round(roc_auc(design_term, control_term), 4),
                }
            )

        # Where the marketed peptides of this class fall inside the design pool. A pool that
        # cannot reach the compounds people already sell is optimising the wrong thing; a
        # pool that buries them is worth looking at twice before believing it.
        reference_calibration = []
        for known in KNOWN_PEPTIDES:
            if known.mechanism_class != programme.key:
                continue
            known_score = programme.score(known.sequence)["composite"]
            below = sum(1 for score in design_scores if score < known_score)
            reference_calibration.append(
                {
                    "sequence": known.sequence,
                    "inci": known.inci,
                    "evidence": known.evidence,
                    "composite": round(known_score, 4),
                    "percentile_in_pool": round(100.0 * below / len(design_scores), 1),
                }
            )

        report["tracks"][track_name] = {
            "reference_calibration": reference_calibration,
            "selection_band": len(band),
            "n_evaluated": len(pool),
            "plateau_size": len(top),
            "plateau_score": round(top[0].composite, 4) if top else None,
            "shortlist": [
                asdict(build_card(candidate, track_name, programme.delivery_form))
                for candidate in shortlist
            ],
            "separations": [asdict(separation) for separation in separations],
            "per_term_vs_scrambled": per_term,
        }

    return report


def run_all(
    *, generations: int = 30, population_size: int = 250, seed: int = 0, out: Path = RESULTS
) -> dict:
    """Every programme, written to ``results/campaign.json``."""
    campaign = {
        "seed": seed,
        "generations": generations,
        "population_size": population_size,
        "evaluation_pool": EVALUATION_POOL,
        "programmes": [
            run_programme(
                programme, generations=generations, population_size=population_size, seed=seed
            )
            for programme in PROGRAMMES
        ],
    }
    out.mkdir(parents=True, exist_ok=True)
    (out / "campaign.json").write_text(json.dumps(campaign, indent=2))
    return campaign
