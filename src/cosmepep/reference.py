"""Cosmetic peptides already in commerce, used as the positive reference set.

This is not a training set and nothing here is fitted to it. It serves two purposes: it
locates the region of sequence space that the industry has converged on, and it gives the
scoring functions something to be sanity-checked against -- a scheme that ranks KTTKS and
GHK below random matrix fragments is measuring the wrong thing.

The ``evidence`` field is deliberate. Cosmetic peptide claims range from split-face
randomised studies to a single supplier brochure, and pooling them would import the
marketing into the model. Only ``clinical`` and ``in_vitro`` entries are used when the
reference set is treated as positives.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class KnownPeptide:
    """A marketed or literature-reported cosmetic peptide."""

    sequence: str
    inci: str
    mechanism_class: str
    origin: str
    #: "clinical" (controlled human study published), "in_vitro" (cell or enzyme data),
    #: "supplier" (claim rests on manufacturer literature).
    evidence: str
    note: str


KNOWN_PEPTIDES: tuple[KnownPeptide, ...] = (
    KnownPeptide(
        "KTTKS", "Palmitoyl pentapeptide-4", "matrix_signal",
        "COL1A1 C-propeptide", "clinical",
        "The archetype. Free KTTKS does not permeate; the palmitoyl ester is the product.",
    ),
    KnownPeptide(
        "GHK", "Tripeptide-1 / copper tripeptide-1", "cu_carrier",
        "COL1A2", "in_vitro",
        "Xaa-His three-nitrogen Cu(II) site, not the Xaa-Xaa-His ATCUN motif it is often called.",
    ),
    KnownPeptide(
        "GQPR", "Palmitoyl tetrapeptide-7", "matrix_signal",
        "immunoglobulin G (rigin fragment)", "in_vitro",
        "Paired with pal-KTTKS in Matrixyl 3000; anti-inflammatory framing.",
    ),
    KnownPeptide(
        "EEMQRR", "Acetyl hexapeptide-8 (Argireline)", "snare_competitor",
        "SNAP-25 N-terminal mimic", "clinical",
        "Reported to compete for a position in the SNARE complex and blunt contraction.",
    ),
    KnownPeptide(
        "YAGFL", "Pentapeptide-18 (Leuphasyl)", "neuro_other",
        "Leu-enkephalin analogue", "supplier",
        "Opioid-receptor route, not SNARE competition -- a different mechanism in the same "
        "marketing category, and scored against the SNARE programme it would fail correctly.",
    ),
    KnownPeptide(
        "KDVY", "Acetyl tetrapeptide-5", "matrix_signal",
        "designed", "supplier", "Anti-oedema/periorbital claim, ACE-inhibitory framing.",
    ),
    KnownPeptide(
        "FVAPFP", "Hexapeptide-11", "matrix_signal",
        "yeast-derived", "supplier", "Firming and senescence claims.",
    ),
    KnownPeptide(
        "KVK", "Palmitoyl tripeptide-5", "matrix_signal",
        "TSP-1 mimetic framing", "supplier", "Claimed TGF-beta activation route to collagen.",
    ),
    KnownPeptide(
        "GEKG", "Tetrapeptide (GEKG)", "matrix_signal",
        "collagen-derived", "in_vitro", "Reported to raise COL1A1 and hyaluronan synthase.",
    ),
    KnownPeptide(
        "VGVAPG", "Palmitoyl hexapeptide-12", "matrix_signal",
        "ELN hexapeptide repeat", "in_vitro",
        "The cautionary case: a real elastokine that is also chemotactic and MMP-inducing.",
    ),
    KnownPeptide(
        "GPKG", "GPKG tetrapeptide matrikine", "matrix_signal",
        "predicted matrix cleavage product", "clinical",
        "From in-silico cleavage prediction through to a split-face study (Jariwala 2024).",
    ),
    KnownPeptide(
        "LSVD", "LSVD tetrapeptide matrikine", "matrix_signal",
        "predicted matrix cleavage product", "clinical",
        "Companion tetrapeptide from the same screen.",
    ),
    KnownPeptide(
        "TKPR", "Tuftsin (tetrapeptide-1)", "matrix_signal",
        "IgG Fc (leukokinin)", "in_vitro", "Immunomodulatory tetrapeptide used for repair claims.",
    ),
    KnownPeptide(
        "YR", "Acetyl dipeptide-1", "neuro_other",
        "designed", "supplier", "Sensory-neuron signalling claim; again not a SNARE competitor.",
    ),
    KnownPeptide(
        "VW", "Dipeptide-2", "matrix_signal",
        "designed", "supplier", "ACE-inhibitory dipeptide, periorbital puffiness claim.",
    ),
)


def positives(min_evidence: str = "in_vitro") -> list[KnownPeptide]:
    """Reference peptides at or above an evidence tier."""
    tiers = {"clinical": 2, "in_vitro": 1, "supplier": 0}
    floor = tiers[min_evidence]
    return [peptide for peptide in KNOWN_PEPTIDES if tiers[peptide.evidence] >= floor]
