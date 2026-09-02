"""Design programmes: what each class of cosmetic peptide is actually trying to do.

Five mechanisms, each written as a set of desirability terms on a common 0-1 scale and
combined as a weighted geometric mean. Geometric rather than arithmetic because these are
conjunctive requirements -- a peptide that cannot be delivered is worthless however well it
fits its target, and an arithmetic mean lets a strong term paper over a zero.

Every term carries an ``order_sensitive`` flag. That flag is the honest part of this module.
A term computed from composition alone assigns the same value to a candidate and to every
shuffle of it, so the weight sitting on permutation-invariant terms is weight that cannot,
even in principle, distinguish a designed sequence from a scrambled one. ``programme_order_weight``
reports that fraction, and it is smaller than the confident framing of this literature
implies.

Where a term encodes chemistry that is genuinely positional -- the Xaa-Xaa-His copper motif
needs the histidine at position three and nowhere else -- it is flagged and weighted
accordingly. Where a mechanism is only understood at the level of composition, no positional
term has been invented to dress it up.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

from .chemistry import profile
from .composition import matrix_enrichment
from .liabilities import assess
from .physchem import describe


def window(value: float, low: float, high: float, slack: float) -> float:
    """1.0 inside [low, high], decaying linearly to 0 over ``slack`` beyond each edge."""
    if low <= value <= high:
        return 1.0
    distance = low - value if value < low else value - high
    return max(0.0, 1.0 - distance / slack)


def rising(value: float, floor: float, ceiling: float) -> float:
    """0 at or below ``floor``, 1 at or above ``ceiling``, linear between."""
    if ceiling == floor:
        return 1.0 if value >= ceiling else 0.0
    return min(1.0, max(0.0, (value - floor) / (ceiling - floor)))


def falling(value: float, floor: float, ceiling: float) -> float:
    """1 at or below ``floor``, 0 at or above ``ceiling``."""
    return 1.0 - rising(value, floor, ceiling)


@dataclass(frozen=True)
class Term:
    """One desirability term in a design programme."""

    name: str
    weight: float
    #: True if the term's value can change when the sequence is shuffled.
    order_sensitive: bool
    rationale: str
    score: Callable[[str], float]


# --------------------------------------------------------------------------------------
# Terms shared by every programme: it has to survive the jar, survive the skin, and get in.
# --------------------------------------------------------------------------------------


#: Potts-Guy anchors per deliverable form. The palmitoyl pair is the 5th and ~75th percentile
#: of the enumerated fragment pool in that form; the free pair is the same two percentiles,
#: which sit 2.55 log units lower because that is exactly what the lipid tail is worth.
#: Ac-peptide-NH2 sits 0.55 log units *below* the free acid, not above it: capping trades a
#: charged terminus for a neutral amide, which Crippen scores as more polar overall. It is a
#: protease-stability modification, not a delivery one, and treating it as though it helped
#: permeation would flatter every acetylated design in the field.
PERMEATION_ANCHORS = {
    "palmitoyl": (-9.7, -4.7),
    "free": (-12.25, -7.25),
    "acetyl_amide": (-12.8, -7.8),
}


def _permeation_for(delivery_form: str):
    """Permeation scored on the form the peptide can actually be shipped in.

    This has to be per-programme rather than global. Palmitoylation acylates the alpha-amine,
    and the alpha-amine is the first of the nitrogens that an ATCUN or GHK-type site uses to
    hold Cu(II) -- so a copper carrier cannot be palmitoylated and remain a copper carrier.
    Scoring one on its lipidated permeability, as an earlier version of this module did,
    credits it for a molecule that would have lost the mechanism the programme selects for.

    Worth stating plainly, because the industry does not: palmitoyl tripeptide-1 is not a
    copper-delivery peptide. Whatever pal-GHK does, it does not do it by carrying Cu(II)
    through an N-terminal site whose N-terminus is capped with a fatty acid.
    """
    floor, ceiling = PERMEATION_ANCHORS[delivery_form]

    def permeation(sequence: str) -> float:
        molecule = profile(sequence, delivery_form)
        if molecule is None:
            return 0.0
        return rising(molecule.log_kp, floor, ceiling)

    return permeation


def _chemical_stability(sequence: str) -> float:
    return falling(assess(sequence).chemical_penalty, 0.0, 1.5)


def _synthesisability(sequence: str) -> float:
    report = assess(sequence)
    length_cost = falling(float(len(sequence)), 6.0, 12.0)
    return 0.5 * falling(report.synthesis_penalty, 0.0, 1.2) + 0.5 * length_cost


def _protease_resistance(sequence: str) -> float:
    """Fewer skin-protease cut sites per residue is better.

    Measured against the only ground truth available -- the marketed peptides -- this term
    has no discriminating power at all. Their site densities run the full 0.0 to 1.0 range:
    GHK, GQPR and TKPR score 0.0, while KVK, LSVD and YR score 1.0, and all six are sold. The
    reason is that products are shipped N-acylated and C-amidated, which defeats the
    exopeptidases these rules stand in for, and the endopeptidase sites are the design
    premise rather than a defect.

    So it is kept at low weight and used as a tie-break among otherwise equal candidates, not
    as a filter. Removing it entirely would also be defensible.
    """
    return falling(assess(sequence).protease_site_density, 0.5, 1.3)


#: A hydrophobic moment is a vector sum around a helical wheel. Below about two turns there
#: is no wheel to sum around, and the number computed for a tripeptide is an artefact of the
#: assumed geometry rather than a property of the molecule. Every term built on the moment is
#: gated on this length.
MIN_HELICAL_LENGTH = 7


def _irritation_margin(sequence: str) -> float:
    """Penalise the cationic-amphipathic combination that permeabilises membranes.

    Inapplicable, not merely absent, below ``MIN_HELICAL_LENGTH``: a tripeptide cannot be
    facially amphipathic, so the term returns 1.0 rather than a number derived from a moment
    that does not exist.
    """
    if len(sequence) < MIN_HELICAL_LENGTH:
        return 1.0
    return falling(describe(sequence).amphipathic_index, 1.0, 4.0)


def developability(delivery_form: str) -> tuple[Term, ...]:
    """Terms every programme carries, with permeation bound to its deliverable form."""
    return (
        Term(
            "permeation", 1.0, False,
            f"Potts-Guy log Kp of the {delivery_form} form", _permeation_for(delivery_form),
        ),
        Term("chemical_stability", 0.8, True, "deamidation, oxidation, DKP, acid-labile bonds", _chemical_stability),
        Term("synthesisability", 0.6, True, "SPPS aggregation, difficult couplings, length", _synthesisability),
        Term("protease_resistance", 0.3, True, "skin degradome cut-site density (uninformative vs the marketed set; tie-break only)", _protease_resistance),
        Term("irritation_margin", 0.5, True, "cationic amphipathicity as a membrane-lysis proxy", _irritation_margin),
    )


# --------------------------------------------------------------------------------------
# Mechanism-specific terms.
# --------------------------------------------------------------------------------------


def _matrix_signal_terms() -> tuple[Term, ...]:
    def fit_length(sequence: str) -> float:
        return window(float(len(sequence)), 3.0, 6.0, 3.0)

    def matrix_like(sequence: str) -> float:
        return rising(matrix_enrichment(sequence), -0.9, 0.6)

    def charge_fit(sequence: str) -> float:
        return window(describe(sequence).charge_surface, 0.0, 2.5, 2.0)

    def hydrophilicity(sequence: str) -> float:
        return window(describe(sequence).gravy, -3.5, -0.3, 1.5)

    def lipidation_handle(sequence: str) -> float:
        """A free alpha-amine that is not immediately lost to diketopiperazine.

        Positional, and the only reason it is: the DKP route needs Pro at position two
        specifically, and a lysine anywhere gives a second acylation site that turns a
        single-product palmitoylation into a mixture.
        """
        if sequence.startswith("P") or (len(sequence) > 1 and sequence[1] == "P"):
            return 0.2
        return 1.0 if sequence.count("K") <= 1 else 0.6

    return (
        Term("length_fit", 1.0, False, "3-6 residues, as every marketed signal peptide is", fit_length),
        Term("matrix_composition", 0.9, False, "per-residue enrichment in the skin ECM panel vs UniProt background", matrix_like),
        Term("charge_fit", 0.7, False, "neutral to +2.5 at skin surface pH", charge_fit),
        Term("hydrophilicity", 0.5, False, "hydrophilic core, delivered by the lipid tail", hydrophilicity),
        Term("lipidation_handle", 0.8, True, "clean single-site palmitoylation, no DKP at position 2", lipidation_handle),
    )


def _cu_carrier_terms() -> tuple[Term, ...]:
    def copper_motif(sequence: str) -> float:
        """Two distinct N-terminal Cu(II) sites, which are positions and not compositions.

        Xaa-Xaa-His (ATCUN) is the four-nitrogen square-planar site: the free alpha-amine,
        the two intervening deprotonated backbone amides, and the His imidazole. Histidine
        has to be third.

        Xaa-His is the GHK site, and it is *not* ATCUN -- a point worth being exact about,
        since GHK is routinely mislabelled as one. With histidine second there is only one
        intervening amide, giving three-nitrogen coordination through the alpha-amine, that
        amide and the imidazole, with the fourth position filled by the His carboxylate or an
        exogenous ligand. Lower affinity than ATCUN, and the motif that the entire copper-
        peptide cosmetic category is actually built on.

        Proline is handled by which nitrogen it removes. At position two of an ATCUN site it
        replaces the backbone NH that has to deprotonate, so the site is abolished. At
        position one it leaves a secondary rather than primary alpha-amine, which still
        coordinates but more weakly, so it is a discount and not a disqualification.
        """
        if len(sequence) < 2 or "H" not in sequence:
            return 0.0
        if len(sequence) >= 3 and sequence[2] == "H":
            if sequence[1] == "P":
                return 0.0
            return 0.75 if sequence[0] == "P" else 1.0
        if sequence[1] == "H":
            return 0.70 if sequence[0] == "P" else 0.9
        return 0.15

    def no_competing_thiol(sequence: str) -> float:
        return 0.0 if "C" in sequence else 1.0

    def fit_length(sequence: str) -> float:
        return window(float(len(sequence)), 3.0, 5.0, 2.0)

    def charge_fit(sequence: str) -> float:
        return window(describe(sequence).charge_surface, 0.5, 2.5, 1.5)

    return (
        Term("copper_motif", 2.0, True, "ATCUN (Xaa-Xaa-His) or GHK-type (Xaa-His) N-terminal Cu(II) site", copper_motif),
        Term("no_competing_thiol", 1.0, False, "Cys reduces Cu(II) and destroys the complex", no_competing_thiol),
        Term("length_fit", 0.8, False, "3-5 residues, as GHK is", fit_length),
        Term("charge_fit", 0.6, False, "cationic at skin surface pH", charge_fit),
    )


def _snare_competitor_terms() -> tuple[Term, ...]:
    def acidic_n_terminus(sequence: str) -> float:
        """An acidic block at the N-terminus, as in the SNAP-25 mimic.

        Positional by construction: the reported mechanism is competition for an N-terminal
        position in the SNARE bundle, so acidic residues at the far end of the peptide are
        not the same molecule at all.
        """
        head = sequence[:2]
        return sum(1 for residue in head if residue in "ED") / 2.0

    def basic_c_terminus(sequence: str) -> float:
        tail = sequence[-2:]
        return sum(1 for residue in tail if residue in "KR") / 2.0

    def fit_length(sequence: str) -> float:
        return window(float(len(sequence)), 6.0, 8.0, 3.0)

    def near_neutral(sequence: str) -> float:
        return window(describe(sequence).charge_dermal, -1.0, 1.0, 2.0)

    def helix_propensity(sequence: str) -> float:
        """Chou-Fasman helix formers net of breakers. SNARE motifs are helical bundles.

        Breakers are subtracted rather than merely not counted. An earlier version counted
        formers only, and the optimiser promptly returned EEPLKK -- acidic head, basic tail,
        and a proline sitting in the middle of a sequence whose entire premise is that it
        forms a helix. Not counting a defect is not the same as penalising it.
        """
        formers = sum(sequence.count(residue) for residue in "AELMQKRH")
        breakers = sum(sequence.count(residue) for residue in "PG")
        return rising((formers - breakers) / len(sequence), 0.0, 0.6)

    return (
        Term("acidic_n_terminus", 1.6, True, "Glu/Asp at positions 1-2, as in EEMQRR", acidic_n_terminus),
        Term("basic_c_terminus", 0.9, True, "basic tail, as in EEMQRR", basic_c_terminus),
        Term("length_fit", 0.8, False, "6-8 residues", fit_length),
        Term("near_neutral", 0.6, False, "near-neutral at dermal pH", near_neutral),
        Term("helix_propensity", 0.7, False, "helical character for SNARE bundle mimicry", helix_propensity),
    )


def _tyrosinase_terms() -> tuple[Term, ...]:
    def copper_chelator(sequence: str) -> float:
        """Tyrosinase is a type-3 dicopper enzyme; His is the ligand that matters."""
        return rising(float(sequence.count("H")), 0.0, 2.0)

    def substrate_mimicry(sequence: str) -> float:
        """Aromatic residues that can occupy the L-tyrosine pocket."""
        return rising(sum(sequence.count(residue) for residue in "YF") / len(sequence), 0.0, 0.4)

    def fit_length(sequence: str) -> float:
        return window(float(len(sequence)), 3.0, 6.0, 3.0)

    def cationic(sequence: str) -> float:
        return window(describe(sequence).charge_surface, 0.0, 3.0, 2.0)

    return (
        Term("copper_chelator", 1.4, False, "His content for the dicopper active site", copper_chelator),
        Term("substrate_mimicry", 1.2, False, "Tyr/Phe for the substrate pocket", substrate_mimicry),
        Term("length_fit", 0.7, False, "small enough to enter the active-site channel", fit_length),
        Term("cationic", 0.5, False, "cationic at skin surface pH", cationic),
    )


def _barrier_amp_terms() -> tuple[Term, ...]:
    def cationicity(sequence: str) -> float:
        return window(describe(sequence).charge_surface, 3.0, 6.0, 2.0)

    def amphipathicity(sequence: str) -> float:
        """Hydrophobic moment: one of only two descriptors here that reads sequence order.

        Zero below ``MIN_HELICAL_LENGTH``, where the moment is undefined in substance -- for
        this programme that is a real disqualification rather than an inapplicable term,
        since a peptide too short to be amphipathic cannot work by this mechanism.
        """
        if len(sequence) < MIN_HELICAL_LENGTH:
            return 0.0
        return rising(describe(sequence).hydrophobic_moment, 0.2, 0.6)

    def selectivity(sequence: str) -> float:
        """Bounded hydrophobicity. Past a point, activity stops being selective for
        bacterial membranes and starts being haemolysis -- on skin, stinging."""
        return window(describe(sequence).gravy, -1.0, 0.8, 1.2)

    def fit_length(sequence: str) -> float:
        return window(float(len(sequence)), 7.0, 12.0, 4.0)

    return (
        Term("cationicity", 1.3, False, "+3 to +6 at skin surface pH", cationicity),
        Term("amphipathicity", 1.3, True, "hydrophobic moment on a helical wheel", amphipathicity),
        Term("selectivity", 1.0, False, "bounded hydrophobicity as a haemolysis proxy", selectivity),
        Term("length_fit", 0.6, False, "7-12 residues for a facially amphipathic turn", fit_length),
    )


@dataclass(frozen=True)
class Programme:
    """One design campaign: a mechanism, its terms, and the developability it shares."""

    key: str
    title: str
    target: str
    reference: str
    terms: tuple[Term, ...]
    #: The chemical form this programme's output can actually be shipped in. Not a
    #: preference -- for ``cu_carrier`` it is forced by the coordination chemistry.
    delivery_form: str = "palmitoyl"

    @property
    def all_terms(self) -> tuple[Term, ...]:
        return self.terms + developability(self.delivery_form)

    def score(self, sequence: str) -> dict[str, float]:
        """Per-term desirabilities plus the weighted geometric mean under ``composite``."""
        values = {term.name: term.score(sequence) for term in self.all_terms}
        total_weight = sum(term.weight for term in self.all_terms)
        # A single zero must propagate. Terms are floored at a small epsilon only so that a
        # near-miss on one axis is ranked above a hard failure on it, rather than both
        # collapsing to an indistinguishable zero.
        logarithm = 0.0
        for term in self.all_terms:
            logarithm += term.weight * _log(max(values[term.name], 1e-6))
        values["composite"] = pow(2.718281828459045, logarithm / total_weight)
        return values

    @property
    def order_weight(self) -> float:
        """Fraction of total weight carried by terms that can see sequence order."""
        total = sum(term.weight for term in self.all_terms)
        ordered = sum(term.weight for term in self.all_terms if term.order_sensitive)
        return ordered / total


def _log(value: float) -> float:
    from math import log

    return log(value)


PROGRAMMES: tuple[Programme, ...] = (
    Programme(
        "matrix_signal",
        "Procollagen-stimulating matrikine",
        "dermal fibroblast; procollagen I/III and fibronectin transcription",
        "KTTKS, GEKG, GPKG",
        _matrix_signal_terms(),
    ),
    Programme(
        "cu_carrier",
        "Copper-carrier tripeptide",
        "Cu(II) delivery via an N-terminal ATCUN or GHK-type chelation site",
        "GHK / GHK-Cu",
        _cu_carrier_terms(),
        delivery_form="free",
    ),
    Programme(
        "snare_competitor",
        "SNARE-competing relaxant",
        "SNAP-25 / SNARE assembly at the neuromuscular junction",
        "EEMQRR (acetyl hexapeptide-8)",
        _snare_competitor_terms(),
        delivery_form="acetyl_amide",
    ),
    Programme(
        "tyrosinase_modulator",
        "Tyrosinase-modulating brightener",
        "tyrosinase type-3 dicopper active site",
        "no well-attested peptide reference; the weakest-grounded programme here",
        _tyrosinase_terms(),
    ),
    Programme(
        "barrier_amp",
        "Selective anti-C. acnes peptide",
        "Cutibacterium acnes membrane, at skin surface pH",
        "host defence peptide chemotype",
        _barrier_amp_terms(),
    ),
)

PROGRAMME_BY_KEY = {programme.key: programme for programme in PROGRAMMES}
