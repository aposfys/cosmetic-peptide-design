"""Two routes into the design space, and the distinction between them.

The first route is enumeration: every 3-8mer the skin degradome could cut out of the panel
proteins. These are *encrypted* peptides, not de novo ones -- they exist in the human
proteome already and the only novelty is that nobody has looked at them. The second route is
a genetic algorithm seeded from those fragments, which mutates, splices and extends them
under the programme objective until the sequences are no longer substrings of any panel
protein. Those are de novo.

Keeping the two labelled apart matters commercially and scientifically. An encrypted peptide
has a prior -- the body makes it, so something has probably evolved to read it -- and a
freedom-to-operate problem, since a natural human sequence is hard to claim. A de novo
sequence has neither.
"""

from __future__ import annotations

import random
from dataclasses import dataclass, field

from .ecm import EcmProtein, fetch_panel
from .objectives import Programme
from .proteases import enumerate_fragments

RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


@dataclass
class Candidate:
    """One designed sequence with its provenance and its scores."""

    sequence: str
    programme: str
    #: "encrypted" if the sequence occurs verbatim in the panel, "de_novo" otherwise.
    provenance: str
    composite: float
    terms: dict[str, float] = field(default_factory=dict)
    source_gene: str | None = None
    generation: int = 0


class PanelIndex:
    """Membership test for "does this sequence occur in the human skin matrix?"."""

    def __init__(self, proteins: list[EcmProtein]) -> None:
        # A single joined string with a separator that cannot appear in a sequence, so a
        # substring search cannot straddle two proteins and invent a fragment.
        self._joined = "|".join(protein.sequence for protein in proteins)
        self._by_gene = {protein.gene: protein.sequence for protein in proteins}

    def contains(self, sequence: str) -> bool:
        return sequence in self._joined

    def source_gene(self, sequence: str) -> str | None:
        for gene, protein in self._by_gene.items():
            if sequence in protein:
                return gene
        return None


def levenshtein(left: str, right: str) -> int:
    """Edit distance, used to keep designs away from sequences that already exist."""
    if len(left) < len(right):
        left, right = right, left
    previous = list(range(len(right) + 1))
    for index, left_residue in enumerate(left, start=1):
        current = [index]
        for position, right_residue in enumerate(right, start=1):
            current.append(
                min(
                    previous[position] + 1,
                    current[position - 1] + 1,
                    previous[position - 1] + (left_residue != right_residue),
                )
            )
        previous = current
    return previous[-1]


def seed_library(programme: Programme, limit: int = 400) -> list[Candidate]:
    """Best enumerated matrix fragments for one programme."""
    panel = fetch_panel()
    index = PanelIndex(panel)
    seen: dict[str, str] = {}
    for protein in panel:
        for fragment in enumerate_fragments(protein.sequence, protein.gene):
            seen.setdefault(fragment.sequence, fragment.gene)

    scored: list[Candidate] = []
    for sequence, gene in seen.items():
        terms = programme.score(sequence)
        scored.append(
            Candidate(
                sequence=sequence,
                programme=programme.key,
                provenance="encrypted",
                composite=terms["composite"],
                terms=terms,
                source_gene=gene,
            )
        )
    scored.sort(key=lambda candidate: candidate.composite, reverse=True)
    _ = index
    return scored[:limit]


def _mutate(sequence: str, rng: random.Random, min_length: int, max_length: int) -> str:
    operation = rng.random()
    # Precedence made explicit, not changed: `and` already bound tighter, so a
    # sequence at the minimum length substitutes rather than shrinking further.
    if operation < 0.65 or (len(sequence) <= min_length and operation < 0.85):
        position = rng.randrange(len(sequence))
        return sequence[:position] + rng.choice(RESIDUES) + sequence[position + 1 :]
    if operation < 0.85 and len(sequence) < max_length:
        position = rng.randrange(len(sequence) + 1)
        return sequence[:position] + rng.choice(RESIDUES) + sequence[position:]
    if len(sequence) > min_length:
        position = rng.randrange(len(sequence))
        return sequence[:position] + sequence[position + 1 :]
    return sequence


def _crossover(left: str, right: str, rng: random.Random, max_length: int) -> str:
    cut_left = rng.randrange(1, len(left))
    cut_right = rng.randrange(0, len(right))
    return (left[:cut_left] + right[cut_right:])[:max_length]


def evolve(
    programme: Programme,
    *,
    generations: int = 40,
    population_size: int = 300,
    seed: int = 0,
    min_length: int = 3,
    max_length: int = 12,
    avoid: tuple[str, ...] = (),
    min_edit_distance: int = 2,
    require_de_novo: bool = False,
) -> list[Candidate]:
    """Genetic algorithm over the programme objective, seeded from matrix fragments.

    ``avoid`` and ``min_edit_distance`` keep the search away from sequences that already
    exist as products: a design one substitution from KTTKS is a KTTKS analogue with a
    freedom-to-operate problem, not a new peptide, and without this constraint the optimiser
    walks straight back onto the marketed compounds -- which is reassuring about the
    objective and useless as an output.

    ``require_de_novo`` forbids sequences that occur in the panel proteome. It is needed
    because the objective is bounded at 1.0 and the enumerated matrix fragments already reach
    it, so an unconstrained search has no gradient pushing it off the natural manifold and
    returns encrypted peptides indefinitely. Running the same objective with and without this
    constraint is the experiment that says what novelty costs.
    """
    rng = random.Random(seed)
    index = PanelIndex(fetch_panel())

    def admissible(sequence: str) -> bool:
        if not min_length <= len(sequence) <= max_length:
            return False
        if require_de_novo and index.contains(sequence):
            return False
        return all(levenshtein(sequence, other) >= min_edit_distance for other in avoid)

    seeds = seed_library(programme, limit=population_size)
    population = [candidate.sequence for candidate in seeds if admissible(candidate.sequence)]
    while len(population) < population_size:
        length = rng.randint(min_length, max_length)
        population.append("".join(rng.choice(RESIDUES) for _ in range(length)))

    scores: dict[str, dict[str, float]] = {}

    def score_of(sequence: str) -> float:
        if sequence not in scores:
            scores[sequence] = programme.score(sequence)
        return scores[sequence]["composite"]

    best: dict[str, Candidate] = {}
    for generation in range(generations):
        ranked = sorted(population, key=score_of, reverse=True)
        for sequence in ranked[: population_size // 4]:
            if sequence not in best:
                best[sequence] = Candidate(
                    sequence=sequence,
                    programme=programme.key,
                    provenance="encrypted" if index.contains(sequence) else "de_novo",
                    composite=score_of(sequence),
                    terms=scores[sequence],
                    source_gene=index.source_gene(sequence),
                    generation=generation,
                )
        elite = ranked[: population_size // 5]
        children: list[str] = list(elite)
        while len(children) < population_size:
            if rng.random() < 0.35:
                child = _crossover(rng.choice(elite), rng.choice(elite), rng, max_length)
            else:
                child = _mutate(rng.choice(elite), rng, min_length, max_length)
            if admissible(child):
                children.append(child)
        population = children

    return sorted(best.values(), key=lambda candidate: candidate.composite, reverse=True)


def plateau(candidates: list[Candidate], tolerance: float = 1e-9) -> list[Candidate]:
    """Every candidate scoring within ``tolerance`` of the best.

    Reported rather than hidden. A bounded desirability objective built from windows and
    ramps does not have a unique optimum, it has a plateau, and the size of that plateau is
    the honest statement of how much the objective actually decides. Presenting the first
    eight rows of a table where four hundred sequences share the top score would be a
    presentation artefact, not a result.
    """
    if not candidates:
        return []
    top = max(candidate.composite for candidate in candidates)
    return [candidate for candidate in candidates if top - candidate.composite <= tolerance]


def diverse_subset(candidates: list[Candidate], count: int) -> list[Candidate]:
    """Greedy max-min edit-distance selection, used to choose within a plateau.

    When the objective cannot rank a set, sequence diversity is a defensible basis for
    choosing among it -- it maximises the chance that at least one member survives contact
    with an assay -- and it is at least an explicit criterion rather than whichever order the
    dictionary happened to be in.
    """
    if not candidates:
        return []
    chosen = [candidates[0]]
    pool = candidates[1:]
    while len(chosen) < count and pool:
        best_index, best_distance = 0, -1
        for index, candidate in enumerate(pool):
            distance = min(levenshtein(candidate.sequence, taken.sequence) for taken in chosen)
            if distance > best_distance:
                best_index, best_distance = index, distance
        chosen.append(pool.pop(best_index))
    return chosen
