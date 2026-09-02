# Methods

## 1. Design space: the skin ECM degradome

Nineteen human proteins are fetched from UniProt and cached (`data/cache/ecm_proteins.json`):
the fibrillar collagens of the dermis, the elastic fibre, the dermal-epidermal junction and
basement membrane, two small leucine-rich proteoglycans, and filaggrin. 38,211 residues.

Ten protease specificity rules stand in for the skin degradome — MMP-1, the gelatinases,
MMP-12, neutrophil elastase, cathepsins G and K, the desquamation kallikreins KLK5 and KLK7,
mast-cell chymase, and the commensal glutamyl endopeptidases (V8/Esp) that are abundant on
healthy skin and cut after Glu. Each rule is a regular expression whose match end is the
scissile bond.

A **fragment** is any 3–8mer bounded at both ends by a predicted cut, with at least two
distinct enzymes bracketing it. That yields **40,026 unique fragments**.

These rules are coarse on purpose. They are not the Manchester group's neural cleavage
predictors (PROSPER, DeepCleave, MPSC), and a regex over four residues knows nothing about
tertiary structure or exosites. What they give is an enumeration with the right flavour and a
site count comparable between two candidates — not a prediction of proteolysis in tissue.

## 2. Two provenance tracks

| Track | Definition | Why it is kept separate |
| --- | --- | --- |
| `encrypted` | Occurs verbatim in the panel proteome | Has a prior — the body makes it — and a freedom-to-operate problem |
| `de_novo` | Panel substrings forbidden during search | No prior, no prior art |

The de novo constraint is structural, not cosmetic. The objective is bounded at 1.0 and the
enumerated fragments already reach it, so an unconstrained optimiser has no gradient pushing
it off the natural manifold and returns encrypted peptides indefinitely. Running the same
objective with and without the constraint is the experiment that says what novelty costs.

## 3. Scoring

Five programmes, each a set of desirability terms on [0, 1] combined as a **weighted
geometric mean**. Geometric because the requirements are conjunctive: a peptide that cannot
be delivered is worthless however well it fits its target, and an arithmetic mean lets a
strong term paper over a zero.

Shared developability terms: permeation (Potts–Guy log Kp), chemical stability (deamidation,
oxidation, diketopiperazine, acid-labile bonds), synthesisability, protease-site density,
and a cationic-amphipathicity irritation proxy.

### Delivery form is programme-specific

Permeation is scored on the form the peptide can actually be shipped in.

- `matrix_signal`, `tyrosinase_modulator`, `barrier_amp` — N-palmitoyl
- `snare_competitor` — acetyl / C-amide, as acetyl hexapeptide-8 is
- `cu_carrier` — **free N-terminus**, forced by the coordination chemistry

That last one is not a preference. Palmitoylation acylates the α-amine, and the α-amine is
the first nitrogen an ATCUN or GHK-type site uses to hold Cu(II). Scoring a copper carrier on
its lipidated permeability credits it for a molecule that has lost the mechanism it was
selected for. The same argument applies outside this repository: **palmitoyl tripeptide-1
cannot be a copper-delivery peptide**, whatever else pal-GHK does.

### Copper coordination, stated correctly

GHK is routinely called an ATCUN motif. It is not.

| Motif | Site | Donors | Score |
| --- | --- | --- | ---: |
| Xaa-Xaa-**His** | ATCUN | α-amine, two deprotonated backbone amides, imidazole (4N) | 1.00 |
| Xaa-**His** | GHK-type | α-amine, one deprotonated amide, imidazole (3N) | 0.90 |
| His elsewhere | — | none of the above | 0.15 |

Proline is handled by which nitrogen it removes. At position two of an ATCUN site it replaces
the backbone NH that must deprotonate, so the site is abolished (score 0). At position one it
leaves a secondary rather than primary α-amine — a discount (0.75), not a disqualification.

## 4. Nulls

Four control families, five decoys per design, length-matched design by design because
several terms are explicit length windows.

| Family | Holds fixed | What separation against it proves |
| --- | --- | --- |
| `uniprot_background` | length | almost nothing — the ECM is 15% Gly and 10% Pro |
| `ecm_background` | length, pooled composition | that the score is not pure ECM-composition |
| `natural_fragment` | length, real residue order | that the pipeline found something rather than re-describing the proteome |
| `scrambled` | length, composition **exactly** | that the score is about the sequence |

`composition_distance` is reported for every family so "composition-matched" is a measurement
rather than an assertion. The scrambled family must come out at exactly 0.000.

## 5. Permutation invariance

Of sixteen descriptors, ten are functions of amino-acid composition alone and take the same
value for a sequence and for any shuffle of it. `tests/test_invariance.py` asserts this
residue by residue rather than trusting the claim.

**Invariant:** molecular weight, charge at pH 5.5 and 7.4, isoelectric point, GRAVY, Boman
index, aliphatic index, length, and the two composition fractions.
**Order-reading:** hydrophobic moment (a vector sum around a helical wheel), instability index
(a sum over dipeptide weights), and the three motif/liability scanners.

Weight sitting on an invariant term cannot distinguish a design from its own scramble, even
in principle. That is measurable in advance, and it is the quantity that predicts what the
pipeline achieves.

### The correlation, and what it is conditional on

| Subset | *r* | n | Reproducible? |
| --- | ---: | ---: | --- |
| Encrypted track (enumeration only) | **0.957** | 5 | Exactly — no stochastic search |
| De novo, 30 generations x 250 | 0.977 | 5 | Seeded; converged |
| De novo, 6 generations x 60 | 0.688 | 5 | Seeded; **under-converged** |
| De novo and encrypted averaged per programme | 0.989 | 5 | — |
| Pooled, full budget | 0.929 | 10 rows, 5 distinct x | — |

Order-sensitive weight predicts an objective's **ceiling**, not a given campaign's result. At a
twentieth of the search budget the de novo correlation drops by nearly 0.3, because the
optimiser has not reached the ceiling in any programme and sits at a different distance below
it in each. The encrypted track is the load-bearing evidence: a ranked enumeration over all
40,026 fragments with no search in it, so it reproduces bit-for-bit from a fresh clone with a
cold cache, and gives *r* = 0.957 on its own.

Quoting the pooled 0.929 without saying which budget produced it would be the same error this
repository is about — a number whose value depends on a choice that was not reported.

The pooled row carries a second problem, and it is a denominator problem too. Order-sensitive
weight is a property of the objective, so a programme's two tracks share an identical x value:
the ten pooled points sit over five distinct weights. `scipy.stats.pearsonr` has no way to know
that and returns p = 0.0001 on the assumption of ten independent observations, against p = 0.0105
for the five encrypted points alone. The correlation coefficient is a fair descriptive summary;
the pooled p-value is not, and this repository does not quote it.
`cosmepep audit` prints all three.

## 6. Known non-results

Three terms were built, measured, and found to contribute nothing. They are documented rather
than deleted, because each is the kind of term that appears in published scoring functions
without its null.

- **Lipidation gain.** Palmitoylation adds 238.4 Da and, since Crippen log P is additive, a
  fixed +5.64 log P. Potts–Guy is linear in both, so the permeability gain is the same 2.55
  log units for *every* sequence. It looks like a delivery term and is arithmetically a
  constant.
- **Protease resistance.** Against the only ground truth available — the marketed peptides —
  it has no discriminating power. Their site densities run the full 0.0–1.0 range: GHK, GQPR
  and TKPR at 0.0, KVK, LSVD and YR at 1.0, and all six are sold. Kept at low weight as a
  tie-break.
- **Irritation margin.** A hydrophobic moment is a vector sum around a helical wheel; below
  about two turns there is no wheel. The term returns 1.0 for anything shorter than seven
  residues, which is four of the five programmes, and its measured AUC against scrambles is
  0.500 in each.

## 7. Two defects found by the nulls

Both were introduced by the author and caught by measurement rather than review.

- **A hand-picked composition term rejected a clinical compound.** "Collagen character" as the
  Gly/Pro/Ala/Lys fraction scores LSVD at exactly zero, and the geometric mean propagates
  that to a composite of 0.051. LSVD is one of two tetrapeptide matrikines taken through a
  split-face clinical study. The term was replaced with per-residue log2 enrichment in the
  panel relative to UniProt background — derived rather than guessed — and LSVD scores 0.72.
- **A helix mimic full of prolines.** The SNARE programme counted Chou–Fasman helix formers
  and did not subtract breakers, so the optimiser returned `EEPLKK`: acidic head, basic tail,
  and a proline in the middle of a sequence whose entire premise is that it forms a helix. Not
  counting a defect is not the same as penalising it.

## 8. What this is not

No structure, no docking, no folding, no learned model, and no wet-lab data. The output is a
prioritised, chemically annotated hypothesis list with its null distributions attached. Every
composite is a desirability aggregate over heuristics, and none of the five mechanisms has
been demonstrated for any sequence here.
