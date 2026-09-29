# cosmetic-peptide-design
How much of a peptide design result survives a control that holds composition fixed?

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Five design programmes against the skin ECM degradome — procollagen matrikines, Cu(II)
carriers, SNARE competitors, tyrosinase modulators and a barrier AMP — each run down two
tracks and scored against four control families. The designs never change between controls.
Only the denominator does.

```
pip install -e ".[dev]"

cosmepep panel        # the skin ECM proteins the design space is drawn from
cosmepep degradome    # protease-liberated fragments
cosmepep design       # the full campaign (seed 0, 30 generations, population 250)
cosmepep evaluate     # separation table from an existing campaign
cosmepep audit        # order-sensitive weight against achieved separation
cosmepep report       # write results/RESULTS.md
make test
```

### The control decides the result

Every AUC below scores the same designs with the same function.

| Programme | Track | vs UniProt background | vs scrambled |
| --- | --- | ---: | ---: |
| matrix_signal | encrypted | 0.999 | **0.807** |
| matrix_signal | de_novo | 0.998 | **0.790** |
| cu_carrier | encrypted | 0.999 | **0.906** |
| cu_carrier | de_novo | 1.000 | **0.912** |
| snare_competitor | encrypted | 0.997 | **0.939** |
| snare_competitor | de_novo | 1.000 | **0.987** |
| tyrosinase_modulator | encrypted | 0.996 | **0.561** |
| tyrosinase_modulator | de_novo | 1.000 | **0.701** |
| barrier_amp | encrypted | 0.991 | **0.784** |
| barrier_amp | de_novo | 0.999 | **0.895** |

Against a background of UniProt sequences every programme looks near-perfect — 0.991 to
1.000, ten times out of ten. Against scrambles of the designs themselves, which hold amino
acid composition exactly fixed (measured distance 0.000) and differ only in residue order,
the same designs score 0.561 to 0.987. **The tyrosinase programme cannot reliably tell a
design from its own anagram**: AUC 0.561, and at a threshold rejecting 95% of scrambles it
keeps 7% of the designs.

### The size of that gap is predictable, and that is the useful result

Each programme's objective is a weighted sum of terms, and some of those terms — length,
composition, net charge, hydrophilicity, permeation — are functions of composition alone.
They return the same value for a sequence and its anagram, so they cannot contribute to
separation from a scrambled control **even in principle**. The fraction of objective weight
carried by the order-sensitive terms is therefore a prediction, made before any control is
run, about how much separation can survive.

It is a good prediction. Per programme, order-sensitive weight against scrambled AUC gives
**r = 0.957** on the encrypted track (p = 0.011) and **0.977** on the de novo track (p = 0.004);
averaging the two tracks gives 0.989.

Pooling all ten cells gives r = 0.93, but that figure needs a caveat rather than a p-value:
order weight is a property of the *programme*, so both tracks of a programme carry an identical
x. The pooled row is ten points over five distinct weights, and `pearsonr` treating it as ten
independent observations returns p = 0.0001 — roughly two orders of magnitude too small. The
n = 5 rows are the ones with an honest denominator.

Quote that number with its search budget. It predicts a *ceiling*, and a campaign reaches its
ceiling only once the search converges: rerun the de novo track at 6 generations × 60
population instead of 30 × 250 and the pooled figure falls to **0.740**, because an
under-converged optimiser sits at a different distance below each programme's limit. The
load-bearing figure is the encrypted track's **r = 0.957** — a ranked enumeration with no
stochastic search in it, which reproduces bit-for-bit from a fresh clone with a cold cache.
`cosmepep audit` prints all four, and the assertion in `tests/` is pinned to the
deterministic track for exactly this reason.

| Programme | Order-sensitive weight | Scrambled AUC |
| --- | ---: | ---: |
| tyrosinase_modulator | 0.31 | 0.561 / 0.701 |
| matrix_signal | 0.42 | 0.807 / 0.790 |
| barrier_amp | 0.47 | 0.784 / 0.895 |
| cu_carrier | 0.55 | 0.906 / 0.912 |
| snare_competitor | 0.60 | 0.939 / 0.987 |

An objective that is mostly composition terms will produce designs that a composition
control cannot distinguish, and the arithmetic says so in advance. `cosmepep audit` computes
this weight per programme, so a scrambled AUC can be predicted before the campaign runs and
compared against what it actually achieves.

### The reference calibration is circular, and is reported as such

Fifteen marketed cosmetic peptides are in the panel, each placed in the pool of the programme
whose mechanism class it belongs to. Eleven of them fall in `matrix_signal`, including KTTKS,
whose clinical evidence is the reason the class exists, and GPKG and LSVD from a published
split-face trial. `cu_carrier` and `snare_competitor` carry one each, two are in neither
class, and `tyrosinase_modulator` and `barrier_amp` carry none, so those two programmes have
no reference calibration at all. Of the twenty-six reference rows in the campaign,
twenty-four land at **percentile 0**: no design scores below them.

That is not evidence the designs are better. The designs were produced by an evolutionary
search *maximising this composite*, and the marketed peptides were not. A pool optimised on
a function will outrank anything not optimised on it, so the percentile restates the search
objective rather than validating it. It is kept in the output because a pipeline that
quietly dropped it would look like it had a calibration. Two rows land *inside* their pool
rather than under it, and those two are the only cells where the comparison carries any
information at all: palmitoyl tetrapeptide-7 (GQPR) at percentile 13 in the `matrix_signal`
de novo pool, and acetyl hexapeptide-8 (EEMQRR) at percentile 34 in the `snare_competitor`
enumeration.

**Nothing here was measured in a laboratory.** Every quantity is a model output over
enumerated fragments and published amino-acid scales. The shortlists are hypotheses.

### Plateaus, not rankings

Several tracks saturate: 29, 39, 40, 69 and 165 candidates tie at the top score. Where that
happens the shortlist is a max-min diversity selection across the plateau, not a ranking,
and the reports say so. A ranked table over tied scores would invent an order the objective
does not support.

### Two corrections the coordination chemistry forced

**GHK is not an ATCUN motif.** ATCUN is Xaa-Xaa-His — four-nitrogen square-planar Cu(II)
coordination. GHK is Xaa-His: one intervening amide, three nitrogens, lower affinity. Both are
scored here, at 1.00 and 0.90. The copper-peptide category rests on the second and is routinely
described as the first.

**Palmitoyl tripeptide-1 cannot be a copper carrier.** Palmitoylation acylates the α-amine, and
the α-amine is the first nitrogen either site uses to hold Cu(II). Permeation is therefore scored
per programme on the form the peptide can actually ship in — free N-terminus for `cu_carrier`,
acetyl/C-amide for `snare_competitor`, palmitoyl for the rest.

### Three terms that measure nothing, kept rather than deleted

The permeability gain from lipidation is arithmetically the same 2.55 log units for every
sequence. Protease-site density spans the full 0.0–1.0 range across marketed peptides, so it
separates nothing. The irritation proxy is a hydrophobic moment, undefined below two helical
turns and therefore inert in four of five programmes. Each is the kind of term that appears in
published scoring functions without its null — see [docs/METHODS.md](docs/METHODS.md) §6, along
with the two defects the nulls caught in this pipeline's own objectives.

### Where this sits in the literature

That the choice of control or negative set drives reported performance is established, and
this repository does not claim it as new. Virtual screening has a decoy-selection literature;
Guo et al. showed that even scaffold splits — long treated as the realistic choice —
overestimate performance, and Fooladi et al. found classical models and graph networks
largely indistinguishable under them. My [`peptide-binder-design`](https://github.com/aposfys/peptide-binder-design)
makes the same point on a published-style filter rather than a hand-built objective.

**The closest precedent for the predictor is AVE bias.** Wallach & Heifets (2017) defined a
computable measure of train–validation redundancy and showed it *strongly correlates with the
performance of ligand-based methods* — the same move made here: a quantity you can compute
about your setup, in advance, whose magnitude tells you how much of your result is going to
be artifact. Anyone reading the order-sensitive weight as a new kind of claim should read it
as AVE bias applied to a different object.

The object is the difference, and it is a narrow one. AVE bias is a property of a **dataset
split** and needs the data to compute. Order-sensitive weight is a property of the
**objective's algebra** and needs nothing at all: a term that is a function of composition
alone returns an identical value for a sequence and its anagram, so weight sitting on it
cannot separate a design from a scrambled control *in principle*. Summing those weights is
arithmetic over the scoring function — no molecules, no split, no training run — and here it
predicts the achieved scrambled AUC at r = 0.957 on the enumerated track. Whether that
specific analogue has been published for hand-specified objectives, I do not know; the
searches I ran are below and came back empty, which is weak evidence and is not the same as
a literature review.

| Query | Verdict |
| --- | --- |
| `permutation invariant descriptors composition-matched controls peptide scoring` | Drifted to MS/MS peptide-spectrum scoring; nothing on point |
| `AVE bias metric predicts inflated benchmark performance` | Returns the precedent above, on dataset splits, not objectives |

Quote it with its search budget. It predicts a ceiling, and a campaign only reaches its
ceiling once the search converges; the caveat is documented above and is the reason the
enumerated track carries the claim.

### More

- [Results](results/RESULTS.md) — every programme, track, control family and per-term AUC
