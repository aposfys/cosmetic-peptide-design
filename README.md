# cosmetic-peptide-design
How much of a peptide design result survives a control that holds composition fixed?

[![CI](https://github.com/aposfys/cosmetic-peptide-design/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/cosmetic-peptide-design/actions/workflows/ci.yml)
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
cosmepep sensitivity  # the same campaign at 6 x 60, for the budget caveat
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
run, about how much separation can survive. [METHODS](docs/METHODS.md) section 5 measures that
invariance instead of assuming it, and records the one place it is only approximate.

It is a good prediction. Per programme, order-sensitive weight against scrambled AUC gives
**r = 0.957** on the encrypted track (p = 0.011) and **0.977** on the de novo track (p = 0.004);
averaging the two tracks gives 0.989.

Pooling all ten cells gives r = 0.93, and that figure takes a caveat rather than a p-value.
Order weight is a property of the *programme*, so both tracks of a programme carry an identical
x and the pooled row is ten points over five distinct weights. The n = 5 rows are the ones with
an honest denominator, and [METHODS](docs/METHODS.md) section 5 sets out why.

Quote any of them with the search budget. Order weight predicts a *ceiling*, and a campaign
reaches its ceiling only once the search converges. Rerun at 6 generations × 60 population
instead of 30 × 250 and the de novo figure falls to **0.689**, the pooled one to **0.742**,
because an under-converged optimiser sits at a different distance below each programme's limit.
`cosmepep sensitivity` writes that rerun to
[results/budget_sensitivity.json](results/budget_sensitivity.json). The load-bearing figure is
the encrypted track's **r = 0.957**, a ranked enumeration with no stochastic search in it, which
reproduces bit-for-bit from a fresh clone with a cold cache, and it is the track the assertion
in `tests/` is pinned to.

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

Fifteen marketed cosmetic peptides sit in the pool of whichever programme matches their
mechanism class. Eleven fall in `matrix_signal`, including KTTKS, whose clinical evidence is
the reason the class exists, and GPKG and LSVD from a published split-face trial. `cu_carrier`
and `snare_competitor` carry one each, two belong to neither class, and `tyrosinase_modulator`
and `barrier_amp` carry none, so those two have no calibration at all.

Twenty-four of the twenty-six reference rows land at **percentile 0**, and that is not evidence
the designs are better. The designs came out of a search *maximising this composite* and the
marketed peptides did not, so the percentile restates the objective rather than validating it.
It stays in the output because a pipeline that quietly dropped it would look like it had a
calibration. Only two rows land *inside* their pool, and they are the only cells where the
comparison carries information at all: palmitoyl tetrapeptide-7 (GQPR) at percentile 13 in the
`matrix_signal` de novo pool, and acetyl hexapeptide-8 (EEMQRR) at percentile 34 in the
`snare_competitor` enumeration.

**Nothing here was measured in a laboratory.** Every quantity is a model output over
enumerated fragments and published amino-acid scales. The shortlists are hypotheses.

### Two corrections the coordination chemistry forced

**GHK is not an ATCUN motif.** ATCUN is Xaa-Xaa-His — four-nitrogen square-planar Cu(II)
coordination. GHK is Xaa-His: one intervening amide, three nitrogens, lower affinity. Both are
scored here, at 1.00 and 0.90. The copper-peptide category rests on the second and is routinely
described as the first.

**Palmitoyl tripeptide-1 cannot be a copper carrier.** Palmitoylation acylates the α-amine, and
the α-amine is the first nitrogen either site uses to hold Cu(II). Permeation is therefore scored
per programme on the form the peptide can actually ship in — free N-terminus for `cu_carrier`,
acetyl/C-amide for `snare_competitor`, palmitoyl for the rest.

### More

- [Results](results/RESULTS.md) covers every programme, track, control family and per-term AUC.
- [Methods](docs/METHODS.md) covers the design space, the protease rules, the four nulls,
  permutation invariance, plateaus, the three terms that measure nothing, the two defects the
  nulls caught, and the prior work the order-weight claim has to be read against.
- [results/budget_sensitivity.json](results/budget_sensitivity.json) holds the reduced-budget
  rerun behind the 0.689 and the 0.742.

Sequence data in `data/cache/ecm_proteins.json` comes from
[UniProt](https://www.uniprot.org) and is used under CC BY 4.0. The CLI resolves `data/cache/`
and `results/` relative to the installed package, so with more than one clone on a machine
`cosmepep audit` reads the campaign belonging to the clone that was installed, whatever the
working directory.
