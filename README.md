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

It is a good prediction: across the ten programme-track cells, order-sensitive weight
against scrambled AUC gives **r = 0.93 (p = 0.0001, R² = 0.86)**, slope 1.11.

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

Eleven marketed cosmetic peptides — including KTTKS, whose clinical evidence is the reason
the class exists, and GPKG and LSVD from a published split-face trial — are scored in the
same pool as the designs. Almost all of them land at **percentile 0**: no design scores
below them.

That is not evidence the designs are better. The designs were produced by an evolutionary
search *maximising this composite*, and the marketed peptides were not. A pool optimised on
a function will outrank anything not optimised on it, so the percentile restates the search
objective rather than validating it. It is kept in the output because a pipeline that
quietly dropped it would look like it had a calibration; the one exception, where a
reference reaches the 34th percentile, is the only cell where the comparison carries any
information at all.

**Nothing here was measured in a laboratory.** Every quantity is a model output over
enumerated fragments and published amino-acid scales. The shortlists are hypotheses.

### Plateaus, not rankings

Several tracks saturate: 29, 39, 40, 69 and 165 candidates tie at the top score. Where that
happens the shortlist is a max-min diversity selection across the plateau, not a ranking,
and the reports say so. A ranked table over tied scores would invent an order the objective
does not support.

### More

- [Results](results/RESULTS.md) — every programme, track, control family and per-term AUC
