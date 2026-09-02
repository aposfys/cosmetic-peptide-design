# Results

Seed 0, 30 generations, population 250, 120 candidates per track into the null comparison.

Every AUC below is the same designs scored by the same function. Only the control
changes. The spread between the first row and the last of each table is the size of
the choice a pipeline makes when it reports a single number.

## Procollagen-stimulating matrikine (`matrix_signal`)

**Target.** dermal fibroblast; procollagen I/III and fibronectin transcription  
**Reference compound.** KTTKS, GEKG, GPKG  
**Deliverable form.** palmitoyl — permeation is scored on this form, not on a lipidated one the mechanism could not survive.  
**Weight on order-sensitive terms.** 42% — the rest cannot distinguish a design from its own scramble even in principle.

### encrypted track

Plateau: 29 candidates tied at 1.0; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.999 | [0.995, 1.000] | 1.87 | 0.511 |
| ecm_background | 0.979 | [0.968, 0.988] | 1.64 | 0.373 |
| natural_fragment | 0.992 | [0.985, 0.997] | 1.73 | 0.368 |
| **scrambled** | **0.807** | [0.775, 0.839] | 1.56 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| KTTKS | Palmitoyl pentapeptide-4 | clinical | 0.782 | 0 |
| GQPR | Palmitoyl tetrapeptide-7 | in_vitro | 0.979 | 0 |
| KDVY | Acetyl tetrapeptide-5 | supplier | 0.854 | 0 |
| FVAPFP | Hexapeptide-11 | supplier | 0.356 | 0 |
| KVK | Palmitoyl tripeptide-5 | supplier | 0.734 | 0 |
| GEKG | Tetrapeptide (GEKG) | in_vitro | 0.962 | 0 |
| VGVAPG | Palmitoyl hexapeptide-12 | in_vitro | 0.369 | 0 |
| GPKG | GPKG tetrapeptide matrikine | clinical | 0.797 | 0 |
| LSVD | LSVD tetrapeptide matrikine | clinical | 0.724 | 0 |
| TKPR | Tuftsin (tetrapeptide-1) | in_vitro | 0.934 | 0 |
| VW | Dipeptide-2 | supplier | 0.281 | 0 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `GLPGPP` | encrypted | COL1A1 | 1.000 | 775 | -4.38 | +0.00 | — |
| `GGGY` | encrypted | COL1A2 | 1.000 | 591 | -3.94 | +0.00 | — |
| `AGPPG` | encrypted | COL3A1 | 1.000 | 636 | -4.26 | +0.00 | — |
| `GGAPP` | encrypted | COL18A1 | 1.000 | 636 | -4.26 | +0.00 | — |
| `GYGPG` | encrypted | ELN | 1.000 | 688 | -4.54 | +0.00 | — |
| `THPG` | encrypted | FN1 | 1.000 | 649 | -4.55 | +0.76 | — |
| `GKPG` | encrypted | COL1A1 | 1.000 | 596 | -3.93 | +1.00 | — |
| `RGPP` | encrypted | COL1A1 | 1.000 | 664 | -4.34 | +1.00 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (length_fit, matrix_composition, charge_fit, hydrophilicity, permeation) sit at 0.500 by construction; the separation such as it is comes from lipidation_handle, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| length_fit | 1.0 | no | 0.500 |
| matrix_composition | 0.9 | no | 0.500 |
| charge_fit | 0.7 | no | 0.500 |
| hydrophilicity | 0.5 | no | 0.500 |
| lipidation_handle | 0.8 | yes | 0.768 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.664 |
| synthesisability | 0.6 | yes | 0.503 |
| protease_resistance | 0.3 | yes | 0.545 |
| irritation_margin | 0.5 | yes | 0.500 |

### de_novo track

Plateau: 39 candidates tied at 1.0; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.998 | [0.996, 1.000] | 1.85 | 0.644 |
| ecm_background | 0.987 | [0.980, 0.994] | 1.69 | 0.524 |
| natural_fragment | 0.993 | [0.987, 0.997] | 1.74 | 0.511 |
| **scrambled** | **0.790** | [0.756, 0.823] | 1.45 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| KTTKS | Palmitoyl pentapeptide-4 | clinical | 0.782 | 0 |
| GQPR | Palmitoyl tetrapeptide-7 | in_vitro | 0.979 | 13 |
| KDVY | Acetyl tetrapeptide-5 | supplier | 0.854 | 0 |
| FVAPFP | Hexapeptide-11 | supplier | 0.356 | 0 |
| KVK | Palmitoyl tripeptide-5 | supplier | 0.734 | 0 |
| GEKG | Tetrapeptide (GEKG) | in_vitro | 0.962 | 0 |
| VGVAPG | Palmitoyl hexapeptide-12 | in_vitro | 0.369 | 0 |
| GPKG | GPKG tetrapeptide matrikine | clinical | 0.797 | 0 |
| LSVD | LSVD tetrapeptide matrikine | clinical | 0.724 | 0 |
| TKPR | Tuftsin (tetrapeptide-1) | in_vitro | 0.934 | 0 |
| VW | Dipeptide-2 | supplier | 0.281 | 0 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `HHPP` | de_novo | — | 1.000 | 725 | -3.97 | +1.52 | — |
| `GGGPY` | de_novo | — | 1.000 | 688 | -4.54 | +0.00 | — |
| `GFSPP` | de_novo | — | 1.000 | 742 | -4.49 | +0.00 | — |
| `HGAPP` | de_novo | — | 1.000 | 716 | -4.51 | +0.76 | — |
| `HTGP` | de_novo | — | 1.000 | 649 | -4.55 | +0.76 | — |
| `GAHPP` | de_novo | — | 1.000 | 716 | -4.51 | +0.76 | — |
| `GGHP` | de_novo | — | 1.000 | 605 | -4.10 | +0.76 | — |
| `FGSPP` | de_novo | — | 1.000 | 742 | -4.49 | +0.00 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (length_fit, matrix_composition, charge_fit, hydrophilicity, permeation) sit at 0.500 by construction; the separation such as it is comes from lipidation_handle, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| length_fit | 1.0 | no | 0.500 |
| matrix_composition | 0.9 | no | 0.500 |
| charge_fit | 0.7 | no | 0.500 |
| hydrophilicity | 0.5 | no | 0.500 |
| lipidation_handle | 0.8 | yes | 0.749 |
| permeation | 1.0 | no | 0.501 |
| chemical_stability | 0.8 | yes | 0.657 |
| synthesisability | 0.6 | yes | 0.491 |
| protease_resistance | 0.3 | yes | 0.523 |
| irritation_margin | 0.5 | yes | 0.500 |

## Copper-carrier tripeptide (`cu_carrier`)

**Target.** Cu(II) delivery via an N-terminal ATCUN or GHK-type chelation site  
**Reference compound.** GHK / GHK-Cu  
**Deliverable form.** free — permeation is scored on this form, not on a lipidated one the mechanism could not survive.  
**Weight on order-sensitive terms.** 55% — the rest cannot distinguish a design from its own scramble even in principle.

### encrypted track

Plateau: 40 candidates tied at 1.0; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.999 | [0.996, 1.000] | 5.98 | 0.345 |
| ecm_background | 0.997 | [0.994, 1.000] | 6.02 | 0.402 |
| natural_fragment | 0.995 | [0.990, 0.999] | 7.24 | 0.407 |
| **scrambled** | **0.906** | [0.885, 0.925] | 1.53 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| GHK | Tripeptide-1 / copper tripeptide-1 | in_vitro | 0.973 | 0 |

| Sequence | Provenance | Source | Composite | MW (free) | log Kp (free) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `GLH` | encrypted | COL1A2 | 1.000 | 325 | -5.39 | +0.76 | — |
| `ISHVI` | encrypted | COL4A3 | 1.000 | 568 | -6.84 | +0.76 | — |
| `KAHPV` | encrypted | LAMA3 | 1.000 | 551 | -7.04 | +1.76 | — |
| `LLHPT` | encrypted | COL3A1 | 1.000 | 580 | -6.75 | +0.76 | — |
| `ANHF` | encrypted | COL6A3 | 1.000 | 488 | -7.12 | +0.76 | — |
| `SYHL` | encrypted | COL7A1 | 1.000 | 519 | -6.71 | +0.76 | — |
| `GIHGL` | encrypted | COL4A3 | 1.000 | 496 | -6.68 | +0.76 | — |
| `TGHY` | encrypted | LAMA3 | 1.000 | 476 | -7.18 | +0.76 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (no_competing_thiol, length_fit, charge_fit, permeation) sit at 0.500 by construction; the separation such as it is comes from copper_motif, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| copper_motif | 2.0 | yes | 0.893 |
| no_competing_thiol | 1.0 | no | 0.500 |
| length_fit | 0.8 | no | 0.500 |
| charge_fit | 0.6 | no | 0.500 |
| permeation | 1.0 | no | 0.501 |
| chemical_stability | 0.8 | yes | 0.518 |
| synthesisability | 0.6 | yes | 0.521 |
| protease_resistance | 0.3 | yes | 0.606 |
| irritation_margin | 0.5 | yes | 0.500 |

### de_novo track

Plateau: 165 candidates tied at 1.0; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 1.000 | [1.000, 1.000] | 6.32 | 0.494 |
| ecm_background | 0.999 | [0.998, 1.000] | 5.61 | 0.453 |
| natural_fragment | 1.000 | [1.000, 1.000] | 8.05 | 0.459 |
| **scrambled** | **0.912** | [0.897, 0.927] | 1.49 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| GHK | Tripeptide-1 / copper tripeptide-1 | in_vitro | 0.973 | 0 |

| Sequence | Provenance | Source | Composite | MW (free) | log Kp (free) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `GFHGI` | de_novo | — | 1.000 | 530 | -6.74 | +0.76 | — |
| `HAH` | de_novo | — | 1.000 | 363 | -6.11 | +1.52 | — |
| `KYHYP` | de_novo | — | 1.000 | 707 | -7.13 | +1.76 | — |
| `AIHHL` | de_novo | — | 1.000 | 590 | -6.74 | +1.52 | — |
| `ATHV` | de_novo | — | 1.000 | 426 | -6.81 | +0.76 | — |
| `FKHGP` | de_novo | — | 1.000 | 585 | -7.11 | +1.76 | — |
| `TYHA` | de_novo | — | 1.000 | 491 | -6.99 | +0.76 | — |
| `KSHF` | de_novo | — | 1.000 | 518 | -7.15 | +1.76 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (no_competing_thiol, length_fit, charge_fit, permeation) sit at 0.500 by construction; the separation such as it is comes from copper_motif, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| copper_motif | 2.0 | yes | 0.870 |
| no_competing_thiol | 1.0 | no | 0.500 |
| length_fit | 0.8 | no | 0.500 |
| charge_fit | 0.6 | no | 0.500 |
| permeation | 1.0 | no | 0.501 |
| chemical_stability | 0.8 | yes | 0.522 |
| synthesisability | 0.6 | yes | 0.542 |
| protease_resistance | 0.3 | yes | 0.678 |
| irritation_margin | 0.5 | yes | 0.500 |

## SNARE-competing relaxant (`snare_competitor`)

**Target.** SNAP-25 / SNARE assembly at the neuromuscular junction  
**Reference compound.** EEMQRR (acetyl hexapeptide-8)  
**Deliverable form.** acetyl_amide — permeation is scored on this form, not on a lipidated one the mechanism could not survive.  
**Weight on order-sensitive terms.** 60% — the rest cannot distinguish a design from its own scramble even in principle.

### encrypted track

Plateau: 1 candidates tied at 0.8878; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.997 | [0.993, 1.000] | 6.63 | 0.349 |
| ecm_background | 0.999 | [0.995, 1.000] | 9.09 | 0.440 |
| natural_fragment | 0.999 | [0.997, 1.000] | 8.12 | 0.429 |
| **scrambled** | **0.939** | [0.921, 0.954] | 2.22 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| EEMQRR | Acetyl hexapeptide-8 (Argireline) | clinical | 0.759 | 34 |

| Sequence | Provenance | Source | Composite | MW (acetyl_amide) | log Kp (acetyl_amide) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `EEIKR` | encrypted | LAMA3 | 0.888 | 715 | -9.13 | +0.11 | — |
| `DEGALK` | encrypted | COL6A3 | 0.865 | 673 | -9.24 | -0.93 | — |
| `DELLIR` | encrypted | HSPG2 | 0.861 | 799 | -8.89 | -0.93 | beta_sheet_run:LLI |
| `EEVRK` | encrypted | COL17A1 | 0.882 | 701 | -9.32 | +0.11 | — |
| `DDVKR` | encrypted | COL17A1 | 0.870 | 673 | -9.71 | +0.03 | — |
| `EEPLMR` | encrypted | COL6A3 | 0.853 | 815 | -9.31 | -0.89 | oxidation_M:M |
| `EDVRR` | encrypted | LAMB3 | 0.856 | 715 | -10.30 | +0.07 | polyarginine:RR |
| `DEFFLK` | encrypted | COL6A3 | 0.854 | 839 | -8.24 | -0.93 | beta_sheet_run:FFL, bulky_pair:FF |

Per-term AUC against the scrambled control. Terms declared composition-only (length_fit, near_neutral, helix_propensity, permeation) sit at 0.500 by construction; the separation such as it is comes from acidic_n_terminus, basic_c_terminus, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| acidic_n_terminus | 1.6 | yes | 0.821 |
| basic_c_terminus | 0.9 | yes | 0.808 |
| length_fit | 0.8 | no | 0.500 |
| near_neutral | 0.6 | no | 0.498 |
| helix_propensity | 0.7 | no | 0.500 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.509 |
| synthesisability | 0.6 | yes | 0.481 |
| protease_resistance | 0.3 | yes | 0.565 |
| irritation_margin | 0.5 | yes | 0.500 |

### de_novo track

Plateau: 1 candidates tied at 0.9623; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 1.000 | [1.000, 1.000] | 7.17 | 0.524 |
| ecm_background | 1.000 | [1.000, 1.000] | 9.71 | 0.600 |
| natural_fragment | 1.000 | [1.000, 1.000] | 10.19 | 0.592 |
| **scrambled** | **0.987** | [0.979, 0.993] | 2.12 | 0.000 |

Marketed peptides of this class, placed in the same pool. A percentile of 0 is **not** a finding: the pool is the top 120 candidates by this composite and the reference compounds were not selected on it, so the comparison restates the objective rather than validating it. The informative case is a reference landing *inside* the pool, which happens once -- acetyl hexapeptide-8 at the 34th percentile of the SNARE programme.

| Peptide | INCI | Evidence | Composite | Percentile in pool |
| --- | --- | --- | ---: | ---: |
| EEMQRR | Acetyl hexapeptide-8 (Argireline) | clinical | 0.759 | 0 |

| Sequence | Provenance | Source | Composite | MW (acetyl_amide) | log Kp (acetyl_amide) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `EEPLKK` | de_novo | — | 0.962 | 784 | -8.94 | +0.11 | — |
| `DDLLKK` | de_novo | — | 0.959 | 772 | -9.04 | +0.03 | — |
| `DEAIKK` | de_novo | — | 0.937 | 744 | -9.32 | +0.07 | — |
| `EDIVKK` | de_novo | — | 0.931 | 772 | -9.04 | +0.07 | — |
| `DEILKK` | de_novo | — | 0.953 | 786 | -8.85 | +0.07 | — |
| `EDYLKK` | de_novo | — | 0.941 | 836 | -9.22 | +0.07 | — |
| `EELIKK` | de_novo | — | 0.940 | 800 | -8.66 | +0.11 | — |
| `DEFAKK` | de_novo | — | 0.935 | 778 | -9.39 | +0.07 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (length_fit, near_neutral, helix_propensity, permeation) sit at 0.500 by construction; the separation such as it is comes from acidic_n_terminus, basic_c_terminus, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| acidic_n_terminus | 1.6 | yes | 0.958 |
| basic_c_terminus | 0.9 | yes | 0.942 |
| length_fit | 0.8 | no | 0.500 |
| near_neutral | 0.6 | no | 0.500 |
| helix_propensity | 0.7 | no | 0.500 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.526 |
| synthesisability | 0.6 | yes | 0.489 |
| protease_resistance | 0.3 | yes | 0.576 |
| irritation_margin | 0.5 | yes | 0.500 |

## Tyrosinase-modulating brightener (`tyrosinase_modulator`)

**Target.** tyrosinase type-3 dicopper active site  
**Reference compound.** no well-attested peptide reference; the weakest-grounded programme here  
**Deliverable form.** palmitoyl — permeation is scored on this form, not on a lipidated one the mechanism could not survive.  
**Weight on order-sensitive terms.** 31% — the rest cannot distinguish a design from its own scramble even in principle.

### encrypted track

Plateau: 3 candidates tied at 0.9226; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.996 | [0.992, 0.999] | 11.03 | 0.498 |
| ecm_background | 1.000 | [1.000, 1.000] | 15.45 | 0.516 |
| natural_fragment | 0.998 | [0.995, 1.000] | 10.41 | 0.514 |
| **scrambled** | **0.561** | [0.509, 0.614] | 0.22 | 0.000 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `HGHY` | encrypted | NID1 | 0.923 | 751 | -4.44 | +1.52 | — |
| `HLYPH` | encrypted | FN1 | 0.888 | 904 | -4.38 | +1.52 | — |
| `HFHTV` | encrypted | HSPG2 | 0.878 | 878 | -5.08 | +1.52 | — |
| `HFHL` | encrypted | HSPG2 | 0.923 | 791 | -3.47 | +1.52 | — |
| `HGHF` | encrypted | FLG | 0.923 | 735 | -4.13 | +1.52 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (copper_chelator, substrate_mimicry, length_fit, cationic, permeation) sit at 0.500 by construction; the separation such as it is comes from chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| copper_chelator | 1.4 | no | 0.500 |
| substrate_mimicry | 1.2 | no | 0.500 |
| length_fit | 0.7 | no | 0.500 |
| cationic | 0.5 | no | 0.500 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.511 |
| synthesisability | 0.6 | yes | 0.519 |
| protease_resistance | 0.3 | yes | 0.599 |
| irritation_margin | 0.5 | yes | 0.500 |

### de_novo track

Plateau: 69 candidates tied at 1.0; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 1.000 | [1.000, 1.000] | 11.00 | 0.760 |
| ecm_background | 1.000 | [1.000, 1.000] | 15.29 | 0.786 |
| natural_fragment | 1.000 | [1.000, 1.000] | 14.17 | 0.784 |
| **scrambled** | **0.701** | [0.650, 0.750] | 0.69 | 0.000 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `FHFHG` | de_novo | — | 1.000 | 882 | -4.51 | +1.52 | — |
| `AYHHF` | de_novo | — | 1.000 | 912 | -4.63 | +1.52 | — |
| `FHHVY` | de_novo | — | 1.000 | 940 | -4.35 | +1.52 | — |
| `HYHY` | de_novo | — | 1.000 | 857 | -4.15 | +1.52 | — |
| `YHPFH` | de_novo | — | 1.000 | 938 | -4.44 | +1.52 | — |
| `HFHFL` | de_novo | — | 1.000 | 938 | -3.85 | +1.52 | — |
| `FHHF` | de_novo | — | 1.000 | 825 | -3.54 | +1.52 | — |
| `IHFHY` | de_novo | — | 1.000 | 954 | -4.16 | +1.52 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (copper_chelator, substrate_mimicry, length_fit, cationic, permeation) sit at 0.500 by construction; the separation such as it is comes from chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| copper_chelator | 1.4 | no | 0.500 |
| substrate_mimicry | 1.2 | no | 0.500 |
| length_fit | 0.7 | no | 0.500 |
| cationic | 0.5 | no | 0.500 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.512 |
| synthesisability | 0.6 | yes | 0.684 |
| protease_resistance | 0.3 | yes | 0.624 |
| irritation_margin | 0.5 | yes | 0.500 |

## Selective anti-C. acnes peptide (`barrier_amp`)

**Target.** Cutibacterium acnes membrane, at skin surface pH  
**Reference compound.** host defence peptide chemotype  
**Deliverable form.** palmitoyl — permeation is scored on this form, not on a lipidated one the mechanism could not survive.  
**Weight on order-sensitive terms.** 47% — the rest cannot distinguish a design from its own scramble even in principle.

### encrypted track

Plateau: 1 candidates tied at 0.8941; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.991 | [0.984, 0.996] | 4.80 | 0.292 |
| ecm_background | 0.997 | [0.992, 1.000] | 5.56 | 0.397 |
| natural_fragment | 0.998 | [0.995, 1.000] | 6.00 | 0.400 |
| **scrambled** | **0.784** | [0.751, 0.819] | 0.94 | 0.000 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `LKQLKKL` | encrypted | LUM | 0.894 | 1109 | -6.10 | +3.00 | — |
| `KIGKLQK` | encrypted | LAMA3 | 0.872 | 1052 | -6.76 | +3.00 | — |
| `VRKMKPL` | encrypted | COL6A3 | 0.883 | 1110 | -6.32 | +3.00 | oxidation_M:M |
| `ITKVRKV` | encrypted | DCN | 0.871 | 1082 | -7.02 | +3.00 | — |
| `KKTTKFL` | encrypted | COL5A1 | 0.880 | 1104 | -7.03 | +3.00 | — |
| `KKALKLM` | encrypted | COL3A1 | 0.876 | 1070 | -5.53 | +3.00 | oxidation_M:M |
| `KGLKGLP` | encrypted | COL4A3 | 0.850 | 950 | -5.61 | +2.00 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (cationicity, selectivity, length_fit, permeation) sit at 0.500 by construction; the separation such as it is comes from amphipathicity, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| cationicity | 1.3 | no | 0.500 |
| amphipathicity | 1.3 | yes | 0.767 |
| selectivity | 1.0 | no | 0.502 |
| length_fit | 0.6 | no | 0.500 |
| permeation | 1.0 | no | 0.500 |
| chemical_stability | 0.8 | yes | 0.527 |
| synthesisability | 0.6 | yes | 0.531 |
| protease_resistance | 0.3 | yes | 0.548 |
| irritation_margin | 0.5 | yes | 0.268 |

### de_novo track

Plateau: 1 candidates tied at 0.9491; the shortlist is a max-min diversity selection from that plateau, not a ranking.

| Control family | AUC | 95% CI | Cohen's *d* | Composition distance |
| --- | ---: | --- | ---: | ---: |
| uniprot_background | 0.999 | [0.996, 1.000] | 5.10 | 0.533 |
| ecm_background | 1.000 | [1.000, 1.000] | 6.80 | 0.589 |
| natural_fragment | 1.000 | [1.000, 1.000] | 6.65 | 0.586 |
| **scrambled** | **0.895** | [0.872, 0.916] | 0.91 | 0.000 |

| Sequence | Provenance | Source | Composite | MW (palmitoyl) | log Kp (palmitoyl) | q(5.5) | Liabilities |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| `KKIKPLL` | de_novo | — | 0.949 | 1078 | -5.03 | +3.00 | — |
| `LRFKKPI` | de_novo | — | 0.927 | 1140 | -5.88 | +3.00 | — |
| `IIRIKKP` | de_novo | — | 0.919 | 1106 | -5.81 | +3.00 | — |
| `KKLKPVF` | de_novo | — | 0.945 | 1098 | -5.29 | +3.00 | — |
| `KKVKPPI` | de_novo | — | 0.932 | 1047 | -5.51 | +3.00 | — |
| `KKIHPII` | de_novo | — | 0.926 | 1087 | -5.20 | +2.76 | — |
| `LKIKKLI` | de_novo | — | 0.923 | 1094 | -4.74 | +3.00 | — |
| `VRIKKPF` | de_novo | — | 0.921 | 1126 | -6.07 | +3.00 | — |

Per-term AUC against the scrambled control. Terms declared composition-only (cationicity, selectivity, length_fit, permeation) sit at 0.500 by construction; the separation such as it is comes from amphipathicity, chemical_stability, synthesisability, protease_resistance, irritation_margin.

| Term | Weight | Order-sensitive | AUC vs scrambled |
| --- | ---: | --- | ---: |
| cationicity | 1.3 | no | 0.500 |
| amphipathicity | 1.3 | yes | 0.726 |
| selectivity | 1.0 | no | 0.500 |
| length_fit | 0.6 | no | 0.500 |
| permeation | 1.0 | no | 0.506 |
| chemical_stability | 0.8 | yes | 0.548 |
| synthesisability | 0.6 | yes | 0.604 |
| protease_resistance | 0.3 | yes | 0.517 |
| irritation_margin | 0.5 | yes | 0.344 |
