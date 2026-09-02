# cosmetic-peptide-design

The same designs score **0.999 or 0.561** depending only on which null they are compared against.

De novo cosmetic peptide design from the human skin ECM degradome — five mechanism
programmes, RDKit + peptide descriptors, and four control families that decide what the
scores are worth.

```
pip install -e .
cosmepep panel        # the 19 skin matrix proteins the design space is drawn from
cosmepep degradome    # 40,026 fragments the skin's own proteases could liberate
cosmepep programmes   # the five objectives, term by term, flagged for order-sensitivity
cosmepep design       # run the campaign -> results/campaign.json + RESULTS.md
cosmepep audit        # the finding below
pytest                # 37 tests, no network after the first fetch
```

### An objective's ceiling is knowable before you run it

Every scoring term is flagged for whether it can see sequence order at all. Charge, GRAVY,
molecular weight, isoelectric point and the Boman index cannot — they are functions of
amino-acid composition and take the same value for a sequence and for any shuffle of it.
Weight sitting on those terms cannot separate a design from its own scramble even in
principle.

| Programme | Order-sensitive weight | AUC vs scrambled | AUC vs naive null |
| --- | ---: | ---: | ---: |
| `snare_competitor` | 60% | 0.939 | 0.997 |
| `cu_carrier` | 55% | 0.906 | 0.999 |
| `barrier_amp` | 47% | 0.784 | 0.991 |
| `matrix_signal` | 42% | 0.807 | 0.999 |
| **`tyrosinase_modulator`** | **31%** | **0.561** | **0.996** |

**Pearson *r* = 0.929** between the two middle columns across ten programme–track runs.

The last row is the point. Against a UniProt-random null the tyrosinase programme looks
indistinguishable from the other four — 0.996, a triumph. Against composition-preserving
scrambles its 95% interval is [0.509, 0.614] and it is barely above chance. Nothing about the
designs changed. Only the denominator did, and the denominator was predictable from the term
flags before a single sequence was scored.

### What that means for this field

Cosmetic peptide scoring is built almost entirely from composition descriptors: net charge,
GRAVY, molecular weight, hydrophobicity, pI. A pipeline reporting a pass rate against random
sequence has reported a fact about amino-acid frequencies. The fix is not a better descriptor,
it is a scramble control and a stated weight budget.

### Designs

Two provenance tracks per programme. `encrypted` peptides occur verbatim in human matrix
proteins — they have a biological prior and a freedom-to-operate problem. `de_novo` sequences
are forbidden from being panel substrings, so they have neither.

| Programme | Track | Shortlist |
| --- | --- | --- |
| Procollagen matrikine | de novo | `GGGPY` `GFSPP` `HGAPP` `HTGP` `GAHPP` `GGHP` `FGSPP` `HHPP` |
| Copper carrier | de novo | `GFHGI` `HAH` `KYHYP` `AIHHL` `ATHV` `FKHGP` `TYHA` `KSHF` |
| SNARE competitor | de novo | `DDLLKK` `DEAIKK` `EDIVKK` `DEILKK` `EDYLKK` `EELIKK` `DEFAKK` |
| Tyrosinase modulator | de novo | `FHFHG` `AYHHF` `FHHVY` `HYHY` `YHPFH` `HFHFL` `FHHF` `IHFHY` |
| Anti-*C. acnes* | de novo | `LRFKKPI` `IIRIKKP` `KKLKPVF` `KKVKPPI` `KKIHPII` `LKIKKLI` |

Each carries a full card in `results/campaign.json`: SMILES, MW and log P in the deliverable
form, Potts–Guy log Kp, charge at skin-surface and dermal pH, per-residue matrix enrichment,
itemised chemical and synthesis liabilities, and edit distance to the nearest marketed
peptide. **Shortlists are diversity selections, not rankings** — where the objective ties 165
candidates at 1.000, presenting the first eight would be a presentation artefact.

### Two corrections the nulls forced

**Palmitoyl tripeptide-1 is not a copper-delivery peptide.** Palmitoylation acylates the
α-amine, and the α-amine is the first nitrogen an ATCUN or GHK-type site uses to hold Cu(II).
A copper carrier scored on its lipidated permeability is credited for a molecule that lost the
mechanism it was selected for. Permeation is therefore scored per programme on the form the
peptide can actually ship in — free N-terminus for `cu_carrier`, acetyl/C-amide for
`snare_competitor`, palmitoyl for the rest.

**GHK is not an ATCUN motif.** ATCUN is Xaa-Xaa-His and gives four-nitrogen square-planar
coordination. GHK is Xaa-His, one intervening amide, three nitrogens. Both are scored, at
1.00 and 0.90; the whole copper-peptide category rests on the second.

### Known non-results, kept rather than deleted

Three terms were built, measured and found to contribute nothing — each the kind that appears
in published scoring functions without its null. The permeability gain from lipidation is
arithmetically the same 2.55 log units for every sequence. Protease-site density runs the full
0.0–1.0 range across marketed peptides, so it separates nothing. The irritation proxy is a
hydrophobic moment, which is undefined below two helical turns and therefore inert in four of
five programmes. See [docs/METHODS.md](docs/METHODS.md) §6.

### What this is not

No structure, no docking, no folding, no learned model, no wet-lab data. The output is a
prioritised, chemically annotated hypothesis list with its null distributions attached. None
of the five mechanisms has been demonstrated for any sequence here.

Full method, protease rules, coordination chemistry and the two author-introduced defects the
nulls caught: [docs/METHODS.md](docs/METHODS.md). Numbers: [results/RESULTS.md](results/RESULTS.md).

MIT.
