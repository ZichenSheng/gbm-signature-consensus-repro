# Reproducibility package: published GBM signatures and their ecological consensus

Code accompanying the manuscript

> **Published glioblastoma prognostic signatures share an ecological consensus that underrepresents a
> driver-associated malignant-programme balance**.

This study-specific package contains selected methods for scoring published glioblastoma signatures
and assessing the ecological and malignant programmes represented by their consensus. It includes
machine-readable definitions, a synthetic example and tests. Patient-level inputs are obtained from
the original repositories.

## Coverage and limitations

The package supports selected computations behind Figs 1A–D, 2A, 3B–C, 5A–E and 6A–C
from independently assembled cohort matrices.

| Manuscript statement (Methods; figure) | Released function(s) | Verified against frozen outputs |
|---|---|---|
| RNA-seq processing: detected protein-coding genes, median-of-ratios size factors, log2, per-gene z (Supplementary Methods) | `signature_analysis.normalize_rnaseq`, `zscore_genes`, `harmonize_symbols` | CPTAC-3 z-matrix identical (15,549 genes × 160 tumours) |
| Membership scores of 201 signatures; weighted formulas (Signatures, scores and reference coordinates; Fig. 1D) | `membership_scores`, `scorable_signatures`, `weighted_formula_scores` | per-formula r(weighted, membership) and collection-level \|r\| identical in CPTAC-3 and CGGA325 |
| Literature consensus: first-component loadings fixed in CGGA693 and applied unchanged (Fig. 1) | `fit_consensus_loadings`, `consensus_score`, `data/consensus_loadings.tsv` | consensus identical (≤ 2×10⁻¹⁵) in all 8 cohort-platforms |
| Reference coordinates derived from GBmap by fixed rules; hypoxia with expression-matched controls; relative developmental programme | `derive_reference_gene_sets`, `panel_from_sets`, `score_reference_coordinates`, `relative_programme_scores` | gene sets re-derived identically from the GBmap pseudobulk; coordinates identical in 8 cohort-platforms |
| Cross-cohort reproducibility of pairwise score correlations (Fig. 1A) | `statistical_validation.structure_reproducibility`, `FeatureMatcher`, `gene_features` | ρ = 0.605 for GLASS vs CGGA693 and its first null draw matched the frozen analysis. The 200-signature sensitivity comparison of CGGA325 vs CGGA693 gave ρ = 0.734 |
| Ecological share = in-sample R² on five ecological coordinates (Fig. 1B, 2A) | `ecological_share`, `collection_scores`, `leading_components` | all 8 consensus values of Table S2 (0.427–0.747); CGGA325 first-component value 0.844 |
| Expression-matched (N1) and co-expression-matched (N2) random collections (Fig. 2, 3E) | `literature_sets_on_pool`, `expression_matched_collection`, `coexpression_matched_collection`, `pooled_standardized_expression`, `within_set_coherence` | 183 pool-restricted sets, N1[0], N1[1] and N2[0] identical gene for gene; within-set correlation 0.126 / 0.032 / 0.122 |
| Leave-one-cohort-out ridge recovery; held-out R² around the held-out mean (Fig. 3B,C) | `loco_ridge_recovery`, `loco_one_variable` | per-fold R² identical for balance (0.822; 0.907 genes outside panel; 0.013 consensus), glial, proliferation, myeloid |
| Driver calls, DNA purity, driver-associated variance (adjusted R² of ranked score), cluster bootstrap, fixed-effect meta-analysis, purity and random-contrast controls (Fig. 5A–E) | `copy_number_driver_calls`, `gistic_driver_calls`, `vaf_purity`, `usable_drivers`, `driver_associated_variance`, `driver_variance_bootstrap`, `purity_adjusted_driver_variance`, `matched_random_contrasts`, `dosage_region_genes`, `fixed_effect_meta` | CPTAC-3 0.228 (0.137–0.367) vs 0.039, Δ 0.189 (0.043–0.317); 68 never-analysed 0.220; Ghannam 0.056; meta 0.102 (0.021–0.184); TCGA arrays; beyond purity 0.225/0.002; 1,000 random contrasts — all identical |
| Between-patient share ICC(1) and its contrast with seven ecological coordinates, patient-cluster bootstrap (Fig. 6A–C) | `icc1`, `icc_contrast`, `bootstrap_icc`, `bootstrap_icc_contrast` | Mathur 0.699 (0.306–0.881), Δ 0.463 (0.097–0.720); GBM-Space 0.802, Δ 0.435 (0.220–0.612) — identical |

70 of 70 numerical checks against frozen outputs passed in Python 3.13.11 with the four pinned
dependencies. The comparisons cover the computations and selected seeded draws specified above.
The checks used local cohort inputs. Those inputs and the private comparison harness are excluded.

### Omitted workflows

The following analyses and end-to-end workflows are outside this release:

* **Cohort assembly.** Download, eligibility filtering and unit conversion for the eight cohort-platforms,
  the TCGA held-out array pairing, the Ghannam libraries and the Mathur and IvyGAP data. The
  GBM-Space pseudobulk preparation step is included. The eligibility rules are stated in the
  Supplementary Methods and summarised below.
* **Single-nucleus analyses** (Fig. 4A,B): site–composition correlations and the exploratory Shapley
  composition–expression attribution.
* **Protein and imaging** (Fig. 4C,D): CD3 immunohistochemistry, MIBI read-outs and the exploratory bulk
  proteome.
* **Constructed signatures** (Fig. 2C; 5,420 LASSO-Cox and comparator models), training-population
  comparisons (Fig. 2D), concentration ratios (Fig. 2B), spectrum position and Tucker congruence
  (Fig. 3D,E), and the platform-agreement read-out.
* **Programme shifts and transcriptome-wide driver signatures** (Fig. 5F,G) and the exploratory
  DNA-methylation class analyses.
* **Spatial anatomy** (Fig. 6D,E): Visium spot annotation and IvyGAP histological variance shares, and
  the sensitivity analyses of Table S5 other than the primary contrast.
* **The Fig 1A CGGA registry comparison and 1,000-panel null.** The figure uses 201 signatures
  (20,100 pairs, observed ρ = 0.730923) from the earlier registry analysis. That pipeline and its
  null summaries are available in the
  [earlier registry software archive](https://github.com/ZichenSheng/GBM-signature-measurement/tree/v1.1.0).
  The package's CGGA sensitivity check uses the 200 signatures also evaluable in GLASS. For the
  separate 500-panel GLASS comparison, the observed statistic and only the first null draw were
  verified with this package.

## Installation and test

Requires Python ≥ 3.11. A clean Python 3.13.11 environment passed 15/15 tests and ran the
synthetic example. `requirements.txt` pins NumPy 2.4.4, pandas 3.0.2, SciPy 1.17.0 and
scikit-learn 1.8.0.

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python -m unittest discover -s tests -v
python scripts/run_example.py
```

The tests check data-file checksums and counts, numerical identities (size factors, z-scores, zero-sum
contrast invariance, closed-form ICC, adjusted R², fixed-effect pooling), deterministic resampling, and
that the expression-matched collection N1 regenerated from the supplied pool features is identical to the
frozen production collection (SHA-256). `run_example.py` exercises every released procedure on
**synthetic** data and writes to `example_output/`.

> The synthetic example (`examples/`, regenerated by `scripts/make_synthetic_example.py`, seed 1) is
> simulated from a latent-factor count model in which published genes are made to lean ecological by
> construction. It demonstrates execution only and is not evidence for any result in the manuscript.

## Data files (`data/`)

| File | Content | Source |
|---|---|---|
| `signatures_201.tsv` | `signature_id`, `gene_symbol` — membership of the 201 signatures (1,277 genes, HGNC symbols, file order is significant) | published articles, curated by the authors |
| `signature_sources.tsv` | `signature_id`, `pmid`, `doi`, `year` (197 publications) | bibliographic metadata |
| `weighted_formulas_104.tsv` | `signature_id`, `gene_symbol`, `coefficient` for the 104 recoverable formulas | published articles |
| `consensus_loadings.tsv` | frozen first-component loadings (201 signatures) | this study (CGGA693) |
| `reference_gene_sets.json` | gene sets of the ten reference coordinates, the pan-malignant reference set and programme sets (AC, developmental, MES, cycling, OPC, NPC) | this study (GBmap rules) |
| `protein_coding_genes.txt` | protein-coding symbol universe (MT- genes excluded) | NCBI Gene via Bioconductor org.Hs.eg.db |
| `formula_gene_reference_sd.tsv` | median SD of log2 expression of formula genes across the five development cohorts | aggregate per-gene statistic |
| `matching_pool_features.tsv` | 9,339-gene pool with median mean, SD and detection across the development cohorts | aggregate per-gene statistic |
| `SHA256SUMS` | checksums of the files above | |

The tables in `data/` contain model definitions and aggregate values. The sample-level records in
`examples/` are simulated. Cite the original publications linked by DOI and PMID in
`signature_sources.tsv` when reusing gene lists or coefficients. The extracted protein-coding
symbol list is attributed to NCBI Gene via [Marc Carlson's Bioconductor annotation package](https://bioconductor.org/packages/org.Hs.eg.db/),
which is distributed under Artistic-2.0. The atlas-derived definitions are attributed to the
[GBmap study](https://pmc.ncbi.nlm.nih.gov/articles/PMC12526130/).

## Input formats

* **Expression**: `pandas.DataFrame`, genes (upper-case symbols) × samples. RNA-seq on a linear scale
  (counts, FPKM or TPM as distributed; a log2(x + 1) distribution is converted back with 2ˣ − 1); arrays
  on their distributed log2 scale, standardized directly with `zscore_genes`.
* **Sample sheet** (driver analyses): one row per sample with `patient_id`, 0/1 indicators
  `EGFR_amp`, `PDGFRA_amp`, `CDK4_amp`, `NF1_loss`, optional `dna_purity`; see
  `examples/synthetic_samples.csv`.
* **Copy number** (`copy_number_driver_calls`): absolute gene-level copy number, genes × samples, with an
  autosome mask; TCGA calls use GISTIC2 thresholded values (`gistic_driver_calls`).
* **Multi-region** (`icc1`, `bootstrap_icc_contrast`): samples × coordinates with a patient vector.
* **GBmap** (`derive_reference_gene_sets`): summed counts per (donor × level-4 annotation) group with
  metadata columns `donor`, `l1`, `l2`, `l4`, `n`.

## Exact parameters and seeds

| Procedure | Parameters | Seed |
|---|---|---|
| Membership score | ≥ 3 genes and ≥ 50% of genes measured (≥ 80% for the loading fit and Fig. 1A) | – |
| Weighted formula | Σ β·sd_ref·z, ≥ max(2, 80%) genes | – |
| Reference coordinate | ≥ 5 genes and ≥ 50% of each set; hypoxia controls: 25 expression bins, 10 controls per gene | 20261107 (fresh generator per cohort) |
| N1 | 10 nearest non-published pool genes, shared one-to-one map | 20261108 + 1000 + *i* |
| N2 | 300 nearest candidates; uniform draw among the 3 closest to the target within-set correlation | 20261108 + 9000 + *i* |
| Coordinate recovery (Fig. 3B,C) | RidgeCV, 25 penalties 10⁻¹–10⁵, efficient leave-one-out on training cohorts | – |
| Collection recovery (Fig. 3E) | RidgeCV, 12 penalties 10⁻¹–10⁴ | – |
| Driver variance | drivers with ≥ 5 events and ≥ 5 non-events; 2,000 patient-cluster resamples, resamples with < 2 events of any driver skipped; 1,000 random contrasts (±0.02 expression percentile; panel and dosage-region genes excluded) drawn after the bootstrap from the same generator | 20261108 |
| Meta-analysis | inverse-variance fixed effect, bootstrap SEs | – |
| ICC(1) | one-way ANOVA, k₀ correction, negative variance → 0; Mathur regions with DNA purity ≥ 0.2 | – |
| ICC intervals | 1,000 patient-cluster resamples per coordinate | 20261108 |
| ICC contrast | 2,000 patient-cluster resamples; Mathur first, then GBM-Space from the same generator | 20261110 |
| Pairwise-structure null | k = 30 same-biotype nearest genes in (mean, SD, detection) space | 20261107 |

## Implementation details that refine the written Methods

Stated here so that the released code and the text can be checked against each other:

1. **Loading fit.** The frozen loadings were estimated on per-gene z-scores of linear-scale CGGA693 RSEM
   FPKM (108 adult primary IDH-wild-type WHO grade IV tumours, without the survival requirement that
   defines the 104-patient analysis cohort), using all 201 signatures (each ≥ 80% covered), with the sign
   oriented to the ecological direction. Development RNA-seq scoring uses size-factor/log2 processing.
2. **Fig 1A** used per-gene z-scores of linear-scale CGGA FPKM (CGGA693 n = 108, CGGA325 n = 73)
   and log2 TPM for GLASS (n = 106). The CGGA registry comparison used 201 signatures and 1,000
   null panels. The separate GLASS comparison used the 200 signatures meeting the 80% common
   coverage threshold and 500 null panels.
3. **N2** draws the first gene uniformly from the 300 eligible feature-matched candidates.
   Each subsequent gene is drawn uniformly from the three candidates whose inclusion gives
   the closest mean within-set correlation to the published set.
4. **"All 201 scores"** in recovery: all 201 signatures are scorable in all five development cohorts.

## GBM-Space pseudobulk input

For the GBM-Space calculations in Fig 6A–C, obtain the raw-count snRNA H5AD from the
[GBM-Space portal](https://www.gbmspace.org/). The production pseudobulk step reads the CSR counts
in `X`, gene symbols in `var/SYMBOL`, and the annotation, donor, site and phase categories in
`obs`. It writes the donor- and site-level sums and three group-count tables to a local directory.

```sh
python -m pip install 'h5py==3.15.1'
python scripts/build_gbmspace_pseudobulk.py /path/to/GBM_space_snRNA.h5ad ../gbmspace_output
```

All eight output arrays matched the frozen pseudobulk values exactly, and all three TSVs matched
byte for byte. Obtain the original data and generate outputs locally under the source terms.

## External datasets

Obtain inputs from the original repositories under their access terms. The loadings are derived
from CGGA693. The two per-gene reference tables summarize CGGA693, CGGA325, GLASS, TCGA and CPTAC
discovery across the five development cohorts.

| Dataset | Use | Source / accession |
|---|---|---|
| CGGA mRNAseq_693, mRNAseq_325 | development cohorts; loading cohort | [CGGA data portal](https://www.cgga.org.cn/download.jsp), Zhao et al. (2021) |
| GLASS | development cohort | [cBioPortal study difg_glass](https://www.cbioportal.org/study/summary?id=difg_glass), Varn et al. (2022) and Malta et al. (2023) |
| TCGA GBM RNA-seq, U133A, G4502A, GISTIC2 | development cohort; held-out arrays | UCSC Xena (TCGA legacy hub) |
| CPTAC GBM (discovery) | development cohort | Proteomic Data Commons PDC000205 |
| CPTAC-3 GBM RNA, ASCAT-NGS copy number, WXS mutations | independent DNA population | Genomic Data Commons (CPTAC-3) |
| Ghannam et al. ICB cohort | second independent DNA population | [Ghannam et al. data and code](https://doi.org/10.5281/zenodo.20517473) |
| Mathur et al. 3-D multi-region | multi-region ICC | supplementary tables of Mathur et al., Cell 2024 |
| GBM-Space single-nucleus | multi-region ICC | [GBM-Space](https://www.gbmspace.org/), de Jong et al. (2026) |
| GBmap core atlas | reference coordinates | [CZ CELLxGENE Discover, core GBmap](https://cellxgene.cziscience.com/collections/999f2a15-3d7e-440b-96ae-2c806799c08c) |
| GENCODE v36 gene positions | dosage regions | UCSC Xena probemap |

Selection (Supplementary Methods): CGGA and GLASS adult primary IDH-wild-type grade 4 with survival;
TCGA non-G-CIMP primary tumours; CPTAC discovery tumours with `is_gcimp = FALSE`; CPTAC-3 GBM tumours
without IDH1 R132 / IDH2 R140/R172 mutations, one sample per case with RNA and copy number from the same
sample.

## Citation and reuse

Use the software citation in `CITATION.cff` and the version-specific DOI of the archived release
you used. Cite the accompanying manuscript once published.

## Licenses

The original Python code in `src/`, `scripts/` and `tests/` is licensed under [MIT](LICENSE).
The authors' original data curation, derived aggregate outputs and synthetic examples are licensed
under CC BY 4.0. See the [data license and attribution](LICENSE-DATA.md) for the grant and source
terms that apply to each file.
