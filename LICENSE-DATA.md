# Data license and source attribution

Copyright © 2026 Zichen Sheng for the original curation, derived study outputs, and synthetic data.
These contributions are licensed under [Creative Commons Attribution 4.0 International](https://creativecommons.org/licenses/by/4.0/).
Cite the versioned software release and give credit to the corresponding original sources.

| Files | Source attribution and scope of the grant |
| --- | --- |
| `data/signatures_201.tsv`, `data/signature_sources.tsv`, `data/weighted_formulas_104.tsv` | CC BY 4.0 applies to the authors' original curation and compilation. Attribute the source gene lists, coefficients and publications using the DOIs and PMIDs in `data/signature_sources.tsv`. Any applicable terms of those sources continue to govern their material. |
| `data/consensus_loadings.tsv` | Study-fitted loadings from CGGA693. Attribute the Chinese Glioma Genome Atlas and Zhao et al. (2021). [CGGA](https://www.cgga.org.cn/download.jsp) |
| `data/formula_gene_reference_sd.tsv`, `data/matching_pool_features.tsv` | Study-derived per-gene aggregates. Attribute CGGA693 and CGGA325, the GLASS Consortium and cBioPortal, the TCGA Research Network and UCSC Xena, and the CPTAC discovery study and Proteomic Data Commons. [CGGA](https://www.cgga.org.cn/download.jsp), [GLASS study](https://www.cbioportal.org/study/summary?id=difg_glass), [cBioPortal data terms](https://docs.cbioportal.org/user-guide/faq/), [UCSC Xena](https://xenabrowser.net/datapages/), [CPTAC PDC000205](https://proteomic.datacommons.cancer.gov/pdc/study/PDC000205) |
| `data/reference_gene_sets.json` | Study-derived coordinate definitions. Attribute the [GBmap study and source atlas](https://pmc.ncbi.nlm.nih.gov/articles/PMC12526130/). See the data provenance in the README. |
| `data/protein_coding_genes.txt` | Gene symbols extracted from NCBI Gene annotations via Marc Carlson's [org.Hs.eg.db package](https://bioconductor.org/packages/org.Hs.eg.db/), distributed upstream under Artistic-2.0. The CC BY grant covers the authors' contributions only. |
| `examples/synthetic_example.csv`, `examples/synthetic_samples.csv` | Author-generated demonstration inputs. |

The data license applies only to rights held by the authors. Original sources and externally obtained data
remain subject to their source terms. The original Python code is distributed separately under [MIT](LICENSE).
