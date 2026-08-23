# CS-NO hydrogel in acute kidney injury

This repository contains the processed single-cell RNA-sequencing data and analysis code supporting the manuscript:

> **Single-Cell Transcriptomics Reveals a Transitional Proximal Tubule State Targeted by a Nitric Oxide Hydrogel to Attenuate AKI-to-CKD Progression**

The study evaluates a thermosensitive, beta-galactosidase-responsive chitosan-nitric oxide hydrogel (CS-NO) delivered to the renal subcapsular space after ischemia-reperfusion injury. The single-cell analysis identifies an injury-associated transitional proximal tubule population (NR-PTC/PTC-2) whose representation among recovered cells is reduced after CS-NO treatment.

## Repository contents

```text
.
├── data/
│   ├── processed/          # Cell-level metadata tracked in Git
│   ├── raw/                # Extract the 10x matrix release asset here
│   └── release/            # Release-asset checksums and manifests
├── docs/                   # Data dictionary and analysis notes
├── environment/            # Package installation and session information
├── results/
│   ├── key_tables/         # Tables supporting reported single-cell results
│   └── Supplementary_Data_1_CellChat.xlsx
└── scripts/                # Numbered R analysis workflow
```

## Data

The `v1.0.0` release contains `cs-no-10x-filtered-matrices-v1.0.0.tar.gz`, which provides the Cell Ranger filtered feature-barcode matrices for the PBS-treated AKI control and CS-NO-treated AKI libraries.

After downloading the matrix asset, extract it from the repository root:

```bash
tar -xzf cs-no-10x-filtered-matrices-v1.0.0.tar.gz
```

This creates `data/raw/Ctrl/filtered_feature_bc_matrix/` and `data/raw/NO/filtered_feature_bc_matrix/`.

The repository itself includes compressed cell metadata and the key composition, bootstrap, macrophage-module, and FDR-corrected CellChat tables. See [the data dictionary](docs/data_dictionary.md) and [the release manifest](data/release/MANIFEST.tsv).

## Experimental design and statistical scope

Kidneys from three mice were pooled for each condition before tissue dissociation, producing one pooled single-cell library per condition. The matrix/QC step retained 16,781 control cells and 8,805 CS-NO cells. After joint clustering and removal of one low-quality cluster (cluster 13), the final annotated object contained 16,638 cells from the PBS-treated control library and 8,737 cells from the CS-NO-treated library.

Cell-type proportions therefore describe the composition of recovered cells and are not estimates of absolute cell abundance in kidney tissue. Fisher's exact tests and the cell-level bootstrap are reported as descriptive sensitivity analyses; neither method provides sample-level biological replication. The NR-PTC/PTC-2 population decreased from 16.8% to 5.3% of recovered cells, corresponding to a 68% reduction in relative representation. The single-cell data cannot distinguish cell loss from transition to another transcriptional state.

## Reproducing the analysis

The workflow was developed in R with Seurat 5. Run commands from the repository root.

```bash
Rscript environment/install_packages.R
Rscript scripts/01_quality_control.R
Rscript scripts/02_integration_clustering.R
Rscript scripts/03_differential_expression.R
Rscript scripts/04_trajectory_analysis.R
Rscript scripts/05_cellchat_analysis.R
Rscript scripts/06_signature_analyses.R
Rscript scripts/07_ptc_injury_scoring.R
Rscript scripts/08_composition_and_nr_ptc.R
Rscript scripts/09_cellchat_visualization.R
Rscript scripts/10_reproduce_manuscript_figures.R
Rscript scripts/11_verify_reported_values.R
```

For a faster check that does not rebuild the Seurat object, begin with script 08 or 11; both operate on the tracked metadata and result tables. Scripts use repository-relative paths and fixed random seeds where resampling is performed.

## Data availability

This public repository provides processed 10x matrices, cell-level annotations, analysis code, and supporting result tables. Raw sequencing reads are not contained in this repository and should be deposited in an appropriate sequence archive; the accession will be added here and to the manuscript when available.

## Citation

The manuscript is under revision. Citation details and the article DOI will be added after acceptance. Until then, cite this repository using [CITATION.cff](CITATION.cff).

## Licensing

Code is released under the [MIT License](LICENSE). Data tables and processed matrices are released under [CC BY 4.0](DATA_LICENSE.md), subject to the terms stated there.

## Contact

Zongjin Li Lab, School of Medicine, Nankai University  
Website: [zongjinlab.com](https://www.zongjinlab.com)  
Email: zongjinli at nankai.edu.cn
