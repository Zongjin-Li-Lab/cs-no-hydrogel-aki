# Data dictionary

## Conditions

| Repository label | Biological condition |
|---|---|
| `Control` / `Ctrl` | PBS-treated kidney at day 3 after ischemia-reperfusion injury |
| `NO_treated` / `NO` | CS-NO-treated kidney at day 3 after ischemia-reperfusion injury |

Kidneys from three mice were pooled before dissociation for each condition, yielding one pooled library per condition.

## Raw matrix release

Each condition is provided in standard 10x Genomics Matrix Market format:

- `matrix.mtx.gz`: sparse gene-by-cell count matrix
- `features.tsv.gz`: Ensembl feature identifier, gene symbol, and feature type
- `barcodes.tsv.gz`: cell barcode identifiers

Both matrices contain 32,285 features. The control matrix contains 24,978 Cell Ranger-filtered barcodes, and the CS-NO matrix contains 16,197.

## `data/processed/cell_metadata.csv.gz`

Cell-level annotations for the 25,375 cells retained after study-specific matrix/QC thresholds, joint clustering, and removal of one low-quality cluster. Important columns include:

- `nCount_RNA`, `nFeature_RNA`, `percent.mt`: quality-control measurements
- `sample`, `condition`: library and treatment labels
- `seurat_clusters`: final clustering identifier
- `cell_type`: study cell-type annotation
- `singler_immgen`, `singler_mouse`: reference-based labels used as supporting evidence

The final downstream Seurat object includes additional module scores described in the scripts and is distributed as a release asset.

## Key result tables

- `library_qc_metrics.csv`: library-level recovery and median QC measurements
- `celltype_composition.csv`: recovered-cell counts and proportions by cell type
- `nr_ptc_bootstrap.csv`: observed and bootstrap sensitivity estimates for PTC-2/NR-PTC representation
- `macrophage_module_scores.csv`: descriptive M1- and M2-associated module scores in recovered macrophages
- `cellchat_all_bh.csv.gz`: CellChat interactions with Benjamini-Hochberg adjusted P values
- `cellchat_ccl_fdr.csv`: CCL-pathway subset of the corrected interaction table

## Interpretation limits

The two libraries are pooled condition-level samples without replicate libraries. Recovered-cell fractions are compositional and may be affected by tissue composition, dissociation survival, loading, capture, and recovery. Cell-level P values must not be interpreted as tests of animal-level treatment effects.
