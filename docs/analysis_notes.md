# Analysis notes

## Reference genome and preprocessing

- Cell Ranger version recorded in the matrix headers: 5.0.0
- Reference genome used in the study: mouse mm10
- Initial matrix/QC retention: 500-6,000 detected genes, 1,000-40,000 UMIs, mitochondrial transcript proportion no greater than 20%, and log10 genes-per-UMI complexity of at least 0.8
- One low-quality cluster (cluster 13) was removed after joint clustering, yielding the final 16,638 control and 8,737 CS-NO cells used for downstream analysis

## Integration and annotation

- Seurat 5 log-normalisation
- 3,000 highly variable genes
- PCA followed by Harmony correction for the joint embedding
- First 30 principal components used downstream
- Differential expression performed on the normalised RNA assay, not the Harmony embedding
- Twelve major cell populations, including PTC-1, PTC-2, and PTC-3

## Compositional analysis

Fisher's exact tests with Benjamini-Hochberg correction are retained as descriptive cell-level comparisons. A 2,000-resample cell-level bootstrap estimates the sensitivity interval for the change in recovered PTC-2 fraction. These procedures do not estimate between-animal variability because there is one pooled library per condition.

## Cell-cell communication

CellChat permutation P values are adjusted within each condition using the Benjamini-Hochberg method. Interactions with FDR below 0.05 are labelled significant. The reported change from 359 to 165 significant source-target ligand-receptor combinations refers specifically to the CCL pathway and is not a global reduction in all inferred communication.
