# R environment

The analysis was developed with R 4.5 and Seurat 5. The exact package versions visible when the final annotated object was checked are recorded in `session-info.txt`.

Run `Rscript environment/install_packages.R` for CRAN and Bioconductor dependencies. The following packages may require their current upstream installation instructions because their release channels change independently:

- CellChat
- monocle3
- SeuratWrappers
- nichenetr

For archival reproduction, use the package versions in `session-info.txt` and record any substitutions in a new release.
