## Package outline

This code supports the manuscript Carr et al. BioRxiv 2021 https://biorxiv.org/cgi/content/short/2021.03.04.433942v1

This repository was created using RStudio and github version control.

R code is provided in Rmarkdown. Each .Rmd is 'knitted' to pdf.

To execute locally, it requires RStudio and R verion 3.6.1. Code has not been tested on later versions of R. (R>4 is likely to give peculiar results as its handling of file import is changed.)

### cohort_2016_17 directory

This contains Rmarkdown files for:

- read counting using `Rsubread`
- `SingleCellExperiment` object creation
- basic QC
- annotation with the index flow dataset and removal of NTC wells
- Differential expression and differential abundance analyses
- Diffusion Map generation

The `data` sub-directory contains the output from these Rmd files.

### Turner dataset reanalysis directory

This contains Rmarkdown files for:

- downloading fastq + cellranger aligning and counting
- Seurat QC & analysis
- VDJ fate mapping with immunarch
- Immcantation suite usage to calculate SHM for a single sample (IgD- PBMC at day 28)

### scripts directory

Helper script to import VDJPuzzle output into R.

### figures directory

Each figure is drawn in turn.
These Rmarkdown files require the cohort_2016_17 data sub-directory.
Where code execution is quick (<1 minute), the code for that figure is self-contained.
Where code execution is longer, a separate processing Rmd has been written in the `cohort_2016_17` directory.
Draft figure legends are included on the second page of each generated pdf.


