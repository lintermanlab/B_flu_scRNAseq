## Repository outline

[![DOI](https://zenodo.org/badge/344757391.svg)](https://zenodo.org/badge/latestdoi/344757391)

This code supports the manuscript Carr _et al._ BioRxiv 2021 https://biorxiv.org/cgi/content/short/2021.03.04.433942v1.
The corresponding GEO accession is GSE167823 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167823.
This repository has been updated (26 September 2022), to include changes after peer review.

- R code is provided in Rmarkdown. Each .Rmd is 'knitted' to pdf.

- To execute locally, it requires RStudio and R verion 3.6.1. Code has not been tested on later versions of R. (R>4 is likely to give peculiar results as its handling of file import is changed. In review, miloR and later Immcantation analysis was specifically performed in R>4 due to dependencies.).


### cohort_2016_17 directory

This contains Rmarkdown files for:

- read counting using `Rsubread`
- `SingleCellExperiment` object creation
- basic QC
- annotation with the index flow dataset and removal of NTC wells
- Differential expression and differential abundance analyses
- Diffusion Map generation

The `data` sub-directory contains the output from these Rmd file. Some of these files are too large to host on github. The final `SingleCellExperiment` object is available.


### Turner dataset reanalysis directory

The Turner _et al._ dataset is described in _Nature_ 2020 https://www.nature.com/articles/s41586-020-2711-0, and available from GEO accession GSE148633 http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148633.

This contains Rmarkdown files for:

- downloading fastq + cellranger aligning and counting
- Seurat QC & analysis
- VDJ fate mapping with immunarch
- Immcantation suite usage to calculate SHM for a single sample (IgD- PBMC at day 28)

The Seurat objects themselves are too large for github. They would need re-creation from via fastq->cellranger->Seurat.

### scripts directory

Helper script to import VDJPuzzle output into R.

### figures directory

Each figure is drawn in turn.
These Rmarkdown files require the `cohort_2016_17/data` directory. Some intermediate, large files from this directory are not duplicated on github - these files would need local regeneration for all figures to be drawn. Code to re-generate these intermediate steps is in `cohort_2016_17` directory.


Where code execution is quick (<1 minute), the code for that figure is self-contained.
Where code execution is longer, a separate processing Rmd has been written in the `cohort_2016_17` directory.
Draft figure legends (see the accepted manuscript for final legends) are included on the second page of each generated pdf.

The spectral cytometry data was analysed separately, using code available elsewhere. 
