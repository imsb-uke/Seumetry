# Seumetry
Recent progress in flow and mass cytometry technologies enables the simultaneous measurement of over 50 parameters for each cell. The concomitant increase in data complexity requires the usage of advanced data analysis, as conventional analysis methods are time-consuming and fail to capture unknown or minor cell populations. Advances in single-cell RNA sequencing (scRNAseq) technologies led to the development of highly sophisticated computational analysis tools, some of which could significantly improve certain aspects of cytometry data analysis. Here, we present Seumetry, a framework that combines flow and mass cytometry data-specific analysis routines with the capabilities of Seurat, a powerful tool for the analysis of scRNAseq data. Seumetry offers advanced quality control, visualizations, and differential abundance and expression analysis. Together, Seumetry provides a scalable framework for the analysis of high-dimensional cytometry data allowing seamless integration into commonly used scRNAseq analysis tools.  

## Installation
Seumetry can be installed via devtools:
```{r}
install.packages("devtools")
devtools::install_github("imsb-uke/Seumetry")
```

## Tutorial
Please see vignette for instructions on how to analyse cytometry data using Seumetry.  
Data used for the vignette is available through zenodo: https://doi.org/10.5281/zenodo.11935872

## Citation
Please use this citation when using Seumetry in your projects:  
Unpublished.
