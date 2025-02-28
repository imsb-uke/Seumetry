# Seumetry
Recent progress in flow and mass cytometry technologies enables the simultaneous measurement of over 50 parameters for each cell. The concomitant increase in data complexity requires the usage of advanced data analysis, as conventional analysis methods are time-consuming and fail to capture unknown or minor cell populations. Advances in single-cell RNA sequencing (scRNAseq) technologies led to the development of highly sophisticated computational analysis tools, some of which could significantly improve certain aspects of cytometry data analysis. Here, we present Seumetry, a framework that combines flow and mass cytometry data-specific analysis routines with the capabilities of Seurat, a powerful tool for the analysis of scRNAseq data. Seumetry offers advanced quality control, visualizations, and differential abundance and expression analysis. Together, Seumetry provides a scalable framework for the analysis of high-dimensional cytometry data allowing seamless integration into commonly used scRNAseq analysis tools.  

## Installation
Seumetry can be installed via devtools:
```{r}
install.packages("devtools")
devtools::install_github("imsb-uke/Seumetry")
```

## Tutorial
Please see our tutorial for instructions on how to analyse cytometry data using Seumetry: https://imsb-uke.github.io/Seumetry/vignettes/vignette.html  
Data used in this tutorial is available through Zenodo: https://doi.org/10.5281/zenodo.11935872

## Docker
A Docker file and image are supplied with in this repository. This image is based on rocker/rstudio and contains rstudio-server and an installation of the current release of Seumetry. Rstudio can be accessed in the browser at http://localhost:8787.  
You can download the Docker image and run a container:
```{bash}
docker pull ghcr.io/imsb-uke/seumetry
docker run --name {CONTAINER_NAME} -ti -e DISABLE_AUTH=true -p 8787:8787 -v {PROJECT_DIR}:/home/rstudio/workspace seumetry
```

## Publication
Seumetry is currently available as a preprint on bioRxiv: https://www.biorxiv.org/content/10.1101/2024.07.23.604747v1

## Citation
Please use this citation when using Seumetry in your projects:  
Borggrewe, M. et al 2024. Seumetry: a versatile and comprehensive R toolkit to accelerate high-dimensional flow and mass cytometry data analysis. bioRxiv doi: 10.1101/2024.07.23.604747
