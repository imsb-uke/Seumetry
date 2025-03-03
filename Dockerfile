# Use an base image with R 4.4.2 and R-Studio Server installed
FROM rocker/rstudio:4.4.2

# Install system libraries
RUN apt-get update && apt-get install -y \
    libboost-all-dev \
    libglpk40 \
    libglpk-dev \
    libzmq3-dev \
    libgsl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libhdf5-dev \
    build-essential \
    && apt-get clean

# Install essentials
RUN Rscript -e "install.packages(c('devtools', 'BiocManager'))"

# Install dependencies that are better installed through Bioconductor directly
RUN Rscript -e "BiocManager::install(c('cytolib', 'CATALYST'))"

# Install Seumetry
RUN Rscript -e "devtools::install_github('imsb-uke/Seumetry')"
