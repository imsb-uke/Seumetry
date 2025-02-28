# Use an base image with R 4.4.2 and R-Studio Server installed
FROM rocker/rstudio:4.4.2

# Install system libraries
RUN apt-get update && apt-get install -y \
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
    && apt-get clean

# Install Seumetry
RUN Rscript -e "install.packages('devtools')"
RUN Rscript -e "devtools::install_github('imsb-uke/Seumetry')"