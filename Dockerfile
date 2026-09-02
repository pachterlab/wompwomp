# Reproducible environment for regenerating the figures in
# "Optimizing alluvial plots" (paper_figures/).
#
# Build (from the wompwomp repo root):
#   docker build -t wompwomp:figures .
#
# Run (RStudio Server on http://localhost:8787, user rstudio):
#   docker run --rm -it -p 8787:8787 -e PASSWORD=changeme wompwomp:figures
#
# The scripts in /home/rstudio/wompwomp/paper_figures write their PDFs/PNGs to
# paper_figures/output/. Copy the finished images into the manuscript repo's
# figures_and_tables/ directory.
#
# NeighborNet is implemented natively in R, so there is no Python/conda setup.

FROM rocker/tidyverse:4.5.1

# System libraries: libpng/jpeg/tiff + freetype/harfbuzz/fribidi/fontconfig are
# needed by png, gridtext, systemfonts, ragg, ggfittext and ggrastr; libglpk for
# igraph; the rest for the usual R build toolchain.
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        libcurl4-openssl-dev \
        libxml2-dev \
        libpng-dev \
        libjpeg-dev \
        libtiff5-dev \
        zlib1g-dev \
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfontconfig1-dev \
        libglpk-dev \
    && rm -rf /var/lib/apt/lists/*

# CRAN dependencies not already in rocker/tidyverse.
RUN R -q -e "install.packages(c('BiocManager', 'TSP', 'ggalluvial', 'ggforce', 'ggfittext', 'ggrastr', 'data.table', 'igraph', 'mclust', 'sessioninfo'))"

# Bioconductor: BiocStyle for the vignette theme, DuoClustering2018 for the
# clustering-comparison dataset.
RUN R -q -e "BiocManager::install(c('BiocStyle', 'DuoClustering2018'), update = FALSE, ask = FALSE)"

# Install wompwomp (which now includes plot_alluvial()) from the build context.
COPY . /home/rstudio/wompwomp
RUN R CMD INSTALL --no-build-vignettes /home/rstudio/wompwomp

# Fail the build if plot_alluvial() cannot render.
RUN Rscript /home/rstudio/wompwomp/paper_figures/_smoke_test.R

RUN chown -R rstudio:rstudio /home/rstudio/wompwomp
WORKDIR /home/rstudio/wompwomp/paper_figures
