FROM rocker/tidyverse:4.5.1

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        # Core R package build dependencies
        build-essential \
        libcurl4-openssl-dev \
        libxml2-dev \
        libpng-dev \
        libjpeg-dev \
        zlib1g-dev \
        # Fonts and text shaping for systemfonts/gridtext/ggfittext
        libfreetype6-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfontconfig1-dev \
        # GLPK for igraph
        libglpk-dev \
        # Python support for reticulate
        python3 python3-dev python3-venv python3-pip && \
    rm -rf /var/lib/apt/lists/*

# Install wompwomp and set up Python environment
 RUN R -e "if (!require('remotes', quietly = TRUE)) install.packages('remotes'); \
            remotes::install_github('pachterlab/wompwomp')"
