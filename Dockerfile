# ---- Multi-Omics Core Pipeline ----
# R 4.5.2 + Bioconductor 3.22 + renv-managed packages
FROM rocker/r-ver:4.5.2

LABEL maintainer="Multiomics Core" \
      description="Automated single- and multi-omics analysis pipeline"

# System libraries required by R packages
# (DESeq2, biomaRt, plotly, openxlsx, rmarkdown, ggiraph, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libglpk-dev \
    libgmp3-dev \
    libcairo2-dev \
    libxt-dev \
    pandoc \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /app

# ---- Layer 1: renv infrastructure (rarely changes) ----
RUN R -e "install.packages('renv', repos = 'https://cloud.r-project.org')"
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
COPY renv/settings.json renv/settings.json

# Restore all packages (cached unless renv.lock changes)
RUN R -e "renv::restore(prompt = FALSE)"

# ---- Layer 2: Project code (changes frequently, fast rebuild) ----
COPY R/ R/
COPY run.R run.R
COPY _targets.R _targets.R

# Create mount points for data and outputs
RUN mkdir -p /app/data /app/outputs /app/config

# Run as non-root user
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Default: interactive wizard
CMD ["Rscript", "run.R", "--new"]
