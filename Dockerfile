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
COPY wizard.html wizard.html
COPY wizard_defaults.template.json wizard_defaults.template.json
COPY tools/ tools/
COPY utils/ utils/
COPY config/templates/ config/templates/
COPY data/example_proteomics/ data/example_proteomics/
COPY data/example_metabolomics/ data/example_metabolomics/
COPY data/example_lipidomics/ data/example_lipidomics/

# Create mount points for user data + outputs (left empty inside the image —
# bind-mount real data at run time with `docker run -v ...`).
RUN mkdir -p /app/data/user /app/outputs /app/config

# Wizard binds to 0.0.0.0 inside the container so `-p 8080:8080` forwarding
# works from the host browser. Outside Docker, default stays 127.0.0.1.
ENV MULTIOMICS_WIZARD_BIND=0.0.0.0

# Expose the wizard port for `docker run -p 8080:8080`
EXPOSE 8080

# Run as non-root user
RUN useradd -m -u 1000 appuser && chown -R appuser:appuser /app
USER appuser

# Default: HTML wizard on port 8080. Override with `docker run ... <cmd>`
# to use the CLI wizard (`Rscript run.R --new`) or run an existing config
# (`Rscript run.R --config /app/config/foo.yaml`).
CMD ["Rscript", "run.R", "--wizard"]
