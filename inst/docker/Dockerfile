# docker build -t jackgisby/reducedexperiment:latest ./inst/docker/pkg/
# docker push jackgisby/reducedexperiment:latest

# Get Bioconductor development image
FROM bioconductor/bioconductor_docker:devel

# Set metadata
LABEL name="jackgisby/reducedexperiment" \
      url="https://github.com/jackgisby/ReducedExperiment" \
      maintainer="jackgisby@gmail.com" \
      description="Docker image for the ReducedExperiment R package" \
      license="GPL-3.0"

# Install the package and dependencies
RUN R -e 'options(repos = BiocManager::repositories()); \
          devtools::install_github("jackgisby/ReducedExperiment@devel", \
                                 dependencies = TRUE, \
                                 upgrade = "never")'

# Verify installation
HEALTHCHECK --interval=1m --timeout=3s \
  CMD R -q -e 'library(ReducedExperiment)' || exit 1
