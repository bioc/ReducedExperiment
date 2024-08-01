# docker build -t jackgisby/reducedexperiment ./inst/docker/pkg/
# docker push jackgisby/reducedexperiment

# Get Bioconductor development image
FROM bioconductor/bioconductor_docker:devel

# Install the package
RUN R -e 'devtools::install_github("jackgisby/ReducedExperiment")'
