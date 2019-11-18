FROM bioconductor/release_base2:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('knitr')"
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "BiocManager::install('biomaRt')"

# Copy package source code
WORKDIR cTRAP
ADD appveyor.yml .
ADD codecov.yml .
ADD CONDUCT.md .
ADD data data
ADD DESCRIPTION .
ADD Dockerfile .
ADD LICENSE .
ADD man man
ADD NAMESPACE .
ADD NEWS.md .
ADD R R
ADD README.md .
ADD tests tests
ADD vignettes vignettes

# Install dependencies
RUN Rscript -e "install.packages(remotes::local_package_deps('.'))"

# Install R package from source
RUN Rscript -e "remotes::install_local()"

# # To start an R session with cTRAP installed:
# docker run -ti [docker image] R
# library(R)
