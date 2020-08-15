FROM bioconductor/bioconductor_docker:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"

# Copy package source code
WORKDIR cTRAP
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

# Install R package from source
RUN Rscript -e "remotes::install_local()"

# # To start an R session with cTRAP installed:
# docker run -ti [docker image] R
# library(cTRAP)

# # To start an RStudio session on http://localhost:8787 and enter with user rstudio and password bioc
# docker run -e PASSWORD=bioc \
#	  -p 8787:8787 \
#	  [docker image]
