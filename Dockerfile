FROM bioconductor/release_base2:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('knitr')"
RUN Rscript -e "install.packages('tidyverse')"
RUN Rscript -e "BiocManager::install('biomaRt')"
RUN Rscript -e "remotes::install_local()"

# # To start a R session with cTRAP installed:
# docker run -ti [docker image] R
# library(R)
