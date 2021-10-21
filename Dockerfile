FROM bioconductor/bioconductor_docker:latest
MAINTAINER Nuno Agostinho <nunodanielagostinho@gmail.com>

RUN apt-get update && apt-get -y upgrade && apt-get -y autoremove

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "remotes::install_github('nuno-agostinho/floweRy')"

# Copy package source code
WORKDIR cTRAP
ADD . .

# Install R package from source
RUN Rscript -e "remotes::install_local()"
RUN rm -rf *
WORKDIR /data

# # To start an R session with cTRAP installed:
# docker run -ti [docker image] R
# library(cTRAP)

# # To start an RStudio session on http://localhost:8787 and enter with user rstudio and password bioc
# docker run -e PASSWORD=bioc \
#	  -p 8787:8787 \
#	  [docker image]
