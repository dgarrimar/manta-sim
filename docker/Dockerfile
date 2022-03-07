FROM rocker/r-ver:3.5.2

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for manta-sim (sim0)"

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# R:devtools dependencies
       libcurl4-gnutls-dev \
       libssl1.0-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
       libgit2-dev \
       pkg-config \
# R:car dependencies
       libnlopt-dev \
# nextflow trace (ps)
       procps \
# other
       wget \
       bc \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## R:copula dependencies
RUN wget ftp://ftp.gnu.org/gnu/gsl/gsl-latest.tar.gz \
    && tar -xf gsl-latest.tar.gz \
    && cd gsl-* \
    && ./configure \
    && make \
    && make install

## Install R packages
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.3.tar.gz", repos=NULL, type="source");'\
'install.packages(c("devtools", "optparse", "car", "CompQuadForm", "vegan", "copula"), repos="http://cloud.r-project.org/")'

## Install MANTA
RUN R -e 'devtools:::install_github("dgarrimar/manta")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash


