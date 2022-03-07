FROM rocker/r-ver:3.5.2

LABEL maintainer "Diego Garrido-Martín <diego.garrido@crg.eu>" \
      version "1.0.0" \
      description "Docker image for manta-sim (sim, sim2)"

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# R:devtools dependencies
       libcurl4-gnutls-dev \
       libssl1.0-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
       libgit2-dev \
       pkg-config \
# bcftools dependencies
       libbz2-dev \
       liblzma-dev \
# nextflow trace (ps)
       procps \
# python3
       python3 \
       python3-pip \
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
RUN R -e 'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan", "BEDMatrix", "copula", "R.utils", "mvnfast"), repos="http://cloud.r-project.org/")'

## Install python3 packages
RUN pip3 install pandas scipy

## Install bcftools and HTSlib

RUN wget --quiet https://github.com/samtools/bcftools/releases/download/1.10.2/bcftools-1.10.2.tar.bz2 \
    && tar jxf bcftools-1.10.2.tar.bz2 \
    && cd bcftools-1.10.2 \
    && ./configure \
    && make \
    && make install

RUN cd bcftools-1.10.2/htslib-1.10.2 \
    && ./configure \
    && make \
    && make install

## Install PLINK2
RUN wget --quiet https://github.com/chrchang/plink-ng/releases/download/v2.00a2.3/plink2_linux_x86_64.zip \
    && unzip plink2_linux_x86_64.zip \
    && mv plink2 /usr/local/bin/

## Install GEMMA
RUN wget --quiet https://github.com/genetics-statistics/GEMMA/releases/download/0.98.3/gemma-0.98.3-linux-static.gz \
    && gunzip gemma-0.98.3-linux-static.gz \
    && chmod +x gemma-0.98.3-linux-static \
    && mv gemma-0.98.3-linux-static /usr/local/bin/gemma

## Install MANTA
RUN R -e 'devtools:::install_github("dgarrimar/manta")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash
