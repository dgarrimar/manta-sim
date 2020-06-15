FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# devtools dependencies
       libcurl4-openssl-dev \
       libssl-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
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
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan"), repos="http://cloud.r-project.org/")'

## Install mlm
RUN R -e 'devtools:::install_github("dgarrimar/mlm")'

## Install python packages
RUN pip3 install pandas

## Change ENTRYPOINT
ENTRYPOINT /bin/bash

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
RUN wget --quiet https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz \
    && gunzip gemma-0.98.1-linux-static.gz \
    && chmod +x gemma-0.98.1-linux-static \
    && mv gemma-0.98.1-linux-static /usr/local/bin/gemma