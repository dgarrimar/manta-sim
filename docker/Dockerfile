FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# R:devtools dependencies
       libcurl4-openssl-dev \
       libssl-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
# R:copula dependencies
       libgsl-dev \ 
# bcftools dependencies
       libbz2-dev \
       liblzma-dev \
# nextflow trace (ps)
       procps \
# python3
       python3 \
       python3-pip \
# python2
       python \
       python-pip \
# other
       wget \
       git \
       bc \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan", "BEDMatrix", "copula"), repos="http://cloud.r-project.org/")'

## Install mlm
RUN R -e 'devtools:::install_github("dgarrimar/mlm")'

## Install python3 packages
RUN pip3 install pandas scipy

## Install python2 packages
RUN pip install numpy scipy

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

## Install pyLMM
RUN git clone https://github.com/michaelbilow/pylmm_zarlab \
    && cd pylmm_zarlab \
    # Fix a couple of bugs in the pyLMM code
    && sed -i 's/True -/True ^/g' pylmm/lmm.py \
    && sed -i 's/print(self.N)//' pylmm/input.py \ 
    && python setup.py install

## Install BOLT-LMM
RUN wget --quiet https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4.tar.gz \
    && tar zxf BOLT-LMM_v2.3.4.tar.gz \
    && cd BOLT-LMM_v2.3.4 \
    && chmod +x bolt \
    && ln -s $PWD/bolt /usr/local/bin/
