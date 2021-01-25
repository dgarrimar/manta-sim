FROM rocker/r-ver:3.5.2

RUN apt-get update --fix-missing -qq && apt-get install -y -q \
# R:devtools dependencies
       libcurl4-openssl-dev \
       libssl1.0-dev \
       libssh2-1-dev \
       libxml2-dev \
       zlib1g-dev \
# libgit2 
       cmake \
       pkg-config \
       git \
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
# other
       wget \
       bc \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

## Install libgit2 (compile from source as in https://github.com/libgit2/pygit2/issues/999#issuecomment-614510556) 
RUN git clone --depth=1 -b v1.0.0 https://github.com/libgit2/libgit2.git ~/libgit2_src \
    && cd ~/libgit2_src \
    && cmake . -DBUILD_CLAR=OFF -DCMAKE_BUILD_TYPE=Release -DEMBED_SSH_PATH=~/libssh2_src -DCMAKE_INSTALL_PREFIX=~/libgit2 \
    && cmake --build . --target install \
    && cp -r ~/libgit2/* /usr/bin \
    && cp -r ~/libgit2/* /usr/local

## Install R packages
RUN R -e 'install.packages(c("devtools", "optparse", "data.table", "seqminer", "vegan", "BEDMatrix", "copula"), repos="http://cloud.r-project.org/")'

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
RUN wget --quiet https://github.com/genetics-statistics/GEMMA/releases/download/0.98.1/gemma-0.98.1-linux-static.gz \
    && gunzip gemma-0.98.1-linux-static.gz \
    && chmod +x gemma-0.98.1-linux-static \
    && mv gemma-0.98.1-linux-static /usr/local/bin/gemma

## Install BOLT-LMM
RUN wget --quiet https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.3.4.tar.gz \
    && tar zxf BOLT-LMM_v2.3.4.tar.gz \
    && cd BOLT-LMM_v2.3.4 \
    && chmod +x bolt \
    && ln -s $PWD/bolt /usr/local/bin/

## Install mlm
RUN R -e 'devtools:::install_github("dgarrimar/mlm")'

## Change ENTRYPOINT
ENTRYPOINT /bin/bash
