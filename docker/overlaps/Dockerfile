FROM r-base:4.1.2

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update -y && apt-get install -y \
    autoconf automake make gcc perl zlib1g-dev \
    libbz2-dev liblzma-dev libcurl4-gnutls-dev \
    libssl-dev libncurses5-dev python3-pip \
    build-essential libnss-sss vim wget curl \
    libxml2-dev

WORKDIR /tmp
RUN wget https://github.com/samtools/samtools/releases/download/1.14/samtools-1.14.tar.bz2
RUN tar xvjf samtools-1.14.tar.bz2
WORKDIR /tmp/samtools-1.14
RUN make
ENV PATH /tmp/samtools-1.14:$PATH

RUN R -e 'install.packages("pacman", repos="http://cran.us.r-project.org")'
RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'install.packages("XML", repos="http://cran.us.r-project.org")'
RUN R -e 'library("pacman")'
RUN R -e 'library("BiocManager")'
RUN R -e 'library("XML")'
RUN R -e 'pacman::p_load(stringr, data.table)'
RUN R -e 'pacman::p_load(rtracklayer, TxDb.Hsapiens.UCSC.hg38.knownGene, Repitools)'
RUN R -e 'pacman::p_load(Organism.dplyr, GenomicFeatures, org.Hs.eg.db)'

RUN apt-get install -y python3.6
RUN pip3 install pandas

WORKDIR /tmp
COPY analyze_srr.r .
COPY filter_srr.py .