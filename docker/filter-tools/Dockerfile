FROM broadinstitute/picard:2.22.9

USER root
WORKDIR /tmp

RUN apt-get update --allow-releaseinfo-change -y && apt-get install -y \
    build-essential \
    bwa \
    git \
    unzip \
    wget

WORKDIR /home
RUN echo "installing prinseq"
RUN wget http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz \
    -O /home/prinseq-lite-0.20.4.tar.gz \
    && tar -xvf prinseq-lite-0.20.4.tar.gz \
    && install -v prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq \
    && rm -rf prinseq-lite-0.20.4*

RUN apt-get update --allow-releaseinfo-change && \
    apt-get install -y samtools bedtools && \
    apt-get clean

WORKDIR /tmp
CMD ["/bin/bash"]
