FROM ubuntu:18.04
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update -y && apt-get install -y \
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

WORKDIR /tmp
CMD ["/bin/bash"]
