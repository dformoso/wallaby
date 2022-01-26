FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

USER root
WORKDIR /tmp

RUN apt-get -y -m update && \
    apt-get install -y \
    libidn11 fuse libxml2-dev \
    curl && \
    rm -rf /var/lib/apt/lists/*

RUN apt-get clean

# Install sra-toolkit
RUN mkdir -p /usr/local/ncbi/sra-toolkit
WORKDIR /usr/local/ncbi/sra-toolkit
RUN curl -s https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz | \
    tar zxf - --strip-components=1

ENV PATH="/usr/local/ncbi/sra-toolkit/bin:${PATH}"

RUN mkdir ~/.ncbi/
COPY user-settings.mkfg /root/.ncbi

# Install seqtk
WORKDIR /seqtk

RUN apt-get update
RUN apt-get install git build-essential libz-dev -y
RUN git clone https://github.com/lh3/seqtk.git;
RUN cd seqtk; make
ENV PATH="/seqtk/seqtk:$PATH"

WORKDIR /tmp

CMD ["/bin/bash"]
