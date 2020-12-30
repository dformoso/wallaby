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

RUN mkdir -p /usr/local/ncbi/sra-toolkit
WORKDIR /usr/local/ncbi/sra-toolkit
RUN curl -s https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.9/sratoolkit.2.10.9-ubuntu64.tar.gz | \
    tar zxf - --strip-components=1

ENV PATH="/usr/local/ncbi/sra-toolkit/bin:${PATH}"

RUN mkdir ~/.ncbi/
COPY user-settings.mkfg /root/.ncbi

WORKDIR /tmp
CMD ["/bin/bash"]
