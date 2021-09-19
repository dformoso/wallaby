FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

ENV bedops_version 2.4.39

RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    vim \
    wget \
    curl \
    libc6-dev \
    less \
    unzip

WORKDIR /usr/local/bin
RUN curl -SL https://github.com/bedops/bedops/archive/v${bedops_version}.tar.gz \
    > v${bedops_version}.tar.gz
RUN tar -xzvf v${bedops_version}.tar.gz
WORKDIR /usr/local/bin/bedops-${bedops_version}
RUN make all
RUN make install_all
RUN cp /usr/local/bin/bedops-${bedops_version}/bin/* /usr/local/bin/
ENV PATH /usr/local/bin/bedops-${bedops_version}/bin:$PATH

WORKDIR /tmp
CMD ["/bin/bash"]