FROM jupyter/r-notebook:latest

RUN python3 -m pip install papermill

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

ENV bedops_version 2.4.39

USER root

RUN apt-get update -y && apt-get install -y \
    build-essential \
    libnss-sss \
    vim \
    wget \
    curl \
    libc6-dev \
    less \
    unzip \
    libfontconfig1-dev

RUN cd /usr/local/bin && \
    curl -SL https://github.com/bedops/bedops/archive/v${bedops_version}.tar.gz \
    > v${bedops_version}.tar.gz && \
    tar -xzvf v${bedops_version}.tar.gz && \
    cd /usr/local/bin/bedops-${bedops_version} && make all && make install_all
RUN cp /usr/local/bin/bedops-${bedops_version}/bin/* /usr/local/bin/
ENV PATH /usr/local/bin/bedops-${bedops_version}/bin:$PATH

USER $NB_UID

COPY install_libraries.R .
RUN Rscript install_libraries.R