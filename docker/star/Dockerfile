FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update -y

RUN apt-get install -y \
 build-essential \
 gcc-multilib \
 apt-utils \
 zlib1g-dev

RUN apt-get install -y git

WORKDIR /usr/local/
RUN pwd
RUN git clone https://github.com/alexdobin/STAR.git

WORKDIR /usr/local/STAR/
RUN pwd
RUN git checkout 2.7.7a 

WORKDIR /usr/local/STAR/source 
RUN pwd
RUN make STAR

ENV PATH /usr/local/STAR/source:$PATH

WORKDIR /tmp
CMD ["/bin/bash"]