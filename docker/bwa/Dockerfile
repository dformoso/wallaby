FROM ubuntu:18.04
    
ENV DEBIAN_FRONTEND noninteractive
ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

RUN apt-get update -y
RUN apt-get install -y bwa

WORKDIR /tmp
CMD ["/bin/bash"]