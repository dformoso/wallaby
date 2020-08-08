FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

USER root
WORKDIR /tmp

RUN apt-get update && \
    apt-get install -y samtools && \
    apt-get clean

WORKDIR /tmp
CMD ["/bin/bash"]
