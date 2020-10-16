FROM ubuntu:18.04

ENV DEBIAN_FRONTEND noninteractive

WORKDIR /tmp

# FASTQC
ENV VERSION 0.11.8
ENV NAME fastqc
ENV URL "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/${NAME}_v${VERSION}.zip"

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        wget \
        unzip \
        binutils \
        libfindbin-libs-perl \
        openjdk-11-jdk

RUN wget -q -O /opt/${NAME}_v${VERSION}.zip $URL && \
    cd /opt && \
    unzip ./${NAME}_v${VERSION}.zip &&\
    chmod 777 /opt/FastQC/${NAME} && \
    ln -s /opt/FastQC/${NAME} /usr/local/bin/ &&\
    rm /opt/${NAME}_v${VERSION}.zip

RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* /usr/share/doc/* && \
    strip /usr/local/bin/*; true

# MULTIQC
RUN apt-get update && apt-get install -y \
    software-properties-common
RUN add-apt-repository universe
RUN apt-get update && apt-get install -y python3-pip
RUN pip3 install multiqc

WORKDIR /tmp
CMD ["/bin/bash"]