FROM ubuntu:16.04

# Set non-interactive mode
ENV DEBIAN_FRONTEND noninteractive

ARG PACKAGE_VERSION=0.2.4
ARG BUILD_PACKAGES="git python3 python3-pkg-resources build-essential"

# Install dependencies
RUN apt update && \
    apt install --yes $BUILD_PACKAGES

WORKDIR /tmp 
RUN git clone https://github.com/rrwick/Porechop.git
WORKDIR /tmp/Porechop 
RUN git checkout v$PACKAGE_VERSION
RUN python3 setup.py install && \
    rm -rf /tmp/Porechop && \
    apt remove --purge --yes git build-essential && \
    apt autoremove --purge --yes
