FROM buildpack-deps:stretch

RUN git clone https://github.com/lh3/minimap2
RUN cd minimap2 && make
RUN cp /minimap2/minimap2 /usr/local/bin

RUN apt-get update && apt-get install -y samtools