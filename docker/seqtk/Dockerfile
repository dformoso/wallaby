FROM ubuntu:18.04
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /seqtk

RUN apt-get update
RUN apt-get install git build-essential libz-dev -y
RUN git clone https://github.com/lh3/seqtk.git;
RUN cd seqtk; make
ENV PATH="/seqtk/seqtk:$PATH"

CMD ["/bin/bash"]


