FROM ubuntu:xenial

RUN apt-get update && apt-get install -y software-properties-common &&\
  add-apt-repository -y ppa:openjdk-r/ppa &&\
  apt-get update && apt-get install -y \
  openjdk-11-jre \
  wget \
  unzip

RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip &&\
  unzip Trimmomatic-0.39.zip &&\
  rm -rf Trimmomatic-0.39.zip &&\
  chmod +x Trimmomatic-0.39/trimmomatic-0.39.jar &&\
  echo "#!/bin/bash" >> /Trimmomatic-0.39/trimmomatic &&\
  echo "exec java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar """"$""@"""" " >> /Trimmomatic-0.39/trimmomatic &&\
  chmod +x /Trimmomatic-0.39/trimmomatic

ENV PATH="${PATH}:/Trimmomatic-0.39/"