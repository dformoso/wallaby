FROM openjdk:13-alpine

USER root
WORKDIR /tmp

ENV TRIMMOMATIC_JAR_PATH /Trimmomatic-0.39
ENV TRIMMOMATIC_ADAPTERS_PATH /Trimmomatic-0.39/adapters

RUN apk add --no-cache wget zip && \
    wget -q http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    mv  Trimmomatic-0.39/trimmomatic-0.39.jar Trimmomatic-0.39/trimmomatic.jar && \
    rm Trimmomatic-0.39.zip && \
    apk del wget zip
    
WORKDIR /tmp
CMD ["/bin/bash"]