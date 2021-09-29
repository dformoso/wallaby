# Wallaby
A set of .wdl pipelines for finding Lateral Gene Transfer from a donor organism into a recipient.

## Overview
This repository is not complete yet. The intention is to have users be able to run both local CPU and or GPU pipelines, as well as running the pipeline in google Cloud Platform.

For running the pipelines locally, it is assumed that the local machine has a GPU, the drivers have been installed correctly, nvidia-docker is funcional and can access CUDA drivers from a docker container, on Ubuntu 18.04.

This project uses the following tools for its deployment:
- Docker
- Docker Hub
- WDL Pipelines
- Cromwell Runner
- Womtool

## Run your first local pipeline
- Clone this repository to a local directory
~~~
git clone https://github.com/dformoso/wallaby.git
~~~

- Change directory to the scripts and binaries folder, and run the following
~~~
./scripts_and_binaries/download_cromwell_and_womtool_binaries.sh
~~~

- Your test pipeline locally
~~~
./scripts_and_binaries/run_donor_recipient.sh
~~~

## About Me
Twitter:
> https://twitter.com/danielmartinezf

Linkedin:
> https://www.linkedin.com/in/danielmartinezformoso/

Email:
> daniel.martinez.formoso@gmail.com
