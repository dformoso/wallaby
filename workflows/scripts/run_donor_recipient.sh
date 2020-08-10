#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

# Need to install MYSQL connector
# sudo apt-get update -y && sudo apt-get install -y libmysql-java

java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../../config/local_cpu.conf \
  -jar ../../binaries/cromwell-52.jar \
    run ../donor_recipient.wdl \
    --inputs ../../inputs/donor_recipient.json \
    --options ../../config/options.json
    