#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

# Need to install MYSQL connector
# sudo apt-get update -y && sudo apt-get install -y libmysql-java

java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../../cromwell_config/local_cpu.conf \
  -jar ../../binaries/cromwell-52.jar \
    run ../multi_donor_recipient.wdl \
    --inputs ../inputs/multi_donor_recipient_hpv16.json \
    --options ../../cromwell_config/options_local.json
    