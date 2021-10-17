#!/bin/bash
# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../cromwell_config/local_cpu.conf \
  -jar ../../binaries/cromwell-69.jar \
    run ../multi_donor_recipient.wdl \
    --inputs ../inputs/test/inputs.json \
    --options ../inputs/test/options.json