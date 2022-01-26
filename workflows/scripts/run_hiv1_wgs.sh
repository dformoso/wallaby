#!/bin/bash
# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../config/local_cpu.conf \
  -jar ../../binaries/cromwell-69.jar \
    run ../fastqs_multi_donor_recipient.wdl \
    --inputs ../inputs/fastqs/hiv1_wgs/inputs.json \
    --options ../inputs/fastqs/hiv1_wgs/options.json