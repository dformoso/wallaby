#!/bin/bash
# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../config/local_cpu.conf \
  -jar ../../binaries/cromwell-69.jar \
    run ../srrs_multi_donor_recipient.wdl \
    --inputs ../inputs/srrs/htlv1_rnaseq/inputs.json \
    --options ../inputs/srrs/htlv1_rnaseq/options.json