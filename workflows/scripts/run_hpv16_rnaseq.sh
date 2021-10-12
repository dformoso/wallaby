#!/bin/bash
# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

if [ -d "../../workflows/cromwell-final-outputs" ]; then
rm -Rf "../../workflows/cromwell-final-outputs"
fi
java \
  -DLOG_MODE=pretty \
  -Dconfig.file=../../cromwell_config/local_cpu.conf \
  -jar ../../binaries/cromwell-69.jar \
    run ../multi_donor_recipient.wdl \
    --inputs ../inputs/hpv16_inputs.json \
    --options ../../cromwell_config/options_local.json
if [ -d "../../workflows/cromwell-final-outputs" ]; then
mv "../../workflows/cromwell-final-outputs" "../../workflows/outputs-hpv16-rnaseq"
fi