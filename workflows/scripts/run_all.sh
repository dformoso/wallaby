#!/bin/bash
# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

./run_hpv16_rnaseq.sh
./run_hpv18_rnaseq.sh