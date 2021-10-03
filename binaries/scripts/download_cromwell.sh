#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt-get install -y wget --quiet

wget -P ../ https://github.com/broadinstitute/cromwell/releases/download/69/cromwell-69.jar
