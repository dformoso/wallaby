#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt-get install -y wget --quiet

wget -P ../ https://github.com/broadinstitute/cromwell/releases/download/52/cromwell-52.jar
wget -P ../ https://github.com/broadinstitute/cromwell/releases/download/52/womtool-52.jar