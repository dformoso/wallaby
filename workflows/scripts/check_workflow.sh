#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt-get install -y --quiet -qq graphviz

echo "Please enter a WDL filename in the pipelines folder, without the .wdl extension. Example \"fast5_to_bam\":"
read filename

java -jar ../../binaries/womtool-49.jar validate ../$filename.wdl
java -jar ../../binaries/womtool-49.jar graph ../$filename.wdl > ../$filename.dot
dot ../$filename.dot -Tpng -o ../$filename.png

