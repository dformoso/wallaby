#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"


# Create a docker container running an unauthenticated Jupyter instance.
# The instance will be able to operate docker containers in the host system if docker installed

pwd=`pwd`
root_dir=`echo $pwd | rev | cut -d'/' -f4- | rev`

sudo docker stop fileserver || true && sudo docker rm fileserver || true
sudo docker run -itd \
    --name fileserver \
    --hostname fileserver \
    --restart always \
    -p 8080:80 \
    -v "${root_dir}/wallaby/workflows/outputs:/usr/local/apache2/htdocs/" \
    dformoso/apache2:latest