#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

# https://github.com/rocker-org/rocker-versioned

pwd=`pwd`
root_dir=`echo $pwd | rev | cut -d'/' -f4- | rev`

sudo docker stop rstudio || true && docker rm rstudio || true
sudo docker run -itd \
  --name rstudio \
  --hostname rstudio \
  --restart always \
  -e USER=rstudio \
  -e PASSWORD=rstudio \
  -p 8787:8787 \
  -v "${root_dir}:/home/rstudio" \
  dformoso/rstudio:latest 

# Access RStudio
  # http://localhost:8787

# Delete / connect to instance
  # docker container rm --force rstudio
  # docker exec -it rstudio bash

# Get logs
  # sudo docker logs --tail 50 --follow --timestamps rstudio
