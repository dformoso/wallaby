#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

# https://jupyter-docker-stacks.readthedocs.io/en/latest/using/selecting.html#jupyter-datascience-notebook

# Create a docker container running an unauthenticated Jupyter instance.
# The instance will be able to operate docker containers 
# in the host system if docker is installed

##### !!!! REPLACE YOUR LOCAL DIRECTORY BELOW !!!! #####

sudo docker stop jupyter || true && docker rm jupyter || true
sudo docker run -itd \
  --name jupyter \
  --hostname jupyter \
  --restart always \
  -e GRANT_SUDO=yes \
  -e CHOWN_HOME=yes \
  --user root \
  -p 8888:8888 \
  -v "/root/research/genomics/projects/wallaby:/home/jovyan/wallaby" \
  dformoso/jupyterlab:latest \
  start-notebook.sh --NotebookApp.token='' 


jupyter/r-notebook
dformoso/jupyterlab


# Access Jupyter
  # http://localhost:8888/lab

# Delete / connect to instance
  # docker container rm --force jupyter
  # docker exec -it jupyter bash

# Get logs
  # sudo docker logs --tail 50 --follow --timestamps jupyter