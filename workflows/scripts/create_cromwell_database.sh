#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

# https://github.com/broadinstitute/cromwell/blob/develop/scripts/docker-compose-mysql/compose/mysql/init/init_user.sql

# Create a docker container running an unauthenticated Jupyter instance.
# The instance will be able to operate docker containers in the host system if docker installed
docker stop mysql || true && docker rm mysql || true
sudo docker run -itd \
    --name mysql \
    --hostname mysql \
    --restart always \
    -e MYSQL_ROOT_PASSWORD=cromwell \
    -e MYSQL_DATABASE=cromwell_db \
    -e MYSQL_USER=cromwell \
    -e MYSQL_PASSWORD=cromwell \
    -p 3306:3306 \
    mysql/mysql-server:5.7


# Delete / connect to instance
  # docker container rm --force mysql
  # docker exec -it mysql bash

# Get logs
  # sudo docker logs --tail 50 --follow --timestamps mysql
  
# Get status
  # docker container list

