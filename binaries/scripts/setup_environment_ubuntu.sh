#!/bin/bash

# Requesting sudo access
[ "$UID" -eq 0 ] || exec sudo bash "$0" "$@"

sudo apt-get update -y && sudo apt-get install -y libmysql-java
sudo apt-get install -y wget --quiet
sudo apt-get install -y samtools
sudo apt-get install -y default-jre

sudo apt-get remove -y docker docker-engine docker.io containerd runc
sudo apt-get update -y
sudo apt-get install -y \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
sudo apt-get update -y
sudo apt-get install -y docker-ce docker-ce-cli containerd.io
