#!/usr/bin/env bash

# Find the directory that this script is inside
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Trigger docker build, tagged rjbohlender/ibdmap:latest
docker build --platform linux/amd64 -t rjbohlender/ibdmap $SCRIPT_DIR/..
