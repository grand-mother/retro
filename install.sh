#!/bin/bash

echo "o Fetching dependencies ..."
git submodule update --init
echo "--> Done"
make
