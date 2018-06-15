#!/bin/bash

# Fetch dependencies
if [ -z "$(ls -A deps/danton)" ]; then
    echo "Fetching dependencies..."
    git submodule update --init
fi

# Bootstrap ninja
if [ ! -e deps/ninja/ninja ]; then
    cd deps/ninja
    ./configure.py --bootstrap
    cd ../..
fi

# Build the binaries"
./deps/ninja/ninja
