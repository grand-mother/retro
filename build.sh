#!/bin/bash

#  Fetch dependencies"
git submodule update --init

# Bootstrap ninja
if [ ! -e deps/ninja/ninja ]; then
    cd deps/ninja ; ./configure.py --bootstrap ; cd ../..
fi

# Build the binaries"
./deps/ninja/ninja

