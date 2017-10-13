#!/bin/bash

# Script root directory.
sim_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# Set the path for dynamic libraries.
lib_dir=$sim_dir/lib
[[ "$LD_LIBRARY_PATH" =~ "${lib_dir}" ]] || export LD_LIBRARY_PATH=${lib_dir}:$LD_LIBRARY_PATH

# Set the PATH for binaries.
bin_dir=$sim_dir/bin
[[ "$PATH" =~ "${bin_dir}" ]] || export PATH=${bin_dir}:$PATH

# Set the PYTHONPATH For python module.
python_dir=$sim_dir/lib/python
[[ "$PYTHONPATH" =~ "${python_dir}" ]] || export PYTHONPATH=${python_dir}:$PYTHONPATH

# Set the materials.
export PUMAS_MDF=$sim_dir/deps/danton/share/materials/materials.xml
export PUMAS_DEDX=$sim_dir/deps/danton/share/materials/dedx
