#!/bin/bash

if [ $# -ne 2 ]; then
    echo "Usage: $0 <version [0, 1]> <data_size [0, 1, 2]>"
    echo "Example: $0 0 1"
    exit 1
fi

version=$1
data_size=$2

gcc -o solver_${version} serial_solver_${version}.c -lm
./solver_${version} < ../dataset/data/dataset_${version}_${data_size}.txt
echo -e "Timing results:\n"
cat timing_serial_solver_${version}.txt
