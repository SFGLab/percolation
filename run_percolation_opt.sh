#!/usr/bin/env bash

data_name="GM12878.CTCF.clusters_all_edges"
output_dir="./results/${data_name}/"

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

chromosome=$1

if ! [[ -d "${output_dir}" ]]; then
	mkdir -p ${output_dir}
fi

python percolation_opt.py "./data/${data_name}.csv_split/${chromosome}" "${output_dir}" ${chromosome} \
    -n 100 \
	--ncores 0 --batch_size 20