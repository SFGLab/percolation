#!/usr/bin/env bash

if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters"
    exit 1
fi

data_name=$1
output_dir="${data_name}_split"

if ! [[ -d "${output_dir}" ]]; then
	mkdir -p ${output_dir}
fi

for i in {1..22} X M
do
	grep "chr${i}" ${data_name} > ${output_dir}/"chr${i}"
done

echo "Number of lines:"
wc -l ${output_dir}/*