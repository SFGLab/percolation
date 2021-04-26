#!/usr/bin/env bash

for i in {1..22} X
do
	./run_percolation_opt.sh "chr${i}"
done