#!/bin/bash

current_dir=$(pwd)
fit_model=$current_dir/cascade_fitting.R
analyze_results=$current_dir/analyze_results.py

k_conversion_const=1

while [ $k_conversion_const -lt 50 ]; do
	echo k_convert_gfp_to_tetR: $k_conversion_const

	Rscript $fit_model --k_convert_gfp_to_tetR=$k_conversion_const
	python3 $analyze_results --param_search --k_convert_gfp_to_tetR=$k_conversion_const

	# increment k
	let k_conversion_const=k_conversion_const+2
done
