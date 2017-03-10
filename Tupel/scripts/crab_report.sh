#!/bin/sh
cwd=$PWD
mulumi_files=()
ellumi_files=()
era=(B_bk C D E F G H_v2 H_v3)
	len=${#era[*]}
	i=0
	while [ $i -lt $len ]; do
	crab report data_mc_grid/crab_Data13TeV_SingleMuon_2016${era[$i]}/
	crab report data_mc_grid/crab_Data13TeV_SingleElectron_2016${era[$i]}/
	mulumi_files+=("data_mc_grid/crab_Data13TeV_SingleMuon_2016"${era[$i]}"/results/processedLumis.json")
	ellumi_files+=("data_mc_grid/crab_Data13TeV_SingleElectron_2016"${era[$i]}"/results/processedLumis.json")
  	echo ${era[$i]}
	echo ${cwd}
  	let i++
   done
	printf 'mergeJSON.py '
	for j in "${ellumi_files[@]}"; do
	printf '%s\t' "$j"
   done
#exit 0
