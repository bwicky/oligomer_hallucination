#!/bin/bash

rm -r rosetta_cmds rosetta_min_and_scored_models
mkdir rosetta_min_and_scored_models

for pdb in af2_models/*pdb 
	do
		basename=`basename $pdb`
		basename="${basename%.*}"
		pdb=`realpath $pdb`
		echo $pdb
		cmd_str="/software/rosetta/versions/v2022.12-dev61855/bin/rosetta_scripts.hdf5.linuxgccrelease @flags -s ${pdb} -parser:protocol ROSETTA_min_and_score.xml -out:path:all rosetta_min_and_scored_models/ -out:file:scorefile ${basename}.sc"
		echo $cmd_str >> rosetta_cmds

	done
