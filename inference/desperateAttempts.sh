#!/bin/bash
set -eux

for rank in `seq 5 100`
do
	for lambda_A in 0 0.01 0.1
	do
		for lambda_R in 0 0.01 0.1
		do
			python desperateAttempts.py --inCancermine /projects/bioracle/jake/pubrunner/workspace/CancerMine/full/cancermine_all.tsv --annotationData annotations.tsv  --rank $rank --lambda_A $lambda_A --lambda_R $lambda_R | tee results/$rank\_$lambda_A\_$lambda_R
		done
	done
done
