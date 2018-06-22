
RANK = list(range(50,500,10))
LAMBDA_A = [0,0.001,0.01,0.1,0.5]
LAMBDA_R = [0,0.001,0.01,0.1,0.5]

rule all:
	input: expand("results/{rank}_{lambda_A}_{lambda_R}", rank=RANK, lambda_A=LAMBDA_A, lambda_R=LAMBDA_R)

rule argh:
	output: "results/{rank}_{lambda_A}_{lambda_R}"
	shell: "python desperateAttempts.py --inCancermine /projects/bioracle/jake/pubrunner/workspace/CancerMine/full/cancermine_all.tsv --annotationData annotations.tsv --rank {wildcards.rank} --lambda_A {wildcards.lambda_A} --lambda_R {wildcards.lambda_R} > {output}"

