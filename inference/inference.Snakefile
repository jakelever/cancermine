
rule infer:
	input:
		training="training.tensor",
		eval="validation.tensor"
	output:
		"data/scores_{rank}_{lambda_A}_{lambda_R}"
	shell:
		"python optimiseInferences.py --trainingTensor {input.training} --evalTensor {input.eval} --outScores {output} --rank {wildcards.rank} --lambda_A {wildcards.lambda_A} --lambda_R {wildcards.lambda_R}"

rule evaluate:
	input:
		"data/scores_{rank}_{lambda_A}_{lambda_R}"
	output:
		"data/auprc_{rank}_{lambda_A}_{lambda_R}"
	shell:
		"python evaluate.py --data {input} > {output}"

