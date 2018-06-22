import argparse
import itertools
import random
import sys
import os
from collections import defaultdict,Counter

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Make inferences on CancerMine using RESCAL')
	parser.add_argument('--inCancermine',required=True,type=str,help='Cancermine TSV to calculate inference scores for')
	parser.add_argument('--outAnnotations',required=True,type=str,help='Cancermine TSV with inferences scores added')
	args = parser.parse_args()

	random.seed(1)

	sentences = defaultdict(list)
	maxpredictprobs = defaultdict(lambda : -1.0)
	inference_scores = {}

	geneID2Name,cancerID2Name = {},{}

	with open(args.inCancermine,'r') as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			lineDict = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			role = lineDict['role']
			gene_id = lineDict['gene_id']
			cancer_id = lineDict['cancer_id']
			gene_normalized = lineDict['gene_normalized']
			cancer_normalized = lineDict['cancer_normalized']
			gene_name = lineDict['gene_name']
			cancer_name = lineDict['cancer_name']
			sentence = lineDict['sentence']
			predictprob = float(lineDict['predictprob'])
			inference_score = float(lineDict['inference_score'])

			geneID2Name[gene_id] = gene_normalized
			cancerID2Name[cancer_id] = cancer_normalized

			key = (role,gene_id,cancer_id)

			sentences[key].append((gene_name,cancer_name,sentence))
			maxpredictprobs[key] = max(maxpredictprobs[key],predictprob)
			inference_scores[key] = inference_score	


	keys = sorted(list(sentences.keys()))

	existingAnnotations = set()
	if os.path.isfile(args.outAnnotations):
		with open(args.outAnnotations) as f:
			for line in f:
				yesOrNo,role,gene_id,cancer_id,sentenceCount,maxpredictprob,inference_score = line.strip().split('\t')
				existingAnnotations.add((role,gene_id,cancer_id))

	while True:
		annotationCount = len(existingAnnotations)

		key = random.choice(keys)
		role,gene_id,cancer_id = key

		if key in existingAnnotations:
			continue

		gene_normalized = geneID2Name[gene_id]
		cancer_normalized = cancerID2Name[cancer_id]


		print('#'*30 + ' (%d)' % annotationCount)

		print()

		for gene_name,cancer_name,sentence in sentences[key]:
			print("%s (%s,%s)" % (sentence,gene_name,cancer_name))
			print("-"*30)

		print()
		print('Proposed cancer gene role:')
		print("role:%s | gene:%s | cancer:%s" % (role,gene_normalized,cancer_normalized))
		print()
		while True:
			yesOrNo = input('Does this proposed role appear in the sentences shown? (y/n/-) ').lower()
			if yesOrNo in ['y','n','-']:
				break
		print()

		sentenceCount = len(sentences[key])
		maxpredictprob = maxpredictprobs[key]
		inference_score = inference_scores[key]

		existingAnnotations.add(key)
		with open(args.outAnnotations,'a+') as f:
			out = [yesOrNo,role,gene_id,cancer_id,sentenceCount,maxpredictprob,inference_score]
			f.write("\t".join(map(str,out)) + "\n")

		#break
	
