import argparse
import logging
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
from rescal import rescal_als
from numpy import dot,zeros
import itertools
import random
import sys
from collections import defaultdict,Counter
import scipy.spatial.distance

def predict_rescal_als(T,rank,lambda_A,lambda_R):
	A, R, _, _, _ = rescal_als(
		T, rank, init='nvecs', conv=1e-4,
		lambda_A=lambda_A, lambda_R=lambda_R
	)
	n = A.shape[0]
	P = zeros((n, n, len(R)))
	for k in range(len(R)):
		P[:, :, k] = dot(A, dot(R[k], A.T))
	return P

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Make inferences on CancerMine using RESCAL')
	parser.add_argument('--inCancermine',required=True,type=str,help='Cancermine TSV to calculate inference scores for')
	parser.add_argument('--outCancermine',required=True,type=str,help='Cancermine TSV with inferences scores added')
	parser.add_argument('--rank',default=10,type=int,help='Parameter for RESCAL')
	parser.add_argument('--lambda_A',default=0.0,type=float,help='Parameter for RESCAL')
	parser.add_argument('--lambda_R',default=0.0,type=float,help='Parameter for RESCAL')
	args = parser.parse_args()

	knowledgebase = defaultdict(lambda : -1.0)

	roles,gene_ids,cancer_ids = Counter(),Counter(),Counter()

	geneID_2_name = {}
	cancerID_2_name = {}

	with open(args.inCancermine,'r') as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			lineDict = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			role = lineDict['role']
			gene_id = lineDict['gene_id']
			cancer_id = lineDict['cancer_id']
			gene_name = lineDict['gene_normalized']
			cancer_name = lineDict['cancer_normalized']
			predictprob = float(lineDict['predictprob'])

			roles[role] += 1
			gene_ids[gene_id] += 1
			cancer_ids[cancer_id] += 1

			geneID_2_name[gene_id] = gene_name
			cancerID_2_name[cancer_id] = cancer_name

			key = (role,gene_id,cancer_id)
			knowledgebase[key] = max(knowledgebase[key],predictprob)

	roles = sorted(list(roles))
	gene_ids = sorted(list(gene_ids))
	cancer_ids = sorted(list(cancer_ids))

	roleMap = { r:i for i,r in enumerate(roles) }
	geneIDMap = { r:i for i,r in enumerate(gene_ids) }
	cancerIDMap = { r:(len(geneIDMap)+i) for i,r in enumerate(cancer_ids) }

	training = { (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]):prob for (role,gene_id,cancer_id),prob in knowledgebase.items() }

	combinedSize = len(gene_ids) + len(cancer_ids)
	trainingTensor = [ lil_matrix((combinedSize,combinedSize)) for _ in roles ]
	for (roleIndex,geneIndex,cancerIndex),prob in training.items():
		trainingTensor[roleIndex][geneIndex,cancerIndex] = prob
			
	print("Running RESCAL with rank=%d, lambda_A=%f and lambda_R=%f" % (args.rank,args.lambda_A,args.lambda_R))
	preds = predict_rescal_als(trainingTensor, args.rank, args.lambda_A, args.lambda_R)
	
	with open(args.inCancermine,'r') as inF, open(args.outCancermine,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		outHeaders = headers + ['inference_score']
		outF.write("\t".join(outHeaders) + "\n")
		for line in inF:
			lineSplit = line.strip('\n').split('\t')
			lineDict = { h:v for h,v in zip(headers,lineSplit) }
			role = lineDict['role']
			gene_id = lineDict['gene_id']
			cancer_id = lineDict['cancer_id']
			predictprob = float(lineDict['predictprob'])

			inference_score = preds[geneIDMap[gene_id],cancerIDMap[cancer_id],roleMap[role]]

			lineSplit.append(str(inference_score))
			outF.write("\t".join(lineSplit) + "\n")

	
