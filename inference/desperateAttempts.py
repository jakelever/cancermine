import argparse
import random

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_curve, auc
from collections import Counter,defaultdict
import logging
from scipy.io.matlab import loadmat
from scipy.sparse import lil_matrix
from rescal import rescal_als
from numpy import dot,zeros
import itertools
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
	parser = argparse.ArgumentParser(description='Fit a linear model for final Cancermine results')
	parser.add_argument('--annotationData',type=str,required=True,help='TSV of data from annotation solicitation')
	parser.add_argument('--inCancermine',required=True,type=str,help='Cancermine TSV to calculate inference scores for')
	parser.add_argument('--rank',default=10,type=int,help='Parameter for RESCAL')
	parser.add_argument('--lambda_A',default=0.0,type=float,help='Parameter for RESCAL')
	parser.add_argument('--lambda_R',default=0.0,type=float,help='Parameter for RESCAL')
	args = parser.parse_args()

	knowledgebase = defaultdict(lambda : -1.0)

	roles,gene_ids,cancer_ids = Counter(),Counter(),Counter()

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
	

	X,y = [],[]

	with open(args.annotationData) as f:
		for line in f:
			yesOrNo,role,gene_id,cancer_id,sentenceCount,maxpredictprob,_ = line.strip().split('\t')
			sentenceCount = float(sentenceCount)
			maxpredictprob = float(maxpredictprob)

			if gene_id in geneIDMap and cancer_id in cancerIDMap and role in roleMap:
				inference_score = preds[geneIDMap[gene_id],cancerIDMap[cancer_id],roleMap[role]]
			else:
				inference_score = 0
			if yesOrNo in ['y','n']:
				y.append(1 if yesOrNo == 'y' else 0)
				#X.append([maxpredictprob])
				#X.append([sentenceCount,maxpredictprob])#,inference_score])
				X.append([sentenceCount,maxpredictprob,inference_score])
				#X.append([sentenceCount,maxpredictprob])

	X = np.array(X)
	#y = np.array(y)

	random.seed(1)

	X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=1)

	#clf = LogisticRegression()
	clf = LinearSVC()
	clf.fit(X_train,y_train)
	y_score = clf.decision_function(X_test)
	fpr, tpr, _ = roc_curve(y_test, y_score)
	roc_auc1 = auc(fpr, tpr)

	clf = LinearSVC()
	clf.fit(X_train[:,:2],y_train)
	y_score = clf.decision_function(X_test[:,:2])
	fpr, tpr, _ = roc_curve(y_test, y_score)
	roc_auc2 = auc(fpr, tpr)

	betterOrWorse = 'BETTER' if roc_auc1 > roc_auc2 else 'WORSE'

	print('RESULT\t%d\t%f\t%f\t|\t%f\t%f\t%s' % (args.rank, args.lambda_A, args.lambda_R, roc_auc1, roc_auc2, betterOrWorse))


