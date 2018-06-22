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
	parser.add_argument('--cancermine',required=True,type=str,help='Cancermine TSV to make inferences with')
	args = parser.parse_args()

	yearSplit = 2010

	#firstPub = defaultdict(lambda : 9999)
	#appearances = defaultdict(list)
	knowledgebase = set()

	roles,gene_ids,cancer_ids = Counter(),Counter(),Counter()

	geneID_2_name = {}
	cancerID_2_name = {}

	with open(args.cancermine,'r') as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			lineDict = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			role = lineDict['role']
			gene_id = lineDict['gene_id']
			cancer_id = lineDict['cancer_id']
			gene_name = lineDict['gene_normalized']
			cancer_name = lineDict['cancer_normalized']
			year = lineDict['year']
			if "mir" in gene_name.lower():
				continue

			roles[role] += 1
			gene_ids[gene_id] += 1
			cancer_ids[cancer_id] += 1

			geneID_2_name[gene_id] = gene_name
			cancerID_2_name[cancer_id] = cancer_name

			key = (role,gene_id,cancer_id)
			#firstPub[key] = min(year,firstPub[key])
			knowledgebase.add(key)

	cutoff  = 5
	#print(roles)
	roles = set([ x for x,c in roles.items() if c > cutoff ])
	gene_ids = set([ x for x,c in gene_ids.items() if c > cutoff ])
	cancer_ids = set([ x for x,c in cancer_ids.items() if c > cutoff ])

	knowledgebase = [ (role,gene_id,cancer_id) for role,gene_id,cancer_id in knowledgebase if role in roles and gene_id in gene_ids and cancer_id in cancer_ids ]
	knowledgebase = sorted(knowledgebase)

	random.seed(1)
	random.shuffle(knowledgebase)

	roles = sorted(list(roles))
	gene_ids = sorted(list(gene_ids))
	cancer_ids = sorted(list(cancer_ids))

	roleMap = { r:i for i,r in enumerate(roles) }
	geneIDMap = { r:i for i,r in enumerate(gene_ids) }
	cancerIDMap = { r:(len(geneIDMap)+i) for i,r in enumerate(cancer_ids) }

	roleInverseMap = { i:r for r,i in roleMap.items() }
	geneIDInverseMap = { i:r for r,i in geneIDMap.items() }
	cancerIDInverseMap = { i:r for r,i in cancerIDMap.items() }

	training = knowledgebase[:int(round(0.75*len(knowledgebase)))]
	testing = knowledgebase[int(round(0.75*len(knowledgebase))):]

	trainingTensor = [ (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]) for role,gene_id,cancer_id in training ]
	testingTensor = [ (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]) for role,gene_id,cancer_id in testing ]

	#print(len(roles),len(gene_ids),len(cancer_ids))

	combinedSize = len(gene_ids) + len(cancer_ids)

	combinedX = [ lil_matrix((combinedSize,combinedSize)) for _ in roles ]

	trainingX = [ lil_matrix((combinedSize,combinedSize)) for _ in roles ]
	testingX = [ lil_matrix((combinedSize,combinedSize)) for _ in roles ]
	for a,b,c in trainingTensor:
		trainingX[a][b,c] = 1
		combinedX[a][b,c] = 1
	for a,b,c in testingTensor:
		testingX[a][b,c] = 1
		combinedX[a][b,c] = 1

	rank = 10
	lambda_A = 0
	lambda_R = 0
	#preds = predict_rescal_als(trainingX, rank, lambda_A, lambda_R)

	#T = trainingX
	T = combinedX

	A, R, _, _, _ = rescal_als(
		T, rank, init='nvecs', conv=1e-4,
		lambda_A=lambda_A, lambda_R=lambda_R
	)
	n = A.shape[0]
	preds = zeros((n, n, len(R)))
	for k in range(len(R)):
		preds[:, :, k] = dot(A, dot(R[k], A.T))

	#cancerDists = {}
	import sklearn.preprocessing as pp
	def cosine_similarities(mat):
		col_normed_mat = pp.normalize(mat, axis=0)
		return col_normed_mat.T * col_normed_mat

	import sklearn.metrics.pairwise

	roleID = roleMap['Tumor_Suppressor']

	#m = preds[:len(gene_ids),-len(cancer_ids):,roleID]
	m = combinedX[roleID][:len(gene_ids),-len(cancer_ids):]
	#print(len(gene_ids),len(cancer_ids))
	#print(m.shape)
	cancerSimilarities = sklearn.metrics.pairwise.cosine_similarity(m.T)
	#print(cancerSimilarities.shape)
	geneSimilarities = sklearn.metrics.pairwise.cosine_similarity(m)
	#print(geneSimilarities.shape)

	#for (i,cancer1),(j,cancer2) in itertools.product(enumerate(cancer_ids),enumerate(cancer_ids)):
	#	cancerName1 = cancerID_2_name[cancer1]
	#	cancerName2 = cancerID_2_name[cancer2]
	#	out = [sim[i,j],cancer1,cancerName1,cancer2,cancerName2]
		#print("\t".join(map(str,out)))
	#sys.exit(0)

	#for cancer1,cID1 in cancerIDMap.items():
	cancer1 = 'DOID:684'
	cID1 = cancerIDMap[cancer1]
	#for cancer1,cID1 in [cancer
	if False:
		cancerName1 = cancerID_2_name[cancer1]
		v1 = preds[cID1,:,roleID]
		for cancer2,cID2 in cancerIDMap.items():
			cancerName2 = cancerID_2_name[cancer2]
			v2 = preds[cID2,:,roleID]
			print(v1.shape,v2.shape,v1.dot(v2))
			cosineDist = scipy.spatial.distance.cosine(v1,v2)
			out = [cosineDist,cancer1,cancerName1,cancer2,cancerName2]
			print("\t".join(map(str,out)))
			break
			
			#cancerDists[

	#sys.exit(0)

	preds[preds < 0.5] = 0
	aboveThresholdIndices = list(zip(*preds.nonzero()))

	#for a,b,c in testingTensor:
	for b,c,a in aboveThresholdIndices:
		score = preds[b,c,a]
		#if score < 0.9:
		#	continue
		#if trainingX[a][b,c] == 1 or testingX[a][b,c] == 1:
		if combinedX[a][b,c] == 1:
			continue

		role = roleInverseMap[a]
		gene_id = geneIDInverseMap[b]
		cancer_id = cancerIDInverseMap[c]
		gene_name = geneID_2_name[gene_id]
		cancer_name = cancerID_2_name[cancer_id]
		out = [score,role,gene_id,gene_name,cancer_id,cancer_name]

		if gene_name != 'APC':
			continue
		print("-"*30)
		print("\t".join(map(str,out)))
		#continue

		c2 = c - len(gene_ids)
		similarCancers = [ (cancerSimilarities[c2,j],cancer2) for j,cancer2 in enumerate(cancer_ids) ]
		similarCancers = [ (score,cancer2) for score,cancer2 in similarCancers if score > 0.3 ]
		similarCancers = sorted(similarCancers,reverse=True)
		#print(similarCancers[:5])
		#print(sim[
		similarGenes = [ (geneSimilarities[b,j],gene2) for j,gene2 in enumerate(gene_ids) ]
		similarGenes = [ (score,gene2) for score,gene2 in similarGenes if score > 0.3 ]
		similarGenes = sorted(similarGenes,reverse=True)
		#print(similarGenes[:5])
		#for score2,gene2 in similarGenes:
	#		print('huh?',score2,gene2)

		#break
		#continue
		for cancerScore,cancer2 in similarCancers:#[:5]:
			for geneScore,gene2 in similarGenes[:1]:
				gene2_index = geneIDMap[gene2]
				cancer2_index = cancerIDMap[cancer2]
				gene_name2 = geneID_2_name[gene2]
				cancer_name2 = cancerID_2_name[cancer2]
				if combinedX[a][gene2_index,cancer2_index] == 1:
					print('Explained!',cancerScore+geneScore,geneScore,gene2,gene_name2,cancerScore,cancer2,cancer_name2)
			


		#break

