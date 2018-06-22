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

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Make inferences on CancerMine using RESCAL')
	parser.add_argument('--cancermine',required=True,type=str,help='Cancermine TSV to make inferences with')
	parser.add_argument('--trainingTensor',required=True,type=str,help='Tensor with training data')
	parser.add_argument('--validationTensor',required=True,type=str,help='Tensor with validation data (with added negative data)')
	parser.add_argument('--testingTensor',required=True,type=str,help='Tensor with training data (with added negative data)')
	args = parser.parse_args()

	knowledgebase = defaultdict(lambda : -1.0)

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
			#year = lineDict['year']
			predictprob = float(lineDict['predictprob'])

			#if predictprob < 0.97:
			#	continue

			roles[role] += 1
			gene_ids[gene_id] += 1
			cancer_ids[cancer_id] += 1

			geneID_2_name[gene_id] = gene_name
			cancerID_2_name[cancer_id] = cancer_name

			key = (role,gene_id,cancer_id)
			#knowledgebase.add(key)
			knowledgebase[key] = max(knowledgebase[key],predictprob)

	roles = sorted(list(roles))
	gene_ids = sorted(list(gene_ids))
	cancer_ids = sorted(list(cancer_ids))

	roleMap = { r:i for i,r in enumerate(roles) }
	geneIDMap = { r:i for i,r in enumerate(gene_ids) }
	cancerIDMap = { r:(len(geneIDMap)+i) for i,r in enumerate(cancer_ids) }

	keys = sorted(list(knowledgebase.keys()))
	random.seed(1)
	random.shuffle(keys)

	firstHalf = keys[:int(round(0.5*len(keys)))]
	secondHalf = keys[int(round(0.5*len(keys))):]

	trainingKeys = firstHalf
	validationKeys = secondHalf[:int(round(0.5*len(secondHalf)))]
	testingKeys = secondHalf[int(round(0.5*len(secondHalf))):]

	training = { (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]):knowledgebase[(role,gene_id,cancer_id)] for role,gene_id,cancer_id in trainingKeys }
	positiveValidation = { (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]) for role,gene_id,cancer_id in validationKeys }
	positiveTesting = { (roleMap[role],geneIDMap[gene_id],cancerIDMap[cancer_id]) for role,gene_id,cancer_id in testingKeys }

	print("Training data # = ", len(training))
	print("Positive validation data # = ", len(positiveValidation))
	print("Positive testing data # = ", len(positiveTesting))

	roleChoices = sorted(list(roleMap.values()))
	geneChoices = sorted(list(geneIDMap.values()))
	cancerChoices = sorted(list(cancerIDMap.values()))

	classBalance = len(training) / float(len(roleChoices)*len(geneChoices)*len(cancerChoices))
	print("Training class balance: %f" % classBalance)
	
	negativeTesting = set()
	negativeValidation = set()

	while len(negativeValidation) < len(positiveValidation):
		roleIndex = random.choice(roleChoices)
		geneIndex = random.choice(geneChoices)
		cancerIndex = random.choice(cancerChoices)
		key = (roleIndex,geneIndex,cancerIndex)
		if key in training:
			continue
		if key in positiveValidation or key in negativeValidation:
			continue
		if key in positiveTesting or key in negativeTesting:
			continue
		negativeValidation.add(key)
	print("Built negative validation data", len(negativeValidation))

	while len(negativeTesting) < len(positiveTesting):
		roleIndex = random.choice(roleChoices)
		geneIndex = random.choice(geneChoices)
		cancerIndex = random.choice(cancerChoices)
		key = (roleIndex,geneIndex,cancerIndex)
		if key in training:
			continue
		if key in positiveValidation or key in negativeValidation:
			continue
		if key in positiveTesting or key in negativeTesting:
			continue
		negativeTesting.add(key)
	print("Built negative testing data", len(positiveTesting))

	combinedSize = len(geneChoices) + len(cancerChoices)
	with open(args.trainingTensor,'w') as outF:
		outF.write('# %d %d %d\n' % (len(roleChoices),combinedSize,combinedSize))
		for (roleIndex,geneIndex,cancerIndex),prob in training.items():
			outF.write("%d\t%d\t%d\t%f\n" % (roleIndex,geneIndex,cancerIndex,prob))
	with open(args.validationTensor,'w') as outF:
		outF.write('# %d %d %d\n' % (len(roleChoices),combinedSize,combinedSize))
		for roleIndex,geneIndex,cancerIndex in positiveValidation:
			outF.write("%d\t%d\t%d\t%d\n" % (roleIndex,geneIndex,cancerIndex,1))
		for roleIndex,geneIndex,cancerIndex in negativeValidation:
			outF.write("%d\t%d\t%d\t%d\n" % (roleIndex,geneIndex,cancerIndex,0))
	with open(args.testingTensor,'w') as outF:
		outF.write('# %d %d %d\n' % (len(roleChoices),combinedSize,combinedSize))
		for roleIndex,geneIndex,cancerIndex in positiveTesting:
			outF.write("%d\t%d\t%d\t%d\n" % (roleIndex,geneIndex,cancerIndex,1))
		for roleIndex,geneIndex,cancerIndex in negativeTesting:
			outF.write("%d\t%d\t%d\t%d\n" % (roleIndex,geneIndex,cancerIndex,0))


