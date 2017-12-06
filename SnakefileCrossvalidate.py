import itertools
import numpy as np

featureChoices = "entityTypes,unigramsBetweenEntities,bigrams,dependencyPathEdges,dependencyPathEdgesNearEntities".split(',')

expectedFiles = []
for relationType in ["Driver","Oncogene","Tumor_Suppressor"]:
	for threshold in np.arange(0.1, 1.0, 0.01):
		#for featureCount in range(1,len(featureChoices)+1):
		#	for features in itertools.combinations(featureChoices,featureCount):
		filename = "crossvalidate/%s_%f_%s.txt" % (relationType,threshold,",".join(featureChoices))
		expectedFiles.append(filename)

rule:
	input: expectedFiles

rule:
	output: 'crossvalidate/{relationtype}_{threshold}_{features}.txt'
	shell: 'python crossvalidateModel.py --relationTypes {wildcards.relationtype} --inTrain cancermine_corpus/train/ --threshold {wildcards.threshold} --features {wildcards.features} --outFile {output}'

