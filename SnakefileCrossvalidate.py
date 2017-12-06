import itertools
import numpy as np

featureChoices = "entityTypes,unigramsBetweenEntities,bigrams,dependencyPathEdges,dependencyPathEdgesNearEntities".split(',')

expectedFiles = []
for threshold in np.arange(0.1, 1.0, 0.1):
	for featureCount in range(1,len(featureChoices)+1):
		for features in itertools.combinations(featureChoices,featureCount):
			filename = "crossvalidate/%f_%s.txt" % (threshold,",".join(features))
			expectedFiles.append(filename)

rule:
	input: expectedFiles

rule:
	output: 'crossvalidate/{threshold}_{features}.txt'
	shell: 'python crossvalidateModel.py --inTrain cancermine_corpus/train/ --threshold {wildcards.threshold} --features {wildcards.features} --outFile {output}'
