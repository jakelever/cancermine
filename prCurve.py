import argparse
import sys
import kindred
import random
import numpy as np
import sklearn.metrics
import os

if __name__ == '__main__':
	entityTypes = {}
	entityTypes['Driver'] = ('cancer','gene')
	entityTypes['Oncogene'] = ('cancer','gene')
	entityTypes['Tumor_Suppressor'] = ('cancer','gene')

	reltypes = ",".join(sorted(list(entityTypes.keys())))

	parser = argparse.ArgumentParser(description='Calculate classifier metrics for a particulate relation type and data set size')
	parser.add_argument('--trainDir',type=str,required=True,help='Directory containing training corpus to load (in standoff format)')
	parser.add_argument('--testDir',type=str,required=True,help='Directory containing testing corpus to load (in standoff format)')
	parser.add_argument('--reltype',type=str,required=True,help='Relation type to analyze. Must be one of %s' % reltypes)
	parser.add_argument('--outCurve',type=str,required=True,help='File to output curve data to')
	args = parser.parse_args()

	with open(args.outCurve,'w') as outF:
		outF.write("%s\t%s\t%s\n" % ('threshold','precision','recall'))
		for threshold in [-0.1] + list(np.arange(0,1,0.01)) + [1.0]:
			train = kindred.loadDir('standoff',args.trainDir)
			gold = kindred.loadDir('standoff',args.testDir)
			
			# Trim back to relation type of choice
			for doc in train.documents:
				doc.relations = [ r for r in doc.relations if r.relationType == args.reltype ]
			for doc in gold.documents:
				doc.relations = [ r for r in doc.relations if r.relationType == args.reltype ]

			entityType = entityTypes[args.reltype]
			entityCount = len(entityType)
			classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,entityCount=entityCount,acceptedEntityTypes=[entityType])

			classifier.train(train)

			predictions = gold.clone()
			predictions.removeRelations()

			classifier.predict(predictions)

			TP,FN,FP = 0,0,0
			for goldDoc,testDoc in zip(gold.documents,predictions.documents):
				goldTuples = [ (r.relationType,tuple(r.entityIDs)) for r in goldDoc.relations ]
				testTuples = [ (r.relationType,tuple(r.entityIDs)) for r in testDoc.relations ]

				totalSet = set(goldTuples + testTuples)
				for relation in totalSet:
					inGold = relation in goldTuples
					inTest = relation in testTuples

					if inGold and inTest:
						TP += 1
					elif inGold:
						FN += 1
					elif inTest:
						FP += 1

			if (TP+FP) != 0 and (TP+FN) != 0:
				precision = TP / float(TP+FP)
				recall = TP / float(TP+FN)

				outF.write("%f\t%f\t%f\n" % (threshold,precision,recall))
				print(threshold,precision,recall)
				sys.stdout.flush()

