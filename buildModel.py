import kindred
import argparse
import pickle
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Build and save a classifier')
	parser.add_argument('--inTrain',type=str,required=True)
	parser.add_argument('--outModel_Driver',type=str,required=True)
	parser.add_argument('--outModel_Oncogene',type=str,required=True)
	parser.add_argument('--outModel_TumorSuppressor',type=str,required=True)
	parser.add_argument('--conservativeThresholds',action='store_true')

	args = parser.parse_args()

	if args.conservativeThresholds:
		thresholds = {'Driver':0.80, 'Oncogene': 0.76, 'Tumor_Suppressor': 0.92}
	else:
		thresholds = {'Driver':0.5, 'Oncogene': 0.5, 'Tumor_Suppressor': 0.5}

	for relationType,outModel in zip(['Driver','Oncogene','Tumor_Suppressor'], [args.outModel_Driver,args.outModel_Oncogene,args.outModel_TumorSuppressor] ):
		print("Building %s model" % relationType)
		print("  Loading training")
		goldDir = 'gold'
		trainCorpus = kindred.loadDir(dataFormat='standoff',directory=args.inTrain)
	
		for doc in trainCorpus.documents:
			doc.relations = [ r for r in doc.relations if r.relationType == relationType ]

		print("  Doing training")
		features = "entityTypes,unigramsBetweenEntities,bigrams,dependencyPathEdges,dependencyPathEdgesNearEntities".split(',')
		threshold = thresholds[relationType]
		classifier = kindred.RelationClassifier(classifierType='LogisticRegression',threshold=threshold,features=features,acceptedEntityTypes=[('cancer','gene')])
		classifier.train(trainCorpus)

		print("  Saving classifer")
		with open(outModel,'wb') as f:
			pickle.dump(classifier,f)

		print("  Output done!")

