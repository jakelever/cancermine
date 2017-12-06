import kindred
import argparse
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Build and save a classifier')
	parser.add_argument('--inTrain',type=str,required=True)
	parser.add_argument('--threshold',type=float,required=True)
	parser.add_argument('--features',type=str,required=True)
	parser.add_argument('--relationTypes',type=str,required=True)
	parser.add_argument('--outFile',type=str,required=True)
	#parser.add_argument('--outModel',type=str,required=True)

	args = parser.parse_args()

	print("Loading training")
	corpus = kindred.loadDir(dataFormat='standoff',directory=args.inTrain)
	#parser = kindred.Parser()

	acceptedRelationTypes = set(args.relationTypes.split(','))
	for doc in corpus.documents:
		doc.relations = [ r for r in doc.relations if r.relationType in acceptedRelationTypes ]

	f01scores,precisions,recalls = [],[],[]
	for trainCorpus,validateCorpusGold in corpus.nfold_split(5):
		print("Doing training")

		validateCorpus = validateCorpusGold.clone()
		validateCorpus.removeRelations()

		features = args.features.split(',')
		classifier = kindred.RelationClassifier(acceptedEntityPairs=[('cancer','gene')],classifierType='LogisticRegression',threshold=args.threshold,features=features)
		classifier.train(trainCorpus)

		classifier.predict(validateCorpus)

		precision,recall,f1score = kindred.evaluate(validateCorpusGold,validateCorpus,metric='all',display=False)

		if (precision+recall) != 0.0:
			beta = 0.1
			f01score = (1+beta*beta) * (precision*recall) / (beta*beta*precision + recall)
		else:
			f01score = 0.0

		f01scores.append(f01score)
		precisions.append(precision)
		recalls.append(recall)
		#break

	txt_f01scores = ",".join(map(str,f01scores))
	txt_precisions = ",".join(map(str,precisions))
	txt_recalls = ",".join(map(str,recalls))

	mean_f01score = sum(f01scores) / float(len(f01scores))
	mean_precision = sum(precisions) / float(len(precisions))
	mean_recall = sum(recalls) / float(len(recalls))

	#print(f01scores)
	outData = [str(mean_f01score),str(mean_precision),str(mean_recall),str(args.threshold),args.features,txt_f01scores,txt_precisions,txt_recalls]
	outLine = "\t".join(outData)
	with open(args.outFile,'w') as f:
		print(outLine)
		f.write(outLine + '\n')

