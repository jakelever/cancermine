import kindred
import argparse
import pickle
import os
from collections import Counter

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Count candidate relations in corpus and different annotation counts')
	parser.add_argument('--inTrain',type=str,required=True)

	args = parser.parse_args()

	corpus = kindred.loadDir(dataFormat='standoff',directory=args.inTrain)

	parser = kindred.Parser()
	parser.parse(corpus)

	candidateBuilder = kindred.CandidateBuilder(acceptedEntityTypes=[('cancer','gene')])
	candidateRelations = candidateBuilder.build(corpus)

	counter = Counter()
	for cr in candidateRelations:
		#assert len(cr.knownTypesAndArgNames) <= 1
		if len(cr.knownTypesAndArgNames) > 0:
			#knownType,argNames = cr.knownTypesAndArgNames[0]
			for knownType,argNames in cr.knownTypesAndArgNames:
				counter[knownType] += 1
		else:
			counter['None'] += 1

	print(counter)
	print('len(candidateRelations)=',len(candidateRelations))
