import sys
import itertools
import kindred
import pickle
import argparse
import codecs
import time
import re
import string
from collections import defaultdict,Counter
import json

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

def parseAndFindEntities(biocFile,filterTermsFile,wordlistPickle,outSentencesFilename):
	print("%s : start" % now())

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)

	with open(filterTermsFile,'r') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	timers = Counter()

	outSentences = []

	currentID = None
	duplicateCheck = set()

	print("%s : processing..." % now())
	parser = kindred.Parser()
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=True,detectMicroRNA=True,acronymDetectionForAmbiguity=True,mergeTerms=True)
	for corpusno,corpus in enumerate(kindred.iterLoadDataFromBioc(biocFile)):
		startTime = time.time()
		parser.parse(corpus)
		timers['parser'] += time.time() - startTime
		print("%s : parsed" % now())

		startTime = time.time()
		ner.annotate(corpus)
		timers['ner'] += time.time() - startTime
		print("%s : ner" % now())

		startTime = time.time()

		for doc in corpus.documents:

			# Reset the duplicate check set for each new PMID
			if doc.metadata['id'] != currentID:
				currentID = doc.metadata['id']
				duplicateCheck = set()

			for sentence in doc.sentences:
				sentenceTextLower = sentence.text.lower()
				containsFilterTerm = any( ft in sentenceTextLower for ft in filterTerms)
				if not containsFilterTerm:
					continue

				entityinSentence = {'cancer':False,'gene':False}

				entityTypesInSentence = sentence.entityIDToType.values()
				entityinSentence['cancer'] = 'cancer' in entityTypesInSentence
				entityinSentence['gene'] = 'gene' in entityTypesInSentence

				if entityinSentence['cancer'] and entityinSentence['gene']:
					sentenceText = sentence.text.strip(string.whitespace + ',')

					if not sentenceText in duplicateCheck:
						tmpData = dict(doc.metadata)
						tmpData['sentence'] = sentenceText
						outSentences.append(tmpData)
						duplicateCheck.add(sentenceText)

		timers['entitiesAdded'] += time.time() - startTime

		print("%s : entities added" % now())
		sys.stdout.flush()

	with open(outSentencesFilename,'w') as f:
		json.dump(outSentences,f,indent=2)

	print("%s : done" % now())
	
	for section,sectiontime in timers.items():
		print("%s\t%f" % (section,sectiontime))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	parser.add_argument('--filterTerms',required=True)
	parser.add_argument('--wordlistPickle',required=True)
	parser.add_argument('--outSentencesFilename',required=True)

	args = parser.parse_args()

	parseAndFindEntities(args.biocFile,args.filterTerms,args.wordlistPickle,args.outSentencesFilename)

