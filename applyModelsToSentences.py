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

def getStandardizedTerm(text,externalID,IDToTerm):
	standardizedTerms = [ IDToTerm[eid] for eid in externalID.split(';') ]
	standardizedTerms = sorted(list(set(standardizedTerms)))

	standardizedTermsLower = [ st.lower() for st in standardizedTerms ]
	textLower = text.lower()

	if textLower in standardizedTermsLower:
		index = standardizedTermsLower.index(textLower)
		standardizedTerm = standardizedTerms[index]
	else:
		standardizedTerm = ";".join(standardizedTerms)
	return standardizedTerm

def standardizeMIRName(externalID):
	assert externalID.startswith('mirna|'), "Unexpected ID: %s" % externalID
	standardName = externalID[4:]

	search = re.search('mirna\|\D*(?P<id>\d+[A-Za-z]*)',externalID)
	if search:
		mirID = search.groupdict()['id']
		if not mirID is None:
			standardName = "miR-%s" % mirID

	return standardName

def cancermine(sentenceFile,modelFilenames,filterTerms,wordlistPickle,genes,cancerTypes,outData):
	print("%s : start" % now())

	models = {}
	assert isinstance(modelFilenames,list)
	for modelFilename in modelFilenames:
		with open(modelFilename,'rb') as f:
			models[modelFilename] = pickle.load(f)

	IDToTerm = {}
	with codecs.open(genes,'r','utf-8') as f:
		for line in f:
			geneid,singleterm,_ = line.strip().split('\t')
			IDToTerm[geneid] = singleterm

	with codecs.open(cancerTypes,'r','utf-8') as f:
		for line in f:
			cancerid,singleterm,_ = line.strip().split('\t')
			IDToTerm[cancerid] = singleterm

	with codecs.open(filterTerms,'r','utf-8') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	with open(wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)
	
	# Truncate the output file
	with codecs.open(outData,'w','utf-8') as outF:
		pass

	timers = Counter()

	print("%s : loading..." % now())
	with open(sentenceFile) as f:
		sentenceData = json.load(f)

	corpus = kindred.Corpus()
	for sentence in sentenceData:
		metadata = dict(sentence)
		del metadata["sentence"]
		doc = kindred.Document(sentence["sentence"],metadata=metadata)
		corpus.addDocument(doc)

	print("%s : loaded..." % now())
	startTime = time.time()
	parser = kindred.Parser()
	parser.parse(corpus)
	timers['parser'] += time.time() - startTime
	print("%s : parsed" % now())

	startTime = time.time()
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=True,detectMicroRNA=True,acronymDetectionForAmbiguity=True,mergeTerms=True)
	ner.annotate(corpus)
	timers['ner'] += time.time() - startTime
	print("%s : ner" % now())

	with codecs.open(outData,'a','utf-8') as outF:
		startTime = time.time()
		for modelname,model in models.items():
			model.predict(corpus)
		timers['predicted'] += time.time() - startTime

		print("%s : predicted" % now())

		startTime = time.time()

		for doc in corpus.documents:
			if len(doc.relations) == 0:
				continue

			eID_to_sentence = {}
			for sentence in doc.sentences:
				for eID in sentence.getEntityIDs():
					eID_to_sentence[eID] = sentence
			eID_to_entity = doc.getEntityIDsToEntities()

			for relation in doc.relations:
				sentence = eID_to_sentence[relation.entityIDs[0]]
				sentenceTextLower = sentence.text.lower()

				hasFilterTerm = any( filterTerm in sentenceTextLower for filterTerm in filterTerms )
				if not hasFilterTerm:
					continue
				#words = [ t.word for t in sentence.tokens ]
				#text = " ".join(words)

				sentenceStart = sentence.tokens[0].startPos

				relType = relation.relationType
				entityData = []
				for eID in relation.entityIDs:
					entity = eID_to_entity[eID]
					entityData.append(entity.externalID)
					entityData.append(entity.text)


					if entity.externalID.startswith('combo'):
						externalIDsplit = entity.externalID.split('|')
						standardizedTerms = [ getStandardizedTerm("",st.replace('&',';'),IDToTerm) for st in externalIDsplit[1:] ]
						standardizedTerm = "|".join(standardizedTerms)
					elif entity.externalID.startswith('mirna|'):
						standardizedTerm = standardizeMIRName(entity.externalID)
					else:
						standardizedTerm = getStandardizedTerm(entity.text,entity.externalID,IDToTerm)

					entityData.append(standardizedTerm)
					
					assert len(entity.position) == 1, "Expecting entities that are contigious and have only one start and end position within the text"
					startPos,endPos = entity.position[0]
					entityData.append(startPos - sentenceStart)
					entityData.append(endPos - sentenceStart)


				if doc.metadata["pmid"]:
					m = doc.metadata
					outData = [m["pmid"],m['title'],m["journal"],m["year"],m['section'],relType] + entityData + [sentence.text]
					outLine = "\t".join(map(str,outData))
					outF.write(outLine+"\n")

		timers['output'] += time.time() - startTime

		print("%s : output" % now())

	sys.stdout.flush()

	print("%s : done" % now())
	
	for section,sectiontime in timers.items():
		print("%s\t%f" % (section,sectiontime))


if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--sentenceFile',required=True,help='BioC XML file to use')
	parser.add_argument('--models',required=True)
	parser.add_argument('--filterTerms',required=True)
	parser.add_argument('--wordlistPickle',required=True)
	parser.add_argument('--genes',required=True)
	parser.add_argument('--cancerTypes',required=True)
	parser.add_argument('--outData',required=True)

	args = parser.parse_args()

	cancermine(args.sentenceFile,args.models.split(','),args.filterTerms,args.wordlistPickle,args.genes,args.cancerTypes,args.outData)