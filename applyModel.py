import sys
import kindred
import pickle
import argparse
import os
from os import listdir
from os.path import isfile, join
import codecs
import re
from collections import defaultdict,Counter

# Given a tokenized bit of text, find all the words that
# are in a lookup dictionary. Find longest terms first.
def getID_FromLongestTerm(np, lookupDict):
	terms = []
	# Lowercase all the tokens
	np = [ w.lower() for w in np ]

	# The length of each search string will decrease from the full length
	# of the text down to 1
	for l in reversed(range(1, len(np)+1)):
		# We move the search window through the text
		for i in range(len(np)-l+1):
			# Extract that window of text
			s = tuple(np[i:i+l])
			# Search for it in the dictionary
			if s in lookupDict:
				# If found, save the ID(s) in the dictionar
				#terms = terms + lookupDict[s]
				found = defaultdict(list)
				for entityType,entityID in lookupDict[s]:
					found[entityType].append(entityID)

				for entityType,entityIDs in found.items():
					terms.append((entityType,entityIDs,i,i+l))
				# And blank it out
				np[i:i+l] = [ "" for _ in range(l) ]

	# Then return the found term IDs
	return terms

def loadWordlist(filename):
	lookup = defaultdict(set)
	with codecs.open(filename,'r','utf-8') as f:
		for line in f:
			termid,terms = line.strip().split('\t')
			for term in terms.split('|'):
				tupleterm = tuple(term.split())
				lookup[tupleterm].add(termid)

	lookup = { k:sorted(list(kset)) for k,kset in lookup.items() }

	return lookup

def loadWordlists(entityTypesWithFilenames):
	lookup = defaultdict(set)
	for entityType,filename in entityTypesWithFilenames.items():
		with codecs.open(filename,'r','utf-8') as f:
			for line in f:
				termid,terms = line.lower().strip().split('\t')
				for term in terms.split('|'):
					tupleterm = tuple(term.split())
					lookup[tupleterm].add((entityType,termid))

	lookup = { k:sorted(list(kset)) for k,kset in lookup.items() }

	return lookup

import time
def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Finds relations in Pubmed file')
	parser.add_argument('--biocFile',required=True,help='BioC XML file to use')
	parser.add_argument('--inModel',required=True)
	parser.add_argument('--filterTerms',required=True)
	parser.add_argument('--outData',required=True)
	parser.add_argument('--genes',required=True)
	parser.add_argument('--cancerTypes',required=True)

	args = parser.parse_args()

	print("%s : start" % now())

	with open(args.inModel,'rb') as f:
		classifier = pickle.load(f)

	with codecs.open(args.filterTerms,'r','utf-8') as f:
		filterTerms = [ line.strip().lower() for line in f ]

	termLookup = loadWordlists({'gene':args.genes,'cancer':args.cancerTypes})
		
	# Truncate the output file
	with codecs.open(args.outData,'w','utf-8') as outF:
		pass

	timers = Counter()

	print("%s : processing..." % now())
	parser = kindred.Parser()
	for corpusno,corpus in enumerate(kindred.iterLoadDataFromBioc(args.biocFile)):
		startTime = time.time()
		parser.parse(corpus)
		timers['parser'] += time.time() - startTime
		print("%s : parsed" % now())

		startTime = time.time()
		for doc in corpus.documents:
			#print(doc.sourceIDs)
			sys.stdout.flush()
			for sentence in doc.sentences:
				words = [ t.word for t in sentence.tokens ]
				termIDs = getID_FromLongestTerm(list(words),termLookup)
				for entityType,externalIDs,startToken,endToken in termIDs:
					text = " ".join(words[startToken:endToken])
					startPos = sentence.tokens[startToken].startPos
					endPos = sentence.tokens[endToken-1].endPos
					loc = list(range(startToken,endToken))
					e = kindred.Entity(entityType,text,[(startPos,endPos)],externalID=externalIDs)
					doc.addEntity(e)
					sentence.addEntityWithLocation(e,loc)
		timers['entitiesAdded'] += time.time() - startTime

		print("%s : entities added" % now())

		startTime = time.time()
		classifier.predict(corpus)
		timers['predicted'] += time.time() - startTime

		print("%s : predicted" % now())

		startTime = time.time()
		with codecs.open(args.outData,'a','utf-8') as outF:
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

					relType = relation.relationType
					entityData = []
					for eID in relation.entityIDs:
						entity = eID_to_entity[eID]
						entityData.append(entity.text)
						entityData.append(";".join(entity.externalID))

					if doc.metadata["pmid"]:
						m = doc.metadata
						outData = [m["pmid"],m['title'],m["journal"],m["year"],m['section'],relType] + entityData + [sentence.text]
						outLine = "\t".join(outData)
						outF.write(outLine+"\n")

		timers['output'] += time.time() - startTime

		print("%s : output" % now())

		#if corpusno > 5:
		#	break
		#break
		sys.stdout.flush()

	print("%s : done" % now())
	
	for section,time in timers.items():
		print("%s\t%f" % (section,time))

