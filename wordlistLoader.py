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
import spacy

nlp = spacy.load('en')

def loadWordlists(entityTypesWithFilenames):
	global nlp

	lookup = defaultdict(set)
	for entityType,filename in entityTypesWithFilenames.items():
		with codecs.open(filename,'r','utf-8') as f:
			tempLookup = defaultdict(set)
			for line in f:
				termid,terms = line.strip().split('\t')
				for term in terms.split('|'):
					#tupleterm = tuple(term.split())
					tupleterm = tuple([ token.text.lower() for token in nlp(term) ])
					#lookup[tupleterm].add((entityType,termid))
					tempLookup[tupleterm].add(termid)

		for tupleterm,idlist in tempLookup.items():
			lookup[tupleterm].add( (entityType,";".join(sorted(list(idlist)))) )

	return lookup

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Loads up a wordlist of genes and cancer types and saves to a Python pickle')
	parser.add_argument('--genes',required=True)
	parser.add_argument('--cancerTypes',required=True)
	parser.add_argument('--wordlistPickle',required=True)

	args = parser.parse_args()

	print("Loading...")

	termLookup = loadWordlists({'gene':args.genes,'cancer':args.cancerTypes})

	with open(args.wordlistPickle,'wb') as f:
		pickle.dump(termLookup,f)

	print("Wordlist with %d terms written to %s" % (len(termLookup),args.wordlistPickle))

