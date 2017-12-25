import pickle
import kindred
import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Loads up a wordlist of genes and cancer types and saves to a Python pickle')
	parser.add_argument('--genes',required=True)
	parser.add_argument('--cancers',required=True)
	parser.add_argument('--conflicting',required=True)
	parser.add_argument('--wordlistPickle',required=True)

	args = parser.parse_args()

	print("Loading...")

	termLookup = kindred.EntityRecognizer.loadWordlists({'gene':args.genes,'cancer':args.cancers,'unused':args.conflicting}, idColumn=0, termsColumn=2)

	with open(args.wordlistPickle,'wb') as f:
		pickle.dump(termLookup,f)

	print("Wordlist with %d terms written to %s" % (len(termLookup),args.wordlistPickle))

