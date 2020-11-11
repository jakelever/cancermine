import argparse
from collections import Counter
import string

import kindred
import pickle

def remove_punctuation(text):
	exclude = set(string.punctuation)
	return ''.join(ch for ch in text if ch not in exclude)

def runPredictionsOnSentence(sentence,models,termLookup):
	assert isinstance(sentence,str)

	corpus = kindred.Corpus()
	doc = kindred.Document(sentence)
	corpus.addDocument(doc)

	parser = kindred.Parser()
	parser.parse(corpus)

	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=False,detectMicroRNA=False,acronymDetectionForAmbiguity=True,mergeTerms=True,removePathways=True)
	ner.annotate(corpus)

	for modelname,model in models.items():
		model.predict(corpus)

	return doc

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Tool to explain differences between versions')
	parser.add_argument('--newSentences',required=True,type=str,help='New version of cancermine_sentences.tsv')
	parser.add_argument('--newUnfiltered',required=True,type=str,help='New version of cancermine_unfiltered.tsv')
	parser.add_argument('--oldSentences',required=True,type=str,help='Old version of cancermine_sentences.tsv')
	parser.add_argument('--oldUnfiltered',required=True,type=str,help='Old version of cancermine_unfiltered.tsv')

	parser.add_argument('--models',required=True,help='Comma-delimited list of model files to use')
	parser.add_argument('--wordlistPickle',required=True,help='Pickle of wordlist')
	args = parser.parse_args()

	with open(args.wordlistPickle,'rb') as f:
		termLookup = pickle.load(f)

	models = {}
	for modelFilename in args.models.split(','):
		with open(modelFilename,'rb') as f:
			models[modelFilename] = pickle.load(f)
	#thresholds = {'Driver':0.80, 'Oncogene': 0.76, 'Tumor_Suppressor': 0.92}

	newSentences = set()
	#oldPMIDs = set()
	newUnfiltered = set()

	reasonCounter = Counter()

	with open(args.newSentences) as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			d = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }

			cleaned_sentence = remove_punctuation(d['sentence']).lower().strip()
			
			record = ( d['pmid'], d['sentence'], d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			altrecord = ( d['pmid'], cleaned_sentence, d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			newSentences.add(record)
			newSentences.add(altrecord)

	with open(args.newUnfiltered) as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			d = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }

			sentence = remove_punctuation(d['sentence']).lower().strip()

			record = ( d['pmid'], d['sentence'], d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			altrecord = ( d['pmid'], cleaned_sentence, d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			newUnfiltered.add(record)
			newUnfiltered.add(altrecord)

	with open(args.oldSentences) as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			d = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }

			sentence = remove_punctuation(d['sentence']).lower().strip()

			record = ( d['pmid'], d['sentence'], d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			altrecord = ( d['pmid'], cleaned_sentence, d['role'], d['cancer_id'], d['gene_hugo_id'], d['gene_entrez_id'] )
			#oldSentences.add(oldSentence)

			if not record in newSentences:
				reason = 'Unknown'
				if record in newUnfiltered:
					reason = 'Missed threshold'
				elif altrecord in newSentences or altrecord in newUnfiltered:
					reason = 'Parsing change'
				else:
					doc = runPredictionsOnSentence(d['sentence'],models,termLookup)
					gene_hugo_ids = [ e.externalID for e in doc.entities if e.entityType == 'gene' ]
					cancer_ids = [ e.externalID for e in doc.entities if e.entityType == 'cancer' ]

					relations = [ (r.relationType,r.entities[0].externalID,r.entities[1].externalID) for r in doc.relations ]

					if (d['role'],d['cancer_id'],d['gene_hugo_id']) in relations:
						reason = 'Unknown'
					else:
						if d['gene_hugo_id'] in gene_hugo_ids and d['cancer_id'] in cancer_ids:
							reason = 'ML change' # no longer predicting so model change
						elif d['gene_hugo_id'] in gene_hugo_ids:
							reason = 'Cancer name change'
						elif d['cancer_id'] in cancer_ids:
							reason = 'Gene name change'
						else:
							reason = 'Cancer/Gene name change'

					if True and reason == 'Unknown':
						relationsWithProbs = [ (r.relationType,r.entities[0].externalID,r.entities[1].externalID,r.probability) for r in doc.relations ]
						print(doc)
						print(relationsWithProbs)
						print(record)
						print(reason)
						print("-"*30)
						#break
					

				#if reason == 'Unknown':


				reasonCounter[reason] += 1


	print(reasonCounter)

