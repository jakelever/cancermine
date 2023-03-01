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
import html

def now():
	return time.strftime("%Y-%m-%d %H:%M:%S")

def getNormalizedTerm(text,externalID,IDToTerm):
	normalizedTerms = [ IDToTerm[eid] for eid in externalID.split(';') ]
	normalizedTerms = sorted(list(set(normalizedTerms)))

	normalizedTermsLower = [ st.lower() for st in normalizedTerms ]
	textLower = text.lower()

	if textLower in normalizedTermsLower:
		index = normalizedTermsLower.index(textLower)
		normalizedTerm = normalizedTerms[index]
	else:
		normalizedTerm = ";".join(normalizedTerms)
	return normalizedTerm

def normalizeMIRName(externalID):
	assert externalID.startswith('mirna|'), "Unexpected ID: %s" % externalID
	normalizedName = externalID[4:]

	search = re.search('mirna\|\D*(?P<id>\d+[A-Za-z]*)',externalID)
	if search:
		mirID = search.groupdict()['id']
		if not mirID is None:
			normalizedName = "miR-%s" % mirID

	return normalizedName

def getFormattedSentence(sentence,entitiesToHighlight):
	charArray = [ html.escape(c) for c in sentence.text ]

	sentenceStart = sentence.tokens[0].startPos
	for e in entitiesToHighlight:
		for startPos,endPos in e.position:
			startPos -= sentenceStart
			endPos -= sentenceStart

			try:
				charArray[startPos] = '<b>' + charArray[startPos]
				charArray[endPos-1] = charArray[endPos-1] + '</b>'
			except:
				print("ERROR in getFormattedSentence")
				print(doc.text)
				print(e.text)
				print(e.position)
				sys.exit(1)

	return "".join(charArray)

headers = ['pmid','title','journal','journal_short','year','month','day','section','subsection','role','predictprob','cancer_id','cancer_name','cancer_normalized','cancer_start','cancer_end','gene_hugo_id','gene_entrez_id','gene_name','gene_normalized','gene_start','gene_end','sentence','formatted_sentence']

def applyFinalFilter(row):
	# Filter out incorrect output with some rules

	assert len(row) == len(headers), "Number of columns in output data (%d) doesn't  match with header count (%d)" % (len(row),len(headers))

	row = { h:v for h,v in zip(headers,row) }

	# Check for the number of semicolons (suggesting a list)
	if row['sentence'].count(';') > 5:
		return False

	if row['section'] == 'back':
		return False

	return True

def cancermine(sentenceFile,modelFilenames,filterTerms,wordlistPickle,genes,cancerTypes,outData):
	print("%s : start" % now())

	models = {}
	assert isinstance(modelFilenames,list)
	for modelFilename in modelFilenames:
		with open(modelFilename,'rb') as f:
			models[modelFilename] = pickle.load(f)

	IDToTerm = {}
	Hugo2Entrez = defaultdict(lambda : 'NA')
	with codecs.open(genes,'r','utf-8') as f:
		for line in f:
			gene_hugo_id,singleterm,_,gene_entrez_id = line.strip().split('\t')
			IDToTerm[gene_hugo_id] = singleterm
			Hugo2Entrez[gene_hugo_id] = gene_entrez_id

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
	ner = kindred.EntityRecognizer(lookup=termLookup,detectFusionGenes=False,detectMicroRNA=False,acronymDetectionForAmbiguity=True,mergeTerms=True,removePathways=True)
	ner.annotate(corpus)
	timers['ner'] += time.time() - startTime
	print("%s : ner" % now())


	with codecs.open(outData,'a','utf-8') as outF:
		outF.write("\t".join(headers) + "\n")

		startTime = time.time()
		for modelname,model in models.items():
			model.predict(corpus)
		timers['predicted'] += time.time() - startTime

		print("%s : predicted" % now())

		startTime = time.time()

		for doc in corpus.documents:
			#print(doc)
			if len(doc.relations) == 0:
				continue
			if not doc.metadata["pmid"]:
				continue

			entity_to_sentence = {}
			for sentence in doc.sentences:
				for entity,tokenIndices in sentence.entityAnnotations:
					assert not entity in entity_to_sentence
					entity_to_sentence[entity] = sentence

			journal_short = doc.metadata['journal']
			if journal_short and len(journal_short) > 50:
				journal_short = journal_short[:50] + '...'

			for relation in doc.relations:
				sentence = entity_to_sentence[relation.entities[0]]
				sentenceTextLower = sentence.text.lower()

				hasFilterTerm = any( filterTerm in sentenceTextLower for filterTerm in filterTerms )
				if not hasFilterTerm:
					continue
				#words = [ t.word for t in sentence.tokens ]
				#text = " ".join(words)

				entitiesAreAmbiguous = any ( [';' in e.externalID for e in relation.entities] )
				if entitiesAreAmbiguous:
					continue

				sentenceStart = sentence.tokens[0].startPos

				skip = False

				relType = relation.relationType
				entityData = []
				for entity in relation.entities:
					entityData.append(entity.externalID)

					startPos,endPos = entity.position[0]
					if entity.entityType == 'gene':
						entityData.append(Hugo2Entrez[entity.externalID])

						afterText = doc.text[endPos:].strip()
						if afterText.startswith('-AS'):
							skip = True

					entityData.append(entity.text)


					if entity.externalID.startswith('combo'):
						externalIDsplit = entity.externalID.split('|')
						normalizedTerms = [ getNormalizedTerm("",st.replace('&',';'),IDToTerm) for st in externalIDsplit[1:] ]
						normalizedTerm = "|".join(normalizedTerms)
					elif entity.externalID.startswith('mirna|'):
						normalizedTerm = normalizeMIRName(entity.externalID)
					else:
						normalizedTerm = getNormalizedTerm(entity.text,entity.externalID,IDToTerm)

					entityData.append(normalizedTerm)
					
					assert len(entity.position) == 1, "Expecting entities that are contigious and have only one start and end position within the text"

					entityData.append(startPos - sentenceStart)
					entityData.append(endPos - sentenceStart)

				if skip:
					continue

				m = doc.metadata
				if not 'subsection' in m:
					m['subsection'] = None

				formattedSentence = getFormattedSentence(sentence,relation.entities)

				prob = relation.probability
				outData = [m['pmid'],m['title'],m["journal"],journal_short,m["year"],m["month"],m["day"],m['section'],m['subsection'],relType,prob] + entityData + [sentence.text, formattedSentence]
				if applyFinalFilter(outData):
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
