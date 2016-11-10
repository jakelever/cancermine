import sys
import fileinput
import argparse
import time
from collections import Counter
import itertools
from textExtractionUtils import *
import re
from AcronymDetect import detectAcronyms

from java.util import *
from edu.stanford.nlp.pipeline import *
from edu.stanford.nlp.ling.CoreAnnotations import *
from edu.stanford.nlp.semgraph.SemanticGraphCoreAnnotations import *

pipeline = None
def getPipeline():
	global pipeline
	if pipeline is None:
		props = Properties()
		props.put("annotators", "tokenize, ssplit, pos, lemma, depparse");
		pipeline = StanfordCoreNLP(props, False)
		
	return pipeline
	
	
minipipeline = None
def getMiniPipeline():
	global minipipeline
	if minipipeline is None:
		props = Properties()
		props.put("annotators", "tokenize");
		minipipeline = StanfordCoreNLP(props, False)
	return minipipeline

def detectFusionTerms(words, lookupDict):
	termtypesAndids,terms,locs = [],[],[]
	origWords = list(words)
	words = [ w.lower() for w in words ]

	for i,word in enumerate(words):
		split = re.split("[-/]",word)
		if len(split) == 1:
			continue
			
		allGenes = True
		
		geneIDs = ['fusion']
		for s in split:
			key = (s,)
			if key in lookupDict:
				isGene = False
				for type,ids in lookupDict[key]:
					if type == 'gene':
						idsTxt = ";".join(map(str,ids))
						geneIDs.append(idsTxt)
						isGene = True
						break
				if not isGene:
					allGenes = False
					break
			else:
				allGenes = False
				break
	
		if allGenes:
			#geneTxt = ",".join(map(str,geneIDs))
			termtypesAndids.append([('gene',geneIDs)])
			terms.append(tuple(origWords[i:i+1]))
			locs.append((i,i+1))
			
	return termtypesAndids,terms,locs
	
def getTermIDsAndLocations(np, lookupDict):
	termtypesAndids,terms,locs = [],[],[]
	# Lowercase all the tokens
	#np = [ unicodeLower(w) for w in np ]
	orignp = np
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
				# If found, save the ID(s) in the dictionary
				termtypesAndids.append(lookupDict[s])
				terms.append(tuple(orignp[i:i+l]))
				locs.append((i,i+l))
				# And blank it out
				np[i:i+l] = [ "" for _ in range(l) ]

	# Then return the found term IDs
	return termtypesAndids,terms,locs

from __builtin__ import zip
def selectSentencesDebug(outFile, textInput, textSourceInfo):
	pipeline = getPipeline()

	#textInput = [ u'She noticed that a deletion of PTEN correlates with sensitivity to Erlotinib.' ]
	#textInput = [u'The V600E mutation in BRAF is known to cause resistance to Gefitinib.' ]
	print textSourceInfo
	print textInput

	pmid = str(textSourceInfo['pmid'])
	pmcid = str(textSourceInfo['pmcid'])

	#print textInput
	assert isinstance(textInput, list)
	for text in textInput:
		text = text.strip().replace('\n', ' ').replace('\r',' ').replace('\t',' ')
		text = text.replace(u'\u2028',' ').replace(u'\u2029',' ').replace(u'\u202F',' ').replace(u'\u2012',' ').replace(u'\u2010',' ')
		text = "".join(ch for ch in text if unicodedata.category(ch)[0]!="C")
		text = text.strip()

		if len(text) == 0:
			continue

		assert isinstance(text, str) or isinstance(text, unicode)

		count = 0
		document = pipeline.process(text)
		for sentence in document.get(SentencesAnnotation):
			count += 1
		print count

def parseWordlistTerm(text):
	minipipeline = getMiniPipeline()
	tokens = minipipeline.process(text)
	return tuple([ token.word() for token in tokens.get(TokensAnnotation) ])
		
from __builtin__ import zip
def selectSentences(outFile, textInput, textSourceInfo):
	pipeline = getPipeline()

	#textInput = [ u'She noticed that a deletion of PTEN correlates with sensitivity to Erlotinib.' ]
	#textInput = [ u'The V600E mutation in BRAF is known to cause resistance to Gefitinib.' ]

	pmid = str(textSourceInfo['pmid'])
	pmcid = str(textSourceInfo['pmcid'])
	pubYear = str(textSourceInfo['pubYear'])

	print "pmid:%s pmcid:%s" % (pmid,pmcid)

	driven1 = re.compile(re.escape('-driven'), re.IGNORECASE)
	driven2 = re.compile(re.escape('- driven'), re.IGNORECASE)

	#print textInput
	assert isinstance(textInput, list)
	for text in textInput:
		text = text.strip().replace('\n', ' ').replace('\r',' ').replace('\t',' ')
		text = text.replace(u'\u2028',' ').replace(u'\u2029',' ').replace(u'\u202F',' ').replace(u'\u2012',' ').replace(u'\u2010',' ')
		text = "".join(ch for ch in text if unicodedata.category(ch)[0]!="C")
		text = text.decode('utf-8','ignore').encode("utf-8")
		text = text.strip()

		#text = text.replace('-driven',' driven')
		#text = text.replace('- driven',' driven')
		#text = text.replace('-Driven',' Driven')
		#text = text.replace('- Driven',' Driven')

		text = driven1.sub(' driven',text)
		text = driven2.sub(' driven',text)

		if len(text) == 0:
			continue

		assert isinstance(text, str) or isinstance(text, unicode)

		document = pipeline.process(text)
		for sentence in document.get(SentencesAnnotation):
			sentenceStart = None
			
			words = []
			positions = []
			for i,token in enumerate(sentence.get(TokensAnnotation)):
				if sentenceStart is None:
					sentenceStart = token.beginPosition()

				word = token.word()
				startPos = token.beginPosition() - sentenceStart
				endPos = token.endPosition() - sentenceStart
				words.append(word)
				positions.append((startPos,endPos))
			
			#print "-"*30
			#print words
			
			
			snvRegex = r'^[A-Z][0-9]+[A-Z]$'
			snvMatches = [ not (re.match(snvRegex,w) is None) for w in words ]

			termtypesAndids,terms,locs = getTermIDsAndLocations(words,lookup)
			fusionTermtypesAndids,fusionTerms,fusionLocs = detectFusionTerms(words,lookup)
			
			termtypesAndids += fusionTermtypesAndids
			terms += fusionTerms
			locs += fusionLocs
			
			for i,(w,snvMatch) in enumerate(zip(words,snvMatches)):
				if snvMatch:
					termtypesAndids.append([('mutation',['snv'])])
					terms.append((w,))
					locs.append((i,i+1))

			for i,w in enumerate(words):
				if w.lower().startswith("mir-") or w.lower().startswith("hsa-mir-") or w.lower().startswith("microrna-"):
					termtypesAndids.append([('gene',['mrna'])])
					terms.append((w,))
					locs.append((i,i+1))

			#print "-"*30
			#print sentence
			#print termtypesAndids
			#print snvMatches
			
			locsToRemove = set()
			
			acronyms = detectAcronyms(words)
			for (wordsStart,wordsEnd,acronymLoc) in acronyms:
				wordIsTerm = (wordsStart,wordsEnd) in locs
				acronymIsTerm = (acronymLoc,acronymLoc+1) in locs
				
				if wordIsTerm and acronymIsTerm:
					# Remove the acronym
					locsToRemove.add((acronymLoc,acronymLoc+1))
				elif acronymIsTerm:
					# Remove any terms that contain part of the spelt out word
					newLocsToRemove = [ (i,j) for i in range(wordsStart,wordsEnd) for j in range(i,wordsEnd+1) ]
					locsToRemove.update(newLocsToRemove)
					

			zipped = zip(locs,terms,termtypesAndids)
			filtered = [ (locs,terms,termtypesAndids) for locs,terms,termtypesAndids in zipped if not locs in locsToRemove]

			#print "len(list(zipped))", len(list(zipped))
			#print "len(filtered)", len(filtered)

			cancerLocs,geneLocs = set(),set()
			for loc,term,x in filtered:
				for t,_ in x:
					if t == 'cancer':
						cancerLocs.add(loc)
					elif t == 'gene':
						geneLocs.add(loc)

			overlap = [ t for t in cancerLocs if t in geneLocs ]
			uniqCancerLocs = [ t for t in cancerLocs if not t in overlap ]
			uniqGeneLocs = [ t for t in geneLocs if not t in overlap ]
		
			hasCancerAndGeneTerms = (len(cancerLocs)>=1 and len(geneLocs)>=1) and not (len(cancerLocs) == 1 and len(geneLocs)==1 and len(overlap)==1)

			#types = set ( [ t for x in termtypesAndids for t,_ in x ] )
			#if len(types) == 3: # All three types: drugs,genes,mutationTypes
			#if "cancer" in types and "gene" in types:
			if hasCancerAndGeneTerms:
				#print words
				print "-"*30
				#print textSourceInfo
				#print sentence
				out = [pmid,pmcid,pubYear,unicode(sentence)]
				#for thesetypesAndIDs,term,(startT,endT) in zip(termtypesAndids,terms,locs):
				for (startT,endT),term,thesetypesAndIDs in filtered:
					for type,termid in thesetypesAndIDs:
						startPos = positions[startT][0]
						endPos = positions[endT-1][1]
						#termTxt = " ".join(term)
						termTxt = sentence.toString()[startPos:endPos]
						data = [ type, ",".join(map(str,termid)), startPos, endPos, termTxt ]
						txt = u"|".join(map(unicode,data))
						out.append(txt)

				outLine = "\t".join(out)
				outFile.write(outLine + "\n")
				#for (type,termid),term,loc in zip(termids,terms,locs):
				#	print type,termid,loc
				#	if type == 0:
				#		print drugs[termid]
				#	elif type == 1:
				#		print genes[termid]
				#	elif type == 2:
				#		print mutationKeywords[termid]

	#print "MOO:"
	#print lookup[(u'erlotinib',)]		

			

# Find all co-occurrences for a list of text, or just a line of text
def extractTermCoOccurrences_OneList(outFile, textInput, textSourceInfo):

	# First check if it is a list of text or just text
	if isinstance(textInput, list):
		textList = textInput
	else:
		textList = [textInput]
	
	# Go through each bit of text
	for text in textList:
		# Remove weird text issues
		text = handleEncoding(text)
		
		# Extract each sentence
		for sentence in sentenceSplit(text):
			
			# Tokenize each sentence
			tokens = tokenize(sentence)
			
			# Get the IDs of terms found in the sentence
			termIDs = getID_FromLongestTerm(tokens, idLookup1)
			
			terms = [ (i,reverseLookup[i]) for i in termIDs ]

			if sentence.lower().startswith("insulin"):
				print "-------------------------"
				print sentence
				print tokens
				#print [ unicodeLower(t) for t in tokens ]
				print terms
				print "-------------------------"
				#print
			
			# Print co-occurrences to output file immediately
			for (i,j) in itertools.product(termIDs,termIDs):
				if i<j:
					outFile.write("%d\t%d\t%d\n" % (i, j, 1))
		
# It's the main bit. Yay!
if __name__ == "__main__":

	# Arguments for the command line
	parser = argparse.ArgumentParser(description='')

	parser.add_argument('--cancerList', required=True, type=str, help='')
	parser.add_argument('--geneList', required=True, type=str, help='')
	parser.add_argument('--mutationKeywords', required=True, type=str, help='')

	parser.add_argument('--stopwordsFile',  type=str, help='A path to a stopwords file that will be removed from term-lists (e.g. the, there)')
	parser.add_argument('--removeShortwords', help='Remove short words from any term lists (<=2 length)', action='store_true')

	parser.add_argument('--abstractsFile', type=argparse.FileType('r'), help='MEDLINE file containing abstract data')
	parser.add_argument('--articleFile', type=argparse.FileType('r'), help='PMC NXML file containing a single article')
	parser.add_argument('--articleFilelist', type=argparse.FileType('r'), help='File containing filenames of multiple PMC NXML files')

	parser.add_argument('--outFile', type=str, help='File to output cooccurrences')

	args = parser.parse_args()

	print "Loading cancer list..."
	with codecs.open(args.cancerList, "r", "utf-8") as f1:
		cancers = {}
		for line in f1:
			split = line.split('\t')
			cancers[split[0]] = split[1].strip().split('|')

	print "Loading gene list..."
	with codecs.open(args.geneList, "r", "utf-8") as f2:
		genes = {}
		for line in f2:
			split = line.split('\t')
			genes[split[0]] = split[1].strip().split('|')

	print "Loading mutation list..."
	with codecs.open(args.mutationKeywords, "r", "utf-8") as f3:
		mutationKeywords = { i : line.strip().split('|') for i,line in enumerate(f3) } 


	#drugLookup = { x:i for i,x in enumerate(drugs) }
	#geneLookup = { x:i for i,x in enumerate(genes) }
	#mutationLookup = { x:i for i,x in enumerate(mutationKeywords) }

	print "Generating lookup table..."
	duplicates = set()
	lookup = defaultdict(list)
	for termType,mainDict in zip(['cancer','gene','mutation'],[cancers,genes,mutationKeywords]):
	#for type,mainList in enumerate([drugs,genes,mutationKeywords]):
		for id,lst in mainDict.iteritems():
			lst = list(set(lst))
			keys = set( [ parseWordlistTerm(x.lower()) for x in lst ] )
			#for x in lst:
			for key in keys:
				#key = tuple(x.lower().split(' '))
				#key = parseWordlistTerm(x.lower())
				
				matching = None
				if key in lookup:
					prevItems = lookup[key]
					
					matching = None
					for i,(prevType,prevIds) in enumerate(prevItems):
						if prevType == termType:
							matching = i
							break
					
				if matching is None:
					item = (termType,[id])
					lookup[key].append(item)
				else:
					prevItems[matching][1].append(id)
					lookup[key] = prevItems
			 
	stopwords = []
	if args.stopwordsFile:
		with codecs.open(args.stopwordsFile, "r", "utf-8") as f4:
			stopwords = [ line.strip() for line in f4 ]

		print "Removing stopwords..."
		for stopword in stopwords:
			key = tuple(stopword.lower().split(' '))
			if key in lookup:
				del lookup[key]

	if args.removeShortwords:
		print "Removing short words..."
		lookup = { key:val for key,val in lookup.iteritems() if not (len(key)==1 and len(key[0]) <= 2) }

	outFile = codecs.open(args.outFile, "w", "utf-8")

	print "Starting processing..."
	startTime = time.time()
	# And now we try to process either an abstract file, single article file or multiple
	# article files
	try:
		if args.abstractsFile:
			processAbstractFile(args.abstractsFile, outFile, selectSentences)
		elif args.articleFile:
			# Just pull the filename and pass that, instead of the object
			filename = args.articleFile.name
			args.articleFile.close()
			processArticleFiles(filename, outFile, selectSentences)
		elif args.articleFilelist:
			# Extract the file list from another file
			fileList = [ f.strip() for f in args.articleFilelist]
		
			processArticleFiles(fileList, outFile, selectSentences)
	except:
		print "Unexpected error:", sys.exc_info()[0]
		print "COMMAND: " + " ".join(sys.argv)
		raise

	endTime = time.time()
	duration = endTime - startTime
	print "Processing Time: ", duration
	
	# Print completion if we were working on an outFile
	if args.outFile:
		print "Finished output to:", args.outFile
	

