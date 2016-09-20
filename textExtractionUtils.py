import re
import HTMLParser
import unicodedata
import xml.etree.cElementTree as etree
import os.path
import codecs
import pickle
from collections import defaultdict

# Wrappers for GeniaTagger and LingPipe scripts. Loaded with loadParsingTools()
tagger = None
sentencesplitter = None

# Instantiate the wrappers for various parse tools
def loadParsingTools():
	global tagger, sentencesplitter
	
	# Get the locations of geniatagger and lingpipe relative to this script
	scriptPath = os.path.dirname(os.path.realpath(__file__))
	geniaPath = scriptPath + "/../Dependencies/geniatagger-3.0.1/geniatagger"
	lingpipePath = scriptPath + "/../Dependencies/LingPipeSentenceSplitter/run.sh"
	
	# Check they are there
	if not os.path.isfile(geniaPath):
		raise RuntimeError("Cannot access GeniaTagger. Tried: " + geniaPath)
	elif not os.path.isfile(lingpipePath):
		raise RuntimeError("Cannot access LingPipe. Tried: " + lingpipePath)
	
	tagger = GeniaTagger(geniaPath)
	sentencesplitter = LingPipe(lingpipePath)

# Simple sentence splitter wrapper that calls Lingpipe
def sentenceSplit(text):
	return sentencesplitter.parse(text)

# Remove control characters from text and some other weird stuff
def handleEncoding(text):
	# Remove some "control-like" characters (left/right separator)
	text = text.replace(u'\u2028',' ').replace(u'\u2029',' ')
	text = "".join(ch for ch in text if unicodedata.category(ch)[0]!="C")
	text = text.encode('utf8')
	return text.strip()
	
# Unescape HTML special characters e.g. &gt; is changed to >
htmlParser = HTMLParser.HTMLParser()
def htmlUnescape(text):
	return htmlParser.unescape(text)
	
# Remove empty brackets (that could happen if the contents have been removed already
# e.g. for citation ( [3] [4] ) -> ( ) -> nothing
def removeBracketsWithoutWords(text):
	fixed = re.sub(r'\([\W\s]*\)', ' ', text)
	fixed = re.sub(r'\[[\W\s]*\]', ' ', fixed)
	fixed = re.sub(r'\{[\W\s]*\}', ' ', fixed)
	return fixed
	
# Some older articles have titles like "[A study of ...]."
# This removes the brackets while retaining the full stop
def removeWeirdBracketsFromOldTitles(titleText):
	titleText = titleText.strip()
	if titleText[0] == '[' and titleText[-2:] == '].':
		titleText = titleText[1:-2] + '.'
	return titleText
	
def getMaxWordLength(text):
	words = text.split()
	wordLengths = [ len(w) for w in words ]
	return max(wordLengths)
	
# Tokenize text into a tuple of the word/tokens
def tokenize(text):
	if text == '':
		return tuple()
		
	# Do a basic test for crazy length words
	maxWordLength = getMaxWordLength(text)
	if maxWordLength >= 900:
		return tuple()

	try:
		parse = tagger.parse(text.strip())
	except IOError, e:
		# Reraise the IO error with more context
		raise IOError('Error passing information to GeniaTagger subprocess. Likely that GeniaTagger has crashed.')
	tokens = []
	for (w,_,_,iob,_) in parse:
		tokens.append(w)
	return tuple(tokens)
	
# Tokenize text into a tuple of the word/tokens
def tokenizeWithPOS(text):
	# Do a basic test for crazy length words
	maxWordLength = getMaxWordLength(text)
	if maxWordLength >= 900:
		return tuple()

	try:
		parse = tagger.parse(text.strip())
	except IOError, e:
		# Reraise the IO error with more context
		raise IOError('Error passing information to GeniaTagger subprocess. Likely that GeniaTagger has crashed.')
	tokens = []
	for (w,_,pos,_,_) in parse:
		tokens.append((w,pos))
	return tuple(tokens)
	
# Extract the noun-phrases from a sentence and return them
# as a list of lists of tokens
def extractNounphrasesFromSentence(sentence):
	# Do a basic test for crazy length words
	maxWordLength = getMaxWordLength(text)
	if maxWordLength >= 900:
		return []
		
	try:
		parse = tagger.parse(sentence.strip())
	except IOError, e:
		# Reraise the IO error with more context
		raise IOError('Error passing information to GeniaTagger subprocess. Likely that GeniaTagger has crashed.')
	nounphrases = []
	current = []
	for (w,_,_,iob,_) in parse:
		if iob == 'B-NP' or iob == 'I-NP':
			current.append(w.lower())
		elif len(current) > 0:
			nounphrases.append(current)
			current = []
	if len(current) > 0:
		nounphrases.append(current)

	return nounphrases
	


def unicodeLower(str):
	return str.decode('utf8').lower().encode('utf8')
	
# Given a tokenized bit of text, find all the words that
# are in a lookup dictionary. Find longest terms first.
def getID_FromLongestTerm(np, lookupDict):
	terms = []
	# Lowercase all the tokens
	np = [ unicodeLower(w) for w in np ]
	
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
				terms = terms + lookupDict[s]
				# And blank it out
				np[i:i+l] = [ "" for _ in range(l) ]
				
	# Then return the found term IDs
	return terms
	
# Code for extracting text from Medline/PMC XML files

# XML elements to ignore the contents of
ignoreList = ['table', 'table-wrap', 'xref', 'disp-formula', 'inline-formula', 'ref-list', 'bio', 'ack', 'graphic', 'media', 'tex-math', 'mml:math', 'object-id', 'ext-link']

# XML elements to separate text between
separationList = ['title', 'p', 'sec', 'break', 'def-item', 'list-item', 'caption']
def extractTextFromElem(elem):
	textList = []
	
	# Extract any raw text directly in XML element or just after
	head = ""
	if elem.text:
		head = elem.text
	tail = ""
	if elem.tail:
		tail = elem.tail
	
	# Then get the text from all child XML nodes recursively
	childText = []
	for child in elem:
		childText = childText + extractTextFromElem(child)
		
	# Check if the tag should be ignore (so don't use main contents)
	if elem.tag in ignoreList:
		return [tail.strip()]
	# Add a zero delimiter if it should be separated
	elif elem.tag in separationList:
		return [0] + [head] + childText + [tail]
	# Or just use the whole text
	else:
		return [head] + childText + [tail]
	

# Merge a list of extracted text blocks and deal with the zero delimiter
def extractTextFromElemList_merge(list):
	textList = []
	current = ""
	# Basically merge a list of text, except separate into a new list
	# whenever a zero appears
	for t in list:
		if t == 0: # Zero delimiter so split
			if len(current) > 0:
				textList.append(current)
				current = ""
		else: # Just keep adding
			current = current + " " + t
			current = current.strip()
	if len(current) > 0:
		textList.append(current)
	return textList
	
# Main function that extracts text from XML element or list of XML elements
def extractTextFromElemList(elemList):
	textList = []
	# Extracts text and adds delimiters (so text is accidentally merged later)
	if isinstance(elemList, list):
		for e in elemList:
			textList = textList + extractTextFromElem(e) + [0]
	else:
		textList = extractTextFromElem(elemList) + [0]

	# Merge text blocks with awareness of zero delimiters
	mergedList = extractTextFromElemList_merge(textList)
	
	# Remove any newlines (as they can be trusted to be syntactically important)
	mergedList = [ text.replace('\n', ' ') for text in mergedList ]
	
	return mergedList
	
# Process a MEDLINE abstract file
# Pass in the file object, the mode to parse it with and whether to merge the output
def processAbstractFile(abstractFile, outFile, processFunction):
	count = 0
	
	# These XML files are huge, so skip through each MedlineCitation element using etree
	for event, elem in etree.iterparse(abstractFile, events=('start', 'end', 'start-ns', 'end-ns')):
		if (event=='end' and elem.tag=='MedlineCitation'):
			count = count + 1
			
			# Find the elements for the PubMed ID, and publication date information
			pmid = elem.findall('./PMID')
			yearFields = elem.findall('./Article/Journal/JournalIssue/PubDate/Year')
			medlineDateFields = elem.findall('./Article/Journal/JournalIssue/PubDate/MedlineDate')

			# Try to extract the pmidID
			pmidText = ''
			if len(pmid) > 0:
				pmidText = " ".join( [a.text.strip() for a in pmid if a.text ] )
			pmcidText = ''
				
			# Try to extract the publication date
			pubYear = 0
			if len(yearFields) > 0:
				pubYear = yearFields[0].text
			if len(medlineDateFields) > 0:
				pubYear = medlineDateFields[0].text[0:4]
				
			# Extract the title of paper
			title = elem.findall('./Article/ArticleTitle')
			titleText = extractTextFromElemList(title)
			titleText = [ removeWeirdBracketsFromOldTitles(t) for t in titleText ]
			
			# Extract the abstract from the paper
			abstract = elem.findall('./Article/Abstract/AbstractText')
			abstractText = extractTextFromElemList(abstract)
			
			# Combine all the text we want to process
			allText = titleText + abstractText
			allText = [ t for t in allText if len(t) > 0 ]
			allText = [ htmlUnescape(t) for t in allText ]
			allText = [ removeBracketsWithoutWords(t) for t in allText ]
			
			# Information about the source of this text
			textSourceInfo = {'pmid':pmidText, 'pmcid':pmcidText, 'pubYear':pubYear}
			
			# Get the co-occurrences using a single list
			processFunction(outFile, allText, textSourceInfo)
			
			# Important: clear the current element from memory to keep memory usage low
			elem.clear()
			
	
def getMetaInfoForPMCArticle(articleElem):
	# Attempt to extract the PubMed ID, PubMed Central IDs and DOIs
	pmidText = ''
	pmcidText = ''
	doiText = ''
	article_id = articleElem.findall('./front/article-meta/article-id') + articleElem.findall('./front-stub/article-id')
	for a in article_id:
		if a.text and 'pub-id-type' in a.attrib and a.attrib['pub-id-type'] == 'pmid':
			pmidText = a.text.strip().replace('\n',' ')
		if a.text and 'pub-id-type' in a.attrib and a.attrib['pub-id-type'] == 'pmc':
			pmcidText = a.text.strip().replace('\n',' ')
		if a.text and 'pub-id-type' in a.attrib and a.attrib['pub-id-type'] == 'doi':
			doiText = a.text.strip().replace('\n',' ')
			
	# Attempt to get the publication date
	pubdates = articleElem.findall('./front/article-meta/pub-date') + articleElem.findall('./front-stub/pub-date')
	pubYear = ""
	if len(pubdates) >= 1:
		pubYear = pubdates[0].find("year").text.strip().replace('\n',' ')
			
	return pmidText,pmcidText,doiText,pubYear
	
# Process a block of PubMed Central files
# Pass in the list of filenames, the mode to parse it with and whether to merge the output
def processArticleFiles(filelist, outFile, processFunction):
	if not isinstance(filelist, list):
		filelist = [filelist]

	# Go through the list of filenames and open each one
	for filename in filelist:
		with open(filename, 'r') as openfile:

			# Skip to the article element in the file
			for event, elem in etree.iterparse(openfile, events=('start', 'end', 'start-ns', 'end-ns')):
				if (event=='end' and elem.tag=='article'):
				
					pmidText,pmcidText,doiText,pubYear = getMetaInfoForPMCArticle(elem)

					# We're going to process the main article along with any subarticles
					# And if any of the subarticles have distinguishing IDs (e.g. PMID), then
					# that'll be used, otherwise the parent article IDs will be used
					subarticles = [elem] + elem.findall('./sub-article')
					
					for articleElem in subarticles:
						if articleElem == elem:
							# This is the main parent article. Just use its IDs
							subPmidText,subPmcidText,subDoiText,subPubYear = pmidText,pmcidText,doiText,pubYear
						else:
							# Check if this subarticle has any distinguishing IDs and use them instead
							subPmidText,subPmcidText,subDoiText,subPubYear = getMetaInfoForPMCArticle(articleElem)
							if subPmidText=='' and subPmcidText == '' and subDoiText == '':
								subPmidText,subPmcidText,subDoiText = pmidText,pmcidText,doiText
							if subPubYear == '':
								subPubYear = pubYear
								
							
						# Information about the source of this text
						textSourceInfo = {'pmid':subPmidText, 'pmcid':subPmcidText, 'doi':subDoiText, 'pubYear':subPubYear}
							
						
						# Extract the title of paper
						title = articleElem.findall('./front/article-meta/title-group/article-title') + articleElem.findall('./front-stub/title-group/article-title')
						assert len(title) <= 1
						titleText = extractTextFromElemList(title)
						titleText = [ removeWeirdBracketsFromOldTitles(t) for t in titleText ]
						
						# Get the subtitle (if it's there)
						subtitle = articleElem.findall('./front/article-meta/title-group/subtitle') + articleElem.findall('./front-stub/title-group/subtitle')
						subtitleText = extractTextFromElemList(subtitle)
						subtitleText = [ removeWeirdBracketsFromOldTitles(t) for t in subtitleText ]
						
						# Extract the abstract from the paper
						abstract = articleElem.findall('./front/article-meta/abstract') + articleElem.findall('./front-stub/abstract')
						abstractText = extractTextFromElemList(abstract)
						
						# Extract the full text from the paper as well as supplementaries and floating blocks of text
						articleText = extractTextFromElemList(articleElem.findall('./body'))
						backText = extractTextFromElemList(articleElem.findall('./back'))
						floatingText = extractTextFromElemList(articleElem.findall('./floats-group'))
						
						# Combine all the text we want to process
						allText = titleText + subtitleText + abstractText + articleText + backText + floatingText
						allText = [ t for t in allText if len(t) > 0 ]
						allText = [ htmlUnescape(t) for t in allText ]
						allText = [ removeBracketsWithoutWords(t) for t in allText ]
						
						# Get the co-occurrences using a single list
						processFunction(outFile, allText, textSourceInfo)
				
					# Less important here (compared to abstracts) as each article file is not too big
					elem.clear()

# Load a word-list file into a dictionary with IDs
# Allow removal of stopwords and short words
# Also can load directly from a pickled file or save to a pickled file		
def loadWordlistFile(wordlistPath, stopwordsFile, removeShortwords, binaryTermsFile, binaryTermsFile_out):
	
	# Load the word-list directly from a pickled file
	if binaryTermsFile:
		print "Loading idLookup from binary file:", binaryTermsFile.name
		wordlist = pickle.load(binaryTermsFile)
		print "Load complete."
	# Load a word-list from a file if needed
	elif os.path.isfile(wordlistPath):
		# Open the file with a unicode object
		f = codecs.open(wordlistPath, encoding='utf-8')
		
		print "Loading terms list with synonyms..."
		# For each term, split it using the delimiter and insert it into the dictionary
		wordlist = defaultdict(list)
		for (i,terms) in enumerate(f):
			for term in terms.split('|'):
				term = handleEncoding(term.strip().lower())
				key = tokenize(term)
				
				# Insert the ID (i) into a list (for that key)
				wordlist[key].append(i)
		print "Completed loading of " + str(len(wordlist)) + " terms."
			
		# Remove words from a set of stop-words, e.g. from the NLTK
		if stopwordsFile:
			print "Removing stopwords..."
			for stopword in stopwordsFile:
				stopword = tokenize(stopword)
				if stopword in wordlist:
					print "  Removing " + str(stopword) + " from termlist 1"
					del wordlist[stopword]
			print "Completed removal of stopwords"

		# Filter out words less than two characters long
		if removeShortwords:
			print "Removing short words..."
			before = len(wordlist)
			wordlist = {word:id for word,id in wordlist.iteritems() if sum(map(len,word)) > 2 }
			after = len(wordlist)
			print "Completed removal of " + str(before-after) + " short words"
	else:
		raise RuntimeError('No appropriate word-list source')

	# Or save it to a pickled file
	if binaryTermsFile_out:
		print "Saving idLookup to binary file:", binaryTermsFile_out.name
		pickle.dump(wordlist, binaryTermsFile_out)
		print "Save complete."
		
	return wordlist
