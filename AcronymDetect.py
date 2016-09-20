import itertools
import sys

def acronymMatch(words,pos,currentAcronym,atStart,subpos=None):
	if len(currentAcronym) == 0:
		if not (subpos is None): # Can't finish acronym mid-word
			return []
		else:
			return [pos+1]

	curWord = words[pos].lower()
	wordSplit = curWord.split('-')
	curLetter = currentAcronym[-1]
	
	#print curWord,curLetter
	
	moves = []
	
	if subpos is None:
		if atStart and curLetter == 's' and curWord[-1] == 's':
			# Possible plural
			moves.append( (words,pos,currentAcronym[:-1],False) )
			
		if curLetter == curWord[0]:
			moves.append( (words,pos-1,currentAcronym[:-1],False) )

	if len(wordSplit) > 1:
		if subpos is None:
			subpos = len(wordSplit)-1
		if len(wordSplit[subpos]) > 0 and curLetter == wordSplit[subpos][0]:
			if subpos == 0:
				moves.append( (words,pos-1,currentAcronym[:-1],False) )
			else:
				moves.append( (words,pos,currentAcronym[:-1],False,subpos-1) )
			
	possibleStarts = []
	for move in moves:
		possibleStarts += acronymMatch(*move)
		
	return possibleStarts

def detectAcronyms(words):
	#print words
	#sys.exit(0)
	LRBs = [i for i, x in enumerate(words) if x == u'-LRB-']
	RRBs = [i for i, x in enumerate(words) if x == u'-RRB-']
	acronyms = []
	for i,j in itertools.product(LRBs,RRBs):
		if j-i == 2:
			acronymLoc = i+1
			possibleAcronym = words[acronymLoc]
			possibleStarts = acronymMatch(words,i-1,possibleAcronym.lower(),True)
			#print possibleStarts
			if len(possibleStarts) > 0:
				start = min(possibleStarts)
				end = i
				acronyms.append((start,end,acronymLoc))
	return acronyms
	
def detectAcronyms_old(words):
	LRBs = [i for i, x in enumerate(words) if x == u'-LRB-']
	RRBs = [i for i, x in enumerate(words) if x == u'-RRB-']
	
	acronyms = []
	for i,j in itertools.product(LRBs,RRBs):
		if j-i == 2:
			acronymLoc = i+1
			possibleAcronym = words[acronymLoc]
			isAcronym = True
			start = i-len(possibleAcronym)
			end = i
			for k,letter in enumerate(possibleAcronym.lower()):
				print k,letter,words[start+k]
				if letter != words[start+k].lower()[0]:
					isAcronym = False
					break

			if isAcronym:
				acronyms.append((start,end,acronymLoc))
	return acronyms

if __name__ == "__main__":
	sentence = u'Here , we report our results in the assessment of the expression of miRNAs in HPR non-small-cell-lung-cancers -LRB- NSCLCs -RRB- .'
	#words = [u'Here', u',', u'we', u'report', u'our', u'results', u'in', u'the', u'assessment', u'of', u'the', u'expression', u'of', u'miRNAs', u'in', u'HPR', u'non', u'small', u'cell', u'lung', u'cancer', u'-LRB-', u'NSCLC', u'-RRB-', u'.']
	words = sentence.split(' ')

	acronyms = detectAcronyms(words)

	for (wordsStart,wordsEnd,acronymLoc) in acronyms:
		print words[wordsStart:wordsEnd]
		print words[acronymLoc:acronymLoc+1]
		print "-"*30
