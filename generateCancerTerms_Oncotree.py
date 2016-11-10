"""
This script is used to build a word-list of relevant cancer specific terms from the Disease Ontology and UMLS Metathesaurus.
"""
import argparse
import sys
import codecs
from collections import defaultdict

def augmentTermList(terms):
	"""
	Adds additional spellings and plurals to a list of cancer terms
	
	Args:
		terms (list of strings): List of strings of terms
		
	Returns:
		list of augmente strings
	"""
	
	# Lower case everything (if not already done anyway)
	terms = [ t.lower() for t in terms ]
	
	# A list of short cancer terms that are acceptable (others like ALL are too general and excluded)
	acceptedShortTerms = ["gbm","aml","crc","hcc","cll"]
	
	# Filter out smaller terms except the allowed ones
	terms = [ t for t in terms if len(t) > 3 or t in acceptedShortTerms ]

	# Filter out terms with a comma
	terms = [ t for t in terms if not ',' in t ]
	
	# Try the British spelling of tumor
	tumourTerms = [ t.replace('tumor','tumour') for t in terms ]
	
	# Terms that we can add an 'S' to pluralise (if not already included)
	pluralEndings = ["tumor", "tumour", "neoplasm", "cancer", "oma"]
	
	# Check if any term ends with one of the plural endings, and then pluralise it
	plurals = []
	for t in terms:
		pluralize = False
		for e in pluralEndings:
			if t.endswith(e):
				pluralize = True
				break

		if pluralize:
			plurals.append(t + "s")

	# Sorted and unique the terms back together
	merged = sorted(list(set(terms + tumourTerms + plurals)))
	return merged

def loadMetathesaurus(filename):
	"""
	Loads the UMLS metathesaurus into a dictionary where CUID relates to a set of terms. Only English terms are included

	Args:
		filename (str): Filename of UMLS Concept file (MRCONSO.RRF)

	Returns:
		Dictionary where each key (CUID) points to a list of strings (terms)
	"""
	meta = defaultdict(list)
	with codecs.open(filename,'r','utf8') as f:
		for line in f:
			split = line.split('|')
			cuid = split[0]
			lang = split[1]
			term = split[14]
			if lang != 'ENG':
				continue
			meta[cuid].append(term)
	return meta
	
if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Generate term list from Oncotree and UMLS Metathesarus for cancer-specific terms')
	parser.add_argument('--oncotreeFile', required=True, type=str, help='Path to the Oncotree file')
	parser.add_argument('--umlsConceptFile', required=True, type=str, help='Path on the MRCONSO.RRF file in UMLS metathesaurus')
	parser.add_argument('--outFile', required=True, type=str, help='Path to output wordlist file')
	args = parser.parse_args()

	print "Loading metathesaurus..."
	metathesaurus = loadMetathesaurus(args.umlsConceptFile)

	cancerTypes = []

	print "Processing"
	with codecs.open(args.oncotreeFile,'r','utf8') as inF:
		# Skip down to the children of the cancer term and then find all their descendents (recursive children)
		header = inF.readline()

		for line in inF:
			# Split the tab-delimited line
			split = line.rstrip('\n\r').split('\t')

			# If some columns are missing. Skip it
			if len(split) < 9:
				continue

			# Pull out the appropriate IDs
			nciID = split[7]
			umlsID = split[8]

			# Check that both are filled
			if nciID == '' or umlsID == '':
				continue

			# Get the English terms for the metathesaurus
			mmterms = metathesaurus[umlsID]

			# Check that there are some terms
			if len(mmterms) > 0:
				# Lowercase everything
				mmterms = [ mmterm.lower() for mmterm in mmterms ]
				
				# Add extra spellings and plurals
				mmterms = augmentTermList(mmterms)

				# Remove any duplicates and sort it
				mmterms = sorted(list(set(mmterms)))

				# Pull out the numeric NCI ID
				numericID = int(nciID[1:])

				# Add this to our list of cancer terms
				cancerType = (numericID,nciID,mmterms)
				cancerTypes.append(cancerType)

	# Sort the list by ID
	cancerTypes = sorted(cancerTypes)
			
	# Output the list to the output file
	with codecs.open(args.outFile,'w','utf8') as outF:
		# Iterate through them and save them
		for _,nciID,mmterms in cancerTypes:
			line = u"%s\t%s" % (nciID, u"|".join(mmterms))
			outF.write(line + "\n")

	print "Successfully output to %s" % args.outFile

		

