"""
This script is used to build a word-list of relevant cancer specific terms from the Disease Ontology and UMLS Metathesaurus.
"""
import argparse
import sys
import codecs
import pronto
from collections import defaultdict

def findTerm(ont,name):
	"""
	Searches an ontology for a specific term name and returns the first hit

	Args:
		ont (pronto Ontology): Ontology to search
		name (str): Search query

	Returns:
		pronto Ontology: Term that matched name or None if not found
	"""
	for term in ont:
		if term.name == name:
			return term
	return None



def getCUIDs(term):
	"""
	Gets all CUIDs for a given pronto Ontology object (from the xrefs)

	Args:
		term (pronto Ontology): Term from ontology to extract CUIDs for

	Returns:
		list of CUIDs
	"""
	cuids = []
	if 'xref' in term.other:
		for xref in term.other['xref']:
			if xref.startswith('UMLS_CUI'):
				cuid = xref[9:]
				cuids.append(cuid)
	return cuids

def loadMetathesaurus(filename):
	"""
	Loads the UMLS metathesaurus into a dictionary where CUID relates to a set of terms. Only English terms are included

	Args:
		filename (str): Filename of UMLS Concept file (MRCONSO.RRF)

	Returns:
		Dictionary where each key (CUID) points to a list of strings (terms)
	"""
	meta = defaultdict(list)
	with open(filename) as f:
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

	parser = argparse.ArgumentParser(description='Generate term list from Disease Ontology and UMLS Metathesarus for cancer-specific terms')
	parser.add_argument('--diseaseOntologyFile', required=True, type=str, help='Path to the Disease Ontology OBO file')
	parser.add_argument('--umlsConceptFile', required=True, type=str, help='Path on the MRCONSO.RRF file in UMLS metathesaurus')
	parser.add_argument('--outFile', required=True, type=str, help='Path to output wordlist file')
	args = parser.parse_args()

	print "Loading metathesaurus..."
	metathesaurus = loadMetathesaurus(args.umlsConceptFile)

	print "Loading disease ontology..."
	ont = pronto.Ontology(args.diseaseOntologyFile)
	cancerTerm = findTerm(ont,'cancer')

	print "Processing"
	with codecs.open(args.outFile,'w','utf8') as outF:
		# Skip down to the grandchildren of the cancer term and then find all their descendents (recursive children)
		for term in cancerTerm.children.children.rchildren():
			# Get the CUIDs for this term
			cuids = getCUIDs(term)

			# Get the English terms for the metathesaurus
			mmterms = [ metathesaurus[cuid] for cuid in cuids ]

			# Merge the lists together
			mmterms = sum(mmterms, [])

			# Add in the Disease Ontology term (in case it's not already in there)
			mmterms.append(term.name)

			# Lowercase everything
			mmterms = [ mmterm.lower() for mmterm in mmterms ]

			# Remove any duplicates and sort it
			mmterms = sorted(list(set(mmterms)))

			# Then output to the file
			line = "%s\t%s\t%s" % (term.id, term.name, "|".join(mmterms))
			outF.write(line + "\n")
	print "Successfully output to %s" % args.outFile

		

