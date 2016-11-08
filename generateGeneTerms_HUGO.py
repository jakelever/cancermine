"""
This script is used to build a word-list of relevant gene terms from the HUGO gene list
"""
import argparse
import sys
import codecs

def cleanupQuotes(text):
	"""
	Removes quotes if text starts and ends with them

	Args:
		text (str): Text to cleanup

	Returns:
		Text with quotes removed from start and end (if they existed) or original string (if not)
	"""
	if text.startswith('"') and text.endswith('"'):
		return text[1:-1]
	else:
	 	return text

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Generate term list from HUGO gene resource')
	parser.add_argument('--hugoGenesFile', required=True, type=str, help='Path to HUGO gene names file')
	parser.add_argument('--outFile', required=True, type=str, help='Path to output wordlist file')
	args = parser.parse_args()

	print "Processing"
	with codecs.open(args.hugoGenesFile,'r','utf8') as hugoF, codecs.open(args.outFile,'w','utf8') as outF:
		for line in hugoF:
			split = line.split('\t')

			# Get the relevant fields for the gene
			hgnc_id = split[0]
			symbol = split[1]
			name = split[2]
			alias_symbol = split[8]
			alias_name = split[9]
			prev_symbol = split[10]
			prev_name = split[11]

			allNames = [symbol,name,alias_symbol,alias_name,prev_symbol,prev_name]
			allNames = [ x.strip() for x in allNames ]
			allNames = [ x for x in allNames if x ]
			allNames = [ cleanupQuotes(x) for x in allNames ]

			# Remove any duplicates
			merged = '|'.join(allNames)
			noDuplicates = sorted(list(set(merged.split('|'))))

			# Then output to the file
			line = u"%s\t%s" % (hgnc_id, u"|".join(noDuplicates))
			outF.write(line + "\n")

	print "Successfully output to %s" % args.outFile

		

