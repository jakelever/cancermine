import argparse
import csv
import hashlib
from collections import defaultdict

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter Cancermine for more conservative predictions')
	parser.add_argument('--inUnfiltered',required=True,type=str,help='Cancermine TSV to filter')
	parser.add_argument('--outCollated',required=True,type=str,help='Output collated and filtered data')
	parser.add_argument('--outSentences',required=True,type=str,help='Output filtered sentences that match collated data')
	args = parser.parse_args()

	thresholds = {'Driver':0.80, 'Oncogene': 0.76, 'Tumor_Suppressor': 0.92}

	collated = defaultdict(set)
	collatedMatchingID = {}

	sentenceOrdering = set()
	sentenceData = {}
	geneLocs = defaultdict(set)
	cancerLocs = defaultdict(set)

	readCount,writeCount = 0,0
	with open(args.inUnfiltered,'r') as inF:
		inTSV = csv.reader(inF,delimiter='\t')
		
		headers = next(inTSV, None)

		for row in inTSV:
			tmp = ' pmid,title,journal,year,month,day,section,subsection,role,predictprob,cancer_id,cancer_name,cancer_normalized,cancer_start,cancer_end,gene_hugo_id,gene_entrez_id,gene_name,gene_normalized,gene_start,gene_end,sentence '

			r = { h:v for h,v in zip(headers,row) }
			predictprob = float(r['predictprob'])
			role = r['role']
			readCount += 1

			ambigiousIDs = False
			for field in [r['cancer_id'],r['cancer_normalized'],r['gene_hugo_id'],r['gene_normalized']]:
				if ';' in field or '|' in field:
					ambigiousIDs = True
					break

			if ambigiousIDs:
				continue

			if predictprob > thresholds[role]:
				#outTSV.writerow(row)
				pmid = r['pmid']
				if pmid == 'None':
					continue # Uncitable so we skip it

				r['journal_short'] = r['journal']
				if len(r['journal_short']) > 50:
					r['journal_short'] = r['journal_short'][:50] + '...'

				collatedKeyFields = 'role,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized'
				collatedKey = tuple( [ r[k] for k in collatedKeyFields.split(',') ] )
				collated[collatedKey].add(pmid)

				# Make a field using the key data that can be used to match between tables
				matchingID = hashlib.md5("|".join(list(collatedKey)).encode('utf-8')).hexdigest()
				collatedMatchingID[collatedKey] = matchingID

				sentenceKeyFields = 'pmid,title,journal,journal_short,year,month,day,section,subsection,role,predictprob,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,sentence'
				sentenceKey = tuple( [ r[k] for k in sentenceKeyFields.split(',') ] )
				sentenceData[sentenceKey] = matchingID
				sentenceOrdering.add( (r['year'],sentenceKey) )
				geneLocs[sentenceKey].add((int(r['gene_start']),int(r['gene_end'])))
				cancerLocs[sentenceKey].add((int(r['cancer_start']),int(r['cancer_end'])))


	with open(args.outCollated,'w') as outF:
		headers = 'matching_id,role,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,citation_count'
		headerCount = len(headers.split(','))
		outF.write(headers.replace(',','\t') + '\n')

		collatedCounts = [ (len(pmids),key) for key,pmids in collated.items() ]
		collatedCounts = sorted(collatedCounts,reverse=True)
		for citation_count,collatedKey in collatedCounts:

			role,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized = collatedKey

			matchingID = collatedMatchingID[collatedKey]

			outData = [matchingID,role,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,str(citation_count)]
			assert len(outData) == headerCount

			outLine = "\t".join(outData)
			outF.write(outLine + "\n")

	with open(args.outSentences,'w') as outF:
		headers = 'matching_id,pmid,title,journal,journal_short,year,month,day,section,subsection,role,predictprob,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,sentence,formatted_sentence'
		headerCount = len(headers.split(','))
		outF.write(headers.replace(',','\t') + '\n')

		sentenceOrdering = sorted(list(sentenceOrdering),reverse=True)

		for _,sentenceKey in sentenceOrdering:
			matchingID = sentenceData[sentenceKey]

			pmid,title,journal,journal_short,year,month,day,section,subsection,role,predictprob,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,sentence = sentenceKey

			charByChar = list(sentence)
			for start,end in geneLocs[sentenceKey]:
				charByChar[start] = '<b>' + charByChar[start]
				charByChar[end-1] += '</b>'
			for start,end in cancerLocs[sentenceKey]:
				charByChar[start] = '<b>' + charByChar[start]
				charByChar[end-1] += '</b>'

			formattedSentence = "".join(charByChar)

			outData = [ matchingID, pmid,title,journal,journal_short,year,month,day,section,subsection,role,predictprob,cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized,sentence, formattedSentence ]
			assert len(outData) == headerCount

			outLine = "\t".join(outData)
			outF.write(outLine + "\n")


	print("%d records filtered to %d sentences and collated to %d cancer gene roles" % (readCount, len(sentenceData), len(collated)))
	print("Written to %s and %s" % (args.outSentences, args.outCollated))

