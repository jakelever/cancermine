import argparse
from collections import defaultdict,Counter
import sys
import random
import pronto
import math

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--inCancermine',required=True,type=str,help='Cancermine TSV to calculate inference scores for')
	parser.add_argument('--tcgaData',required=True,type=str,help='TCGA data fun')
	args = parser.parse_args()

	cancerGeneCitations = defaultdict(lambda : defaultdict(set))
	cancermineGenes = set()
	genes = set()

	cancerCounts = Counter()

	with open(args.inCancermine,'r') as f:
		headers = f.readline().strip('\n').split('\t')
		for line in f:
			lineDict = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			pmid = lineDict['pmid']
			role = lineDict['role']
			gene_id = lineDict['gene_hugo_id']
			cancer_id = lineDict['cancer_id']
			gene_name = lineDict['gene_normalized']
			cancer_name = lineDict['cancer_normalized']
			predictprob = float(lineDict['predictprob'])

			if role != 'Tumor_Suppressor':
				continue
			if gene_name in ['TP53']:
				continue

			if ";" in gene_name or '|' in gene_name:
				continue

			cancerGeneCitations[cancer_name][gene_name].add(pmid)
			cancermineGenes.add(gene_name)

			cancerCounts[cancer_name] += 1

	samples = defaultdict(set)

	cancerGeneLists = {}
	for cancer_name,geneNameWithCitations in cancerGeneCitations.items():
		#geneNameWithCitations = { gene_name:citations for gene_name,citations in geneNameWithCitations.items() if len(citations) >= 3 }
		#gene_names = [ gene_name for gene_name,citations in geneNameWithCitations.items() if len(citations) >= 3 ]
		#total_citations = sum( [ len(citations) for gene_name,citations in geneNameWithCitations.items() ] )
		#gene_names = { gene_name:len(citations)/total_citations for gene_name,citations in geneNameWithCitations.items() }
		gene_names = { gene_name:math.log10(len(citations)+1) for gene_name,citations in geneNameWithCitations.items() }
		#total = sum(gene_names.values())
		total = max(gene_names.values())
		gene_names = { gene_name:score/total for gene_name,score in gene_names.items() }
		#gene_names = { gene_name:len(citations) for gene_name,citations in geneNameWithCitations.items() }

		gene_names_defaultDict = Counter()
		gene_names_defaultDict.update(gene_names)
		cancerGeneLists[cancer_name] = gene_names_defaultDict

	with open(args.tcgaData) as f:
		for line in f:
			gene,sample = line.strip('\n').split('\t')

			#if gene in cancermineGenes:
			samples[sample].add(gene)
			genes.add(gene)

	cancerFilter = ['colorectal cancer','breast cancer','hepatocellular carcinoma','prostate cancer','lung cancer','malignant glioma','stomach cancer']
	cancerGeneLists = { cancer:geneList for cancer,geneList in cancerGeneLists.items() if cancer in cancerFilter }

	#groupCounter = Counter()
	#sampleGenes = [ gene for sampleID,sample in samples.items() for gene in sample ]
	#for cancerType,cancerGeneList in cancerGeneLists.items():
	#	count = len([ g for g in sampleGenes if g in cancerGeneList ])
	#	print("%s\t%d\t%d\t%f\t%d" % (cancerType,count,len(cancerGeneList),count/len(cancerGeneList),len(sampleGenes)))

	#fingerprint = cancerGeneLists['lung cancer']
	#print(sorted( [ (score,gene) for gene,score in fingerprint.items() ] ))
	#assert False
			
	#sys.exit(0)

	#print( [ (cancerType,len(cancerGeneList)) for cancerType,cancerGeneList in cancerGeneLists.items() ] )
	#assert False
	#for gene,score in cancerGeneLists['malignant glioma'].items():
	#	print("%s\t%f" % (gene,score))
	#assert False
	#sys.exit(0)


	def scoreIt(sampleGenes,fingerprint,allGenes):
		#overlapped = [ g for g in sampleGenes if g in cancerGenes ]
		#return len(set(overlapped))#/len(cancerGenes)
		return sum( fingerprint[g] for g in sampleGenes )
	#	selected = sorted([ (score,g) for g,score in fingerprint.items() ],reverse=True)
	#	selected = set([ g for score,g in selected if score > 0.001 ])
	#	return sum( 1 for g in sampleGenes if g in selected )

	for sampleID in samples:
		#for cancerType,cancerGeneList in cancerGeneLists.items():
		#	print("%s\t%s\t%f" % (sampleID,cancerType,enrich(samples[sampleID],cancerGeneList,genes)))

		#scores = [ [cancerType,overlap(samples[sampleID],cancerGeneList,genes)] for cancerType,cancerGeneList in cancerGeneLists.items() ]
		#scores = sorted(scores)
		#score,bestCancerType = scores[0]
		#scoresTxt = "\t".join(map(str,sum(scores,[])))
		#if score < 0.05:
		#print("%s\t%d\t%s" % (sampleID,len(samples[sampleID]),scoresTxt))

		scores = [ (scoreIt(samples[sampleID],cancerGeneList,genes),cancerType) for cancerType,cancerGeneList in cancerGeneLists.items() ]
		scores = sorted(scores,reverse=True)
		#print(scores)
		if scores[0][0] == 0:
			score,bestCancerType = -1,"X no TS"
		elif scores[0][0] > scores[1][0]:
			score,bestCancerType = scores[0]
		else:
			score,bestCancerType = -1,"unclear"
		#scoresTxt = "\t".join(map(str,sum(scores,[])))
		#if score < 0.05:
		print("%s\t%s\t%f" % (sampleID,bestCancerType,score))

		#for score,cancerType in scores:
		#	print("%s\t%s\t%s\t%f" % (args.tcgaData,sampleID,cancerType,score))
		
		sys.stdout.flush()


