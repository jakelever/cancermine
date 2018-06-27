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
		gene_names = { gene_name:math.log10(len(citations)+1) for gene_name,citations in geneNameWithCitations.items() }
		total = max(gene_names.values())
		gene_names = { gene_name:score/total for gene_name,score in gene_names.items() }

		gene_names_defaultDict = Counter()
		gene_names_defaultDict.update(gene_names)
		cancerGeneLists[cancer_name] = gene_names_defaultDict

	with open(args.tcgaData) as f:
		for line in f:
			gene,sample = line.strip('\n').split('\t')

			samples[sample].add(gene)
			genes.add(gene)

	cancerFilter = ['colorectal cancer','breast cancer','hepatocellular carcinoma','prostate cancer','lung cancer','malignant glioma','stomach cancer']
	cancerGeneLists = { cancer:geneList for cancer,geneList in cancerGeneLists.items() if cancer in cancerFilter }

	def scoreIt(sampleGenes,fingerprint,allGenes):
		return sum( fingerprint[g] for g in sampleGenes )

	for sampleID in samples:
		scores = [ (scoreIt(samples[sampleID],cancerGeneList,genes),cancerType) for cancerType,cancerGeneList in cancerGeneLists.items() ]
		scores = sorted(scores,reverse=True)

		if scores[0][0] == 0:
			score,bestCancerType = -1,"X no TS"
		elif scores[0][0] > scores[1][0]:
			score,bestCancerType = scores[0]
		else:
			score,bestCancerType = -1,"unclear"

		print("%s\t%s\t%f" % (sampleID,bestCancerType,score))

		sys.stdout.flush()


