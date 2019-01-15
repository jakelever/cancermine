import argparse

if __name__ == '__main__':
	parser = argparse.ArgumentParser('Combine CancerMine data with genome coordinates from a downloaded NCBI gene list file')
	parser.add_argument('--inCancerMineCollated',type=str,required=True,help='cancermine_collated.tsv file for input')
	parser.add_argument('--genomeAnnotationGFF',type=str,required=True,help='GFF3 file from ftp://ftp.ncbi.nih.gov/genomes/H_sapiens/GFF/ (unzipped)')
	parser.add_argument('--outCancerMineWithCoordinates',type=str,required=True,help='Output file with added coordinates')
	args = parser.parse_args()

	entrezID_to_coordinates = {}

	with open(args.genomeAnnotationGFF) as f:
		headers = ['seqid','source','type','start','end','score','strand','phase','attributes']
		for line in f:
			if line.startswith('#'):
				continue
			values = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }

			if values['type'] == 'gene' and values['seqid'].startswith('NC_000'):
				chrom = values['seqid'].replace('NC_','').split('.')[0]
				if chrom == '23':
					chrom = 'X'
				elif chrom == '24':
					chrom = 'Y'
				else:
					chrom = str(int(chrom))

				start = values['start']
				end = values['end']
				attributes = [ tuple(a.split('=',1)) for a in values['attributes'].split(';') ]
				attributes = { name:val for name,val in attributes }

				geneIDs = [ tuple(gid.split(':',1)) for gid in attributes['Dbxref'].split(',') ]
				geneIDs = { name:val for name,val in geneIDs }

				entrezID = geneIDs['GeneID']
				entrezID_to_coordinates[entrezID] = (chrom,start,end)

	with open(args.inCancerMineCollated) as inF, open(args.outCancerMineWithCoordinates,'w') as outF:
		headers = inF.readline().strip('\n').split('\t')
		newheaders = headers + ['gene_chrom','gene_start_coord','gene_end_coord']
		outF.write("\t".join(newheaders) + "\n")

		for line in inF:
			values = { h:v for h,v in zip(headers,line.strip('\n').split('\t')) }
			entrezID = values['gene_entrez_id']
			if entrezID in entrezID_to_coordinates:
				chrom,start,end = entrezID_to_coordinates[entrezID]
			else:
				chrom,start,end = 'NA','NA','NA'

			newline = "%s\t%s\t%s\t%s\n" % (line.strip('\n'),chrom,start,end)
			outF.write(newline)
