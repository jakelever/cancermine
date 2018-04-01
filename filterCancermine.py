import argparse
import csv

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Filter Cancermine for more conservative predictions')
	parser.add_argument('--inCancermine',required=True,type=str,help='Cancermine TSV to filter')
	parser.add_argument('--outFiltered',required=True,type=str,help='Output filtered Cancermine')
	args = parser.parse_args()

	thresholds = {'Driver':0.80, 'Oncogene': 0.76, 'Tumor_Suppressor': 0.92}

	readCount,writeCount = 0,0
	with open(args.inCancermine,'r') as inF, open(args.outFiltered,'w') as outF:
		inTSV = csv.reader(inF,delimiter='\t')
		outTSV = csv.writer(outF,delimiter='\t')
		
		headers = next(inTSV, None)
		outTSV.writerow(headers)

		for row in inTSV:
			rowDict = { h:v for h,v in zip(headers,row) }
			predictprob = float(rowDict['predictprob'])
			role = rowDict['role']
			readCount += 1

			if predictprob > thresholds[role]:
				outTSV.writerow(row)
				writeCount += 1

	print("%d of %d rows written out to %s" % (writeCount,readCount,args.outFiltered))

