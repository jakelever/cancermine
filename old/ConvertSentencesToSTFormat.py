import sys
import argparse
import codecs
import json

def attemptIntegerConversion(text):
	try:
		number = int(text)
		return number
	except:
		pass
	return None

def integerOrNone(text):
	if len(text) == 0:
		return None
	else:
		return int(text)

def file_len(fname):
	with open(fname) as f:
		for i, l in enumerate(f):
			pass
	return i + 1

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Converts the output of the SentenceSelector to the ST format expected by the annotator and by VERSE')
	parser.add_argument('--sentencesFile',type=str,required=True,help='File containing sentences with PMID and entity information')
	parser.add_argument('--outDir',type=str,required=True,help='Output directory to save converted files')

	args = parser.parse_args()

	outDir = args.outDir
	if not outDir.endswith('/'):
		outDir += '/'

	totalLineCount = file_len(args.sentencesFile)
	nextPerc = 0.0

	with codecs.open(args.sentencesFile,'r','utf-8') as inFile:
		for sentenceid,line in enumerate(inFile):
			perc = 100.0 * sentenceid / float(totalLineCount)
			if perc > nextPerc:
				print "%.1f%%" % perc
				nextPerc += 1.0


			split = line.rstrip().split('\t')
			assert len(split) > 3
			pmid = split[0]
			pmcid = split[1]
			text = split[2]
			
			base = "%s%08d" % (outDir,sentenceid)
			txtFilename = base + '.txt'
			a1Filename = base + '.a1'
			jsonFilename = base + '.json'

			with codecs.open(txtFilename,'w','utf-8') as txtFile:
				txtFile.write(text)

			#print text
			#print split[3:]
			jsonData = {}
			jsonData['pmid'] = integerOrNone(pmid)
			jsonData['pmcid'] = integerOrNone(pmcid)
			jsonData['text'] = text
			jsonData['entities'] = {}
			
			with codecs.open(a1Filename,'w','utf-8') as a1File:
				for entityid,entityData in enumerate(split[3:]):
					entitytype,entityids,start,end,entitytxt = entityData.split('|')
					entityidtxt = "T%d" % (entityid+1)
					start,end = int(start),int(end)
					#print start,end, text[start:end]
					assert entitytxt == text[start:end], u"%s != %s (%s)" % (entitytxt, text[start:end], base)
					line = u"%s\t%s %d %d\t%s\n" % (entityidtxt,entitytype,start,end,entitytxt)
					a1File.write(line)

					entityDict = {}
					entityDict['entitytype'] = entitytype
					entityDict['entityids'] = entityids
					entityDict['start'] = start
					entityDict['end'] = end
					entityDict['entitytxt'] = entitytxt

					jsonData['entities'][entityidtxt] = entityDict



			with open(jsonFilename,'w') as jsonFile:
				json.dump(jsonData,jsonFile,indent=4)

