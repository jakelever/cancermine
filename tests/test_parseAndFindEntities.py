#from .. import applyModel
import os
import sys

def test_basicrun():
	thisDir = os.path.dirname(os.path.abspath(__file__))
	parentDir = os.path.dirname(thisDir)
	sys.path.append(parentDir)
	import parseAndFindEntities

	biocFile = 'tests/testDoc.bioc.xml'
	filterTerms = 'filterTerms.txt'
	wordlistPickle = 'terms.pickle'
	outJSON = 'tests/testDoc.bioc.json'

	os.chdir(parentDir)
	parseAndFindEntities.parseAndFindEntities(biocFile,filterTerms,wordlistPickle,outJSON)

if __name__ == '__main__':
	test_basicrun()

