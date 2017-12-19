#from .. import applyModel
import os
import sys

def test_basicrun():
	thisDir = os.path.dirname(os.path.abspath(__file__))
	parentDir = os.path.dirname(thisDir)
	sys.path.append(parentDir)
	import applyModel

	biocFile = 'tests/testDoc.bioc.xml'
	inModel_Driver = 'cancermine.driver.model'
	inModel_Oncogene = 'cancermine.oncogene.model'
	inModel_TumorSuppressor = 'cancermine.tumorsuppressor.model'
	filterTerms = 'filterTerms.txt'
	wordlistPickle = 'terms.pickle'
	genes = 'terms_genes.txt'
	cancerTypes = 'terms_cancer.txt'
	outData = 'tests/testDoc.bioc.out'

	os.chdir(parentDir)
	applyModel.cancermine(biocFile,inModel_Driver,inModel_Oncogene,inModel_TumorSuppressor,filterTerms,wordlistPickle,genes,cancerTypes,outData)

if __name__ == '__main__':
	test_basicrun()

