#!/bin/bash
set -eux

# Some clean up
rm -fr exampledata/intermediate
rm -f exampledata/mini_terms.pickle exampledata/intermediate_sentences.json 
rm -f exampledata/out_unfiltered.tsv exampledata/out_collated.tsv exampledata/out_sentences.tsv
mkdir exampledata/intermediate

sh buildModelsIfNeeded.sh

python wordlistLoader.py --genes exampledata/mini_terms_genes.tsv --cancers exampledata/mini_terms_cancers.tsv --drugs exampledata/mini_terms_drugs.tsv --conflicting exampledata/mini_terms_conflicting.tsv --wordlistPickle exampledata/mini_terms.pickle

python parseAndFindEntities.py --biocFile exampledata/input.bioc.xml --filterTerms filterTerms.txt --wordlistPickle exampledata/mini_terms.pickle --outSentencesFilename exampledata/intermediate_sentences.json

python applyModelsToSentences.py --models models/cancermine.driver.model,models/cancermine.oncogene.model,models/cancermine.tumorsuppressor.model --filterTerms filterTerms.txt --wordlistPickle exampledata/mini_terms.pickle --genes exampledata/mini_terms_genes.tsv --cancerTypes exampledata/mini_terms_cancers.tsv --sentenceFile exampledata/intermediate_sentences.json --outData exampledata/intermediate/pubmed_0000.tsv

python filterAndCollate.py --inData exampledata/intermediate --outUnfiltered exampledata/out_unfiltered.tsv --outCollated exampledata/out_collated.tsv --outSentences exampledata/out_sentences.tsv

