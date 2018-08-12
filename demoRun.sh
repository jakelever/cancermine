#!/bin/bash
set -eux

sh buildModelsIfNeeded.sh

python wordlistLoader.py --genes exampledata/mini_terms_genes.tsv --cancers exampledata/mini_terms_cancers.tsv --drugs exampledata/mini_terms_drugs.tsv --conflicting exampledata/mini_terms_conflicting.tsv --wordlistPickle exampledata/mini_terms.pickle

python parseAndFindEntities.py --biocFile exampledata/input.bioc.xml --filterTerms filterTerms.txt --wordlistPickle exampledata/mini_terms.pickle --outSentencesFilename exampledata/intermediate_sentences.json

python applyModelsToSentences.py --models models/cancermine.driver.model,models/cancermine.oncogene.model,models/cancermine.tumorsuppressor.model --filterTerms filterTerms.txt --wordlistPickle exampledata/mini_terms.pickle --genes exampledata/mini_terms_genes.tsv --cancerTypes exampledata/mini_terms_cancers.tsv --sentenceFile exampledata/intermediate_sentences.json --outData exampledata/intermediate_relations.json

cat header.tsv exampledata/intermediate_relations.json > exampledata/out_unfiltered.tsv

python filterAndCollate.py --inUnfiltered exampledata/out_unfiltered.tsv --outCollated exampledata/out_collated.tsv --outSentences exampledata/out_sentences.tsv

