#!/bin/bash

biowordlistDir=$1

if [ ! -f cancermine_terms.pickle ]; then
	python wordlistLoader.py --genes $biowordlistDir/terms_genes.tsv --cancers $biowordlistDir/terms_cancers.tsv --drugs $biowordlistDir/terms_drugs.tsv --conflicting $biowordlistDir/terms_conflicting.tsv --wordlistPickle cancermine_terms.pickle
fi

