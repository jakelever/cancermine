#!/bin/bash
set -eux

if [ ! -d models ]; then
	mkdir models

	mkdir cancermine_corpus
	cd cancermine_corpus
	unzip -j -qq ../data/cancermine_corpus.zip
	cd -

	python buildModels.py --inTrain cancermine_corpus --outModel_Driver models/cancermine.driver.model --outModel_Oncogene models/cancermine.oncogene.model --outModel_TumorSuppressor models/cancermine.tumorsuppressor.model
fi

