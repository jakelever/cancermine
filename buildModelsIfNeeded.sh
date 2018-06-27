#!/bin/bash
set -eux

if [ ! -d models ]; then
	mkdir models
	unzip -qq data/cancermine_corpus.zip
	mv cancermine_corpus/train/ cancermine_corpus/combined/
	mv cancermine_corpus/test/* cancermine_corpus/combined/

	python buildModels.py --inTrain cancermine_corpus/combined/ --outModel_Driver model/cancermine.driver.model --outModel_Oncogene model/cancermine.oncogene.model --outModel_TumorSuppressor model/cancermine.tumorsuppressor.model
fi

