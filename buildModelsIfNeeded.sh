#!/bin/bash

if [ ! -d models ]; then
	mkdir models
	unzip data/cancermine_corpus.zip
	python buildModels.py --inTrain cancermine_corpus/combined/ --outModel_Driver model/cancermine.driver.model --outModel_Oncogene model/cancermine.oncogene.model --outModel_TumorSuppressor model/cancermine.tumorsuppressor.model
fi

