# Pseudocode and Explanation for CancerMine

This file gives an overview for the following core files. Most of the functionality is managed by the [Kindred Python package](https://github.com/jakelever/kindred) which is described in the [documentation](http://kindred.readthedocs.io) and [associated paper](http://aclweb.org/anthology/W17-2322).

- **buildModels.py** : Use training data to build Kindred relation classifier models
- **wordlistLoader.py** : Preprepare parsed wordlists using gene names, cancer types and more
- **parseAndFindEntities.py** : Parse documents and find sentences that mention cancer types and gene names with additional filtering
- **applyModelsToSentences.py** : Apply Kindred relation classifiers to find mentions of drivers, oncogenes and tumor suppressors
- **filterAndCollate.py** : Filter for mentions with higher certainty and collate them for counts

## buildModels.py

- For each relation type (Driver, Oncogene, Tumor_Suppressor)
   - Load the Kindred corpus (1500 annotated sentences)
   - Strip all relations that do not match the relation type of interest
   - Create a Kindred classifier with a logistic regression model and threshold of 0.5
   - Train it on the filtered corpus
   - Save the classifier to a file

## wordlistLoader.py

- Take in wordlists for genes, cancers, drugs (to identify ambigiuity) and conflicting terms
- Get Kindred to parse them and prepare a data structure ready for matching
- Save it to a file as a Python pickle

## parseAndFindEntities.py

- Create a Kindred parser and EntityRecognizer with the terms prepared by wordlistLoader.py
- Read in a BioC corpus file (of abstracts or articles) in chunks:
   - Filter the corpus by removing documents that don't contain keywords (in filterTerms.txt)
   - Parse the documents
   - Annotate them with the EntityRecognizer (for cancer types, genes, etc)
   - For each sentence in each document
      - Ignore if it doesn't contain any of the keywords (in filterTerms.txt)
      - Check if a gene and cancer are mentioned in the sentence and add to output with metadata if so
- Dump all matching sentences to output JSON with metadata of the source of the sentence

## applyModelsToSentences.py

- Load all the models created by buildModels.py
- Create a Kindred parser and EntityRecognizer with the terms prepared by wordlistLoader.py
- Open the JSON file with sentences
- Parse them and annotated with EntityRecognizer (for cancer types, genes, etc)
- Apply the Kindred relation classifier models to this corpus
- Iterate over every relation extracted
    - Normalize gene names and cancer names where possible
    - Output the relation with all metadata and normalized terms

## filterAndCollate.py

- Define thresholds for relations:
   - Driver = 0.80
   - Oncogene = 0.76
   - Tumor_Suppressor = 0.92
- Iterate over the combined outputs of all runs of applyModelsToSentences.py
   - Check the probability of the relation (as the output for the model) and see if it is above the required thresholds
   - Get the core relation info of relation type, cancer type and gene name
   - Create a matching ID key that can link back this core relation info
   - Add the PubMed ID to the number of citations for this core relation info
- Output all the core relations with the number of citations
- Output all the sentences with an additional field of the sentence with simple HTML formatting for the location of entities

