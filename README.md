# CancerMine

<p>
<a href="https://travis-ci.org/jakelever/cancermine">
   <img src="https://travis-ci.org/jakelever/cancermine.svg?branch=master" />
</a>
<a href="http://bionlp.bcgsc.ca/cancermine/">
   <img src="https://img.shields.io/badge/data-viewer-9e42f4.svg" />
</a>
<a href="https://doi.org/10.5281/zenodo.1156241">
   <img src="https://zenodo.org/badge/DOI/10.5281/zenodo.1156241.svg" />
</a>
<a href="https://doi.org/10.1101/364406">
   <img src="https://img.shields.io/badge/bioRxiv-preprint-67baea.svg" />
</a>
</p>

The CancerMine resource is a text-mined knowledgebase of drivers, oncogenes and tumor suppressors in cancer. Abstracts from PubMed and full-text articles from PubMed Central Open Access subset and Author Manuscript Collections are processed to find references to genes as drivers, oncogenes and tumor suppressors in different cancer types.

## Using CancerMine

CancerMine is an automatically updated dataset. You can navigate the data using the [web viewer](http://bionlp.bcgsc.ca/cancermine/) or you can download the latest data from [Zenodo](https://doi.org/10.5281/zenodo.1156241) (or through the [web viewer](http://bionlp.bcgsc.ca/cancermine/)). You likely would not have to run any of the code in this repository.

## Installing and Running The Code

This project relies on text mining using [Kindred](https://github.com/jakelever/kindred) and resource management with [PubRunner](https://github.com/jakelever/pubrunner). These can be installed through pip. Remember to install the English language model for Spacy.

```
pip install kindred pubrunner
python -m spacy download en
```

To do a test run and check everything is okay, run:

```
pubrunner --test .
```

Then to do the full run which may take a long time, run:
```
pubrunner .
```

This will download all the corpora files, build and apply models. PubRunner can be setup to use a cluster (using SnakeMake). This is highly recommended.

## Inputs

The text inputs for processing are:

 - [PubMed](https://www.nlm.nih.gov/databases/download/pubmed_medline.html) abstracts
 - [PubMed Central Open Access subset (PMCOA)](https://www.ncbi.nlm.nih.gov/pmc/tools/openftlist/) full-text articles
 - [PubMed Central Author Manuscript Collection (PMCAMC)](https://www.ncbi.nlm.nih.gov/pmc/about/mscollection/) full-text articles

The text is scanned for references of genes and cancer types. These are based on [HUGO gene names](http://genenames.org/) and cancer types from the [Disease Ontology](http://www.disease-ontology.org/). These are managed through the [BioWordlists project](https://github.com/jakelever/biowordlists) which can be downloaded at https://doi.org/10.5281/zenodo.1286661.

The training data used to build the machine learning models can be found at [data/cancermine_corpus.zip](https://github.com/jakelever/cancermine/blob/master/data/cancermine_corpus.zip). This is stored in [BioNLP Shared Task format](http://2011.bionlp-st.org/home/file-formats) and has one file per sentence. The raw annotations from the three annotators can be found at [data/raw_annotations](https://github.com/jakelever/cancermine/tree/master/data/raw_annotations) 

## Outputs

There are three final results files from CancerMine. These are hosted on [Zenodo](https://doi.org/10.5281/zenodo.1156241) and can also be downloaded through the web viewer. Each file is a tab-delimited file with a header, no comments and no quoting.

You likely want **cancermine\_collated.tsv** if you just want the list of cancer gene roles. If you want the supporting sentences, look at **cancermine\_sentences.tsv**. You can use the *matching\_id* column to connect the two files. If you want to dig further and are okay with a higher false positive rate, look at **cancermine\_unfiltered.tsv**.

**cancermine\_collated.tsv:** This contains the cancer gene roles with citation counts supporting them. It contains the normalized cancer and gene names along with IDs for HUGO, Entrez Gene and the Disease Ontology.

**cancermine\_sentences.tsv:** This contains the supporting sentences for the cancer gene roles in the collated file. Each row is a single supporting sentence for one cancer gene role. This file contains information on the source publication (e.g. journal, publication date, etc), the actual sentence and the cancer gene role extracted.

**cancermine\_unfiltered.tsv:** This is the raw output of the applyModelsToSentences.py script across all of PubMed, Pubmed Central Open Access and PubMed Central Author Manuscript Collection. It contains every predicted relation with a prediction score above 0.5. So this may contain many false positives. Each row contain information on the publication (e.g. journal, publication date, etc) along with the sentence and the specific cancer gene role extracted (with HUGO, Entrez Gene and Disease Ontology IDs). This file is further processed to create the other two.

## Shiny App

The code in [shiny/](https://github.com/jakelever/cancermine/tree/master/shiny) is the Shiny code used for the [web viewer](http://bionlp.bcgsc.ca/cancermine/). If it is helpful, please use the code for your own projects. The list of dependencies is found at the top of the [app.R](https://github.com/jakelever/cancermine/blob/master/shiny/app.R) file.

## Paper

The code to generate all the figures and text for the paper can be found in [paper/](https://github.com/jakelever/cancermine/tree/master/paper). This may be useful for generating an up-to-date version of the plots for a newer version of CancerMine.

## Citing Us

The paper is now up at [bioRxiv](https://doi.org/10.1101/364406). It will be submitted to a journal in due course. It'd be wonderful if you would cite the paper if you use the methods or data set.

```
@article{lever2018cancermine,
  title={CancerMine: A literature-mined resource for drivers, oncogenes and tumor suppressors in cancer},
  author={Lever, Jake and Zhao, Eric Y and Grewal, Jasleen and Jones, Martin R and Jones, Steven JM},
  journal={bioRxiv},
  pages={364406},
  year={2018},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Issues

If you encounter any problems, please [file an issue](https://github.com/jakelever/cancermine/issues) along with a detailed description.
