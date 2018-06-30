This describes the output files for the [CancerMine](https://github.com/jakelever/cancermine) project. These files are loaded directly by the [CancerMine viewer](http://bionlp.bcgsc.ca/cancermine/). The code for this viewer is available in the CancerMine Github repo if you want to run it independently. Each file is a tab-delimited file with a header, no comments and no quoting.

You likely want **cancermine\_collated.tsv** if you just want the list of cancer gene roles. If you want the supporting sentences, look at **cancermine\_sentences.tsv**. You can use the *matching\_id* column to connect the two files. If you want to dig further and are okay with a higher false positive rate, look at **cancermine\_unfiltered.tsv**.

**cancermine\_collated.tsv:** This contains the cancer gene roles with citation counts supporting them. It contains the normalized cancer and gene names along with IDs for HUGO, Entrez Gene and the Disease Ontology.

**cancermine\_sentences.tsv:** This contains the supporting sentences for the cancer gene roles in the collated file. Each row is a single supporting sentence for one cancer gene role. This file contains information on the source publication (e.g. journal, publication date, etc), the actual sentence and the cancer gene role extracted.

**cancermine\_unfiltered.tsv:** This is the raw output of the applyModelsToSentences.py script across all of PubMed, Pubmed Central Open Access and PubMed Central Author Manuscript Collection. It contains every predicted relation with a prediction score above 0.5. So this may contain many false positives. Each row contain information on the publication (e.g. journal, publication date, etc) along with the sentence and the specific cancer gene role extracted (with HUGO, Entrez Gene and Disease Ontology IDs). This file is further processed to create the other two.

