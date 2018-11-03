## @knitr fig4

source('cancermine/dependencies.R')


collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)

unfilteredFilename <- 'cancermine/cancermine_unfiltered.tsv'
unfilteredFilename <- normalizePath(unfilteredFilename)


cancermine <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
#cancermine <- cancermine[order(cancermine$citation_count,decreasing=T),]

cancermineUnfiltered <- read.table(unfilteredFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


cgcData <- read.table('cancermine/cgc_mapped.tsv',sep='\t')
colnames(cgcData) <- c('role','gene_entrez_id','gene_normalized','cancer_id')

paper.cgcGeneCount <- length(unique(cgcData$gene_entrez_id))
paper.cgcCancerCount <- length(unique(cgcData$cancer_id))

paper.cgcGeneCount <- prettyNum(paper.cgcGeneCount,big.mark=',')
paper.cgcCancerCount <- prettyNum(paper.cgcCancerCount,big.mark=',')

#cancermine$cancer_id <- as.character(cancermine$cancer_id)
#cgcData$cancer_id <- as.character(cgcData$cancer_id)
#cancermine <- cancermine[cancermine$role!='Driver',]
#cancermine[cancermine$cancer_id == 'DOID:3908','cancer_id'] <- 'DOID:1324'
#cgcData[cgcData$cancer_id == 'DOID:3908','cancer_id'] <- 'DOID:1324'

cancermine$combined <- paste(cancermine$role,as.character(cancermine$gene_entrez_id),cancermine$cancer_id, sep='_')
cgcData$combined <- paste(cgcData$role,cgcData$gene_entrez_id,cgcData$cancer_id,sep='_')

cancermineUnfiltered$combined <- paste(cancermineUnfiltered$role,as.character(cancermineUnfiltered$gene_entrez_id),cancermineUnfiltered$cancer_id, sep='_')

#cancermine$inCGC <- cancermine$combined %in% cgcData$combined
#missing <- cancermine[!cancermine$inCGC,]

cgcNotInCancerMine <- cgcData[!(cgcData$combined %in% cancermine$combined),]
paper.cgcFoundInLessConservative <- length(unique(cgcNotInCancerMine$combined[cgcNotInCancerMine$combined %in% cancermineUnfiltered$combined]))
#missing <- missing[!duplicated(missing),]

head(cgcData[!(cgcData$combined %in% cancermineUnfiltered$combined),])

matching <- cancermine[cancermine$combined %in% cgcData$combined,]
matching <- matching[order(matching$citation_count),]

noDrivers <- cancermine[cancermine$role!='Driver',]
paper.moreGenesThanCGC <- length(setdiff(noDrivers$gene_entrez_id,cgcData$gene_entrez_id))
paper.moreGenesThanCGC <- prettyNum(paper.moreGenesThanCGC,big.mark=',')


listInput <- list(CancerMine=unique(noDrivers$combined),CGC=unique(cgcData$combined))
cgcPlot <- venn.diagram(
  x = listInput,
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
cgcPlot <- gTree(children=cgcPlot)
cgcPlot <- grid.arrange(cgcPlot,top='(a)')





intogen2014 <- read.table('cancermine/Mutational_drivers_per_tumor_type.tsv',header=T,sep='\t')
intogen2014$Tumor_type <- as.character(intogen2014$Tumor_type)
intogen2014[intogen2014$Tumor_type=='ALL','Tumor_type'] <- 'acute lymphocytic leukemia'
intogen2014[intogen2014$Tumor_type=='GBM','Tumor_type'] <- 'glioblastoma multiforme'
intogen2014[intogen2014$Tumor_type=='LUSC','Tumor_type'] <- 'Lung squamous cell carcinoma'
intogen2014[intogen2014$Tumor_type=='UCEC','Tumor_type'] <- 'endometrial adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='HNSC','Tumor_type'] <- 'Head and neck squamous cell carcinoma'
intogen2014[intogen2014$Tumor_type=='LUAD','Tumor_type'] <- 'Lung adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='STAD','Tumor_type'] <- 'gastric adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='BRCA','Tumor_type'] <- 'Breast cancer'
intogen2014[intogen2014$Tumor_type=='CM','Tumor_type'] <- 'melanoma'
intogen2014[intogen2014$Tumor_type=='COREAD','Tumor_type'] <- 'Colorectal adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='ESCA','Tumor_type'] <- 'Esophageal carcinoma'
intogen2014[intogen2014$Tumor_type=='LGG','Tumor_type'] <- 'malignant glioma'
intogen2014[intogen2014$Tumor_type=='OV','Tumor_type'] <- 'ovary adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='RCCC','Tumor_type'] <- 'clear cell renal cell carcinoma'
intogen2014[intogen2014$Tumor_type=='BLCA','Tumor_type'] <- 'Bladder carcinoma'
intogen2014[intogen2014$Tumor_type=='DLBC','Tumor_type'] <- 'diffuse large B-cell lymphoma'
intogen2014[intogen2014$Tumor_type=='CLL','Tumor_type'] <- 'Chronic lymphocytic leukemia'
intogen2014[intogen2014$Tumor_type=='PAAD','Tumor_type'] <- 'pancreatic adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='HC','Tumor_type'] <- 'liver carcinoma'
intogen2014[intogen2014$Tumor_type=='PRAD','Tumor_type'] <- 'Prostate adenocarcinoma'
intogen2014[intogen2014$Tumor_type=='SCLC','Tumor_type'] <- 'lung small cell carcinoma'
intogen2014[intogen2014$Tumor_type=='THCA','Tumor_type'] <- 'Thyroid carcinoma'
intogen2014[intogen2014$Tumor_type=='NB','Tumor_type'] <- 'Neuroblastoma'
intogen2014[intogen2014$Tumor_type=='NSCLC','Tumor_type'] <- 'non-small cell lung carcinoma'
intogen2014[intogen2014$Tumor_type=='MM','Tumor_type'] <- 'Multiple myeloma'
intogen2014[intogen2014$Tumor_type=='MB','Tumor_type'] <- 'medulloblastoma'
intogen2014[intogen2014$Tumor_type=='AML','Tumor_type'] <- 'Acute myeloid leukemia'
intogen2014[intogen2014$Tumor_type=='PA','Tumor_type'] <- 'Pylocytic astrocytoma'

intogen2014$Tumor_type <- tolower(intogen2014$Tumor_type)

paper.geneCancerCount <- length(unique(intogen2014$Tumor_type))
paper.intogenGeneCount <- length(unique(intogen2014$geneHGNCsymbol))

paper.geneCancerCount <- prettyNum(paper.geneCancerCount,big.mark=',')
paper.intogenGeneCount <- prettyNum(paper.intogenGeneCount,big.mark=',')

intogenGeneCount <- plyr::count(intogen2014$geneHGNCsymbol)
intogenGeneCount <- intogenGeneCount[order(intogenGeneCount$freq,decreasing=T),]
head(intogenGeneCount)
unique(cancermine[cancermine$gene_normalized=='ARID1A' & cancermine$role=='Driver','combined'])

#cancermine <- read.table(cancermineFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
#cancermine <- cancermine[grep(";",cancermine$cancer_normalized,fixed=T,invert=T),]
#cancermine <- cancermine[grep(";",cancermine$gene_normalized,fixed=T,invert=T),]
#cancermine <- cancermine[grep("|",cancermine$gene_normalized,fixed=T,invert=T),]

#cancermap <- cancermine[,c('cancer_name','cancer_normalized')]
#cancermap$cancer_name <- tolower(cancermap$cancer_name)
#cancermap <- cancermap[!duplicated(cancermap),]
#rownames(cancermap) <- cancermap$cancer_name

#cancermineCancers <- tolower(unique(cancermine$cancer_normalized))
#intogenCancers <- tolower(unique(intogen2014$Tumor_type))
#intogen2014$cancer_normalized <- cancermap[tolower(intogen2014$Tumor_type),'cancer_normalized']
#mismatch <- data.frame(mismatch=intogenCancers[!intogenCancers %in% cancermineCancers])
#mismatch$fix <- cancermap[tolower(mismatch$mismatch),'cancer_normalized']

#intogen2014$role <- 'Driver'

cancermine$combined <- tolower(paste(as.character(cancermine$gene_normalized),cancermine$cancer_normalized, sep='_'))
intogen2014$combined <- tolower(paste(as.character(intogen2014$geneHGNCsymbol),intogen2014$Tumor_type, sep='_'))

intersection <- unique(cancermine$combined[cancermine$combined %in% intogen2014$combined])
uniqueCombo <- unique(intogen2014$combined[!intogen2014$combined %in% cancermine$combined])
uniqueGenes <- unique(intogen2014$geneHGNCsymbol[!intogen2014$geneHGNCsymbol %in% cancermine$gene_normalized])

#grep("arid1a", intersection, value=T)
#grep("arid1a", uniqueCombo, value=T)

paper.intogenGenesInCancermine <- sum(unique(intogen2014$geneHGNCsymbol) %in% unique(cancermine$gene_normalized))
paper.intogenGenesInCancermine <- prettyNum(paper.intogenGenesInCancermine,big.mark=',')

#cancermine <- cancermine[tolower(cancermine$cancer_normalized) %in% tolower(intogen2014$Tumor_type),]

#missing <- intogen2014[!(intogen2014$combined %in% cancermine$combined),]

listInput <- list(CancerMine=unique(cancermine$combined),IntOGen=unique(intogen2014$combined))
intogenPlot <- venn.diagram(
  x = listInput,
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
intogenPlot <- gTree(children=intogenPlot)
intogenPlot <- grid.arrange(intogenPlot,top='(b)')

paper.moreGenesThanIntogen <- length(setdiff(unique(cancermine$gene_normalized),unique(intogen2014$geneHGNCsymbol)))
paper.moreGenesThanIntogen <- prettyNum(paper.moreGenesThanIntogen,big.mark=',')

#length(unique(onlyDrivers[onlyDrivers$role=='Driver','gene_normalized']))
#length(unique(cancermine[cancermine$role=='Driver','gene_normalized']))

cancermineTSIDs <- cancermine[cancermine$role=='Tumor_Suppressor','gene_entrez_id',drop=F]
tsgeneData <- read.table('cancermine/TSgene.txt',header=T,sep='\t')

listInput <- list(CancerMine=unique(cancermineTSIDs$gene_entrez_id),TSGene=unique(tsgeneData$GeneID))
tsgenePlot <- venn.diagram(
  x = listInput,
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
tsgenePlot <- gTree(children=tsgenePlot)
tsgenePlot <- grid.arrange(tsgenePlot,top='(c)')

paper.tsgeneCount <- length(unique(tsgeneData$GeneID))
paper.tsgeneOverlapCount <- sum(unique(tsgeneData$GeneID) %in% unique(cancermineTSIDs$gene_entrez_id))
paper.tsgeneOverlapPerc <- round(100.0 * paper.tsgeneOverlapCount / paper.tsgeneCount)

cancermineOncogeneIDs <- cancermine[cancermine$role=='Oncogene','gene_entrez_id',drop=F]
ongeneData <- read.table('cancermine/ongene.txt',header=T,sep='\t',quote='')

listInput <- list(CancerMine=unique(cancermineOncogeneIDs$gene_entrez_id),ONGene=unique(ongeneData$OncogeneID))
ongenePlot <- venn.diagram(
  x = listInput,
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
ongenePlot <- gTree(children=ongenePlot)
ongenePlot <- grid.arrange(ongenePlot,top='(d)')

paper.ongeneCount <- length(unique(ongeneData$OncogeneID))
paper.ongeneOverlapCount <- sum(unique(ongeneData$OncogeneID) %in% unique(cancermineOncogeneIDs$gene_entrez_id))
paper.ongeneOverlapPerc <- round(100.0 * paper.ongeneOverlapCount / paper.ongeneCount)



civicDB <- read.table('cancermine/nightly-ClinicalEvidenceSummaries.tsv',header=T,sep='\t',quote='',comment='')
listInput <- list(CancerMine=unique(cancermine$gene_entrez_id),CIViC=unique(civicDB$entrez_id))
civicPlot <- venn.diagram(
  x = listInput,
  scaled=F,
  fill = c("grey", "white"),
  cat.fontface = 2,
  cat.pos = 0,
  filename=NULL)
civicPlot <- gTree(children=civicPlot)
civicPlot <- grid.arrange(civicPlot,top='(e)')

paper.civicCount <- length(unique(civicDB$entrez_id))
paper.civicOverlapCount <- sum(unique(civicDB$entrez_id) %in% unique(cancermine$gene_entrez_id))
paper.civicOverlapPerc <- round(100.0 * paper.civicOverlapCount / paper.civicCount)





fig_comparisons <- arrangeGrob(cgcPlot,intogenPlot,tsgenePlot,ongenePlot, civicPlot, ncol=2)

grid.arrange(fig_comparisons)

vennjunk <- dir(path=".", pattern="VennDiagram*") # ?dir
removed <- file.remove(vennjunk) # ?file.remove
