## @knitr fig_textminedOverview

source('cancermine/dependencies.R')

collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermineCollated <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


sentencesFilename <- 'cancermine/cancermine_sentences.tsv'
sentencesFilename <- normalizePath(sentencesFilename)
cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
cancermineSentences$journal_extra_short <- tolower(strtrim(cancermineSentences$journal,21))




paper.titleCount <- sum(cancermineSentences$section=='title')
paper.abstractCount <- sum(cancermineSentences$section=='abstract')
paper.articleCount <- sum(cancermineSentences$section=='article')

paper.geneDriverCount <- length(unique(cancermineSentences[cancermineSentences$role=='Driver','gene_normalized']))
paper.geneOncogeneCount <- length(unique(cancermineSentences[cancermineSentences$role=='Oncogene','gene_normalized']))
paper.geneTumorSuppressorCount <- length(unique(cancermineSentences[cancermineSentences$role=='Tumor_Suppressor','gene_normalized']))
paper.geneCount <- length(unique(cancermineSentences$gene_normalized))
paper.cancerCount <- length(unique(cancermineSentences$cancer_normalized))
paper.sentenceCount <- length(unique(cancermineSentences$sentence))
paper.pmidCount <-length(unique(cancermineSentences$pmid))

paper.driverSentenceCount <- length(unique(cancermineSentences[cancermineSentences$role=='Driver','sentence']))
paper.oncogeneSentenceCount <- length(unique(cancermineSentences[cancermineSentences$role=='Oncogene','sentence']))
paper.tumorSuppressorSentenceCount <- length(unique(cancermineSentences[cancermineSentences$role=='Tumor_Suppressor','sentence']))


#cancermineUniqueTriplesWithSection <- cancermine[,c('role','gene_normalized','cancer_normalized')]
#cancermineUniqueTriples <- cancermineUniqueTriples[!duplicated(cancermineUniqueTriples),]

paper.diseaseOntologyCancerCount <- 2044
paper.diseaseOntologyCoverage <- round(100*paper.cancerCount / paper.diseaseOntologyCancerCount)


paper.titleCount <- prettyNum(paper.titleCount,big.mark=',')
paper.abstractCount <- prettyNum(paper.abstractCount,big.mark=',')
paper.articleCount <- prettyNum(paper.articleCount,big.mark=',')

paper.geneDriverCount <- prettyNum(paper.geneDriverCount,big.mark=',')
paper.geneOncogeneCount <- prettyNum(paper.geneOncogeneCount,big.mark=',')
paper.geneTumorSuppressorCount <- prettyNum(paper.geneTumorSuppressorCount,big.mark=',')
paper.geneCount <- prettyNum(paper.geneCount,big.mark=',')
paper.cancerCount <- prettyNum(paper.cancerCount,big.mark=',')
paper.sentenceCount <- prettyNum(paper.sentenceCount,big.mark=',')
paper.pmidCount <- prettyNum(paper.pmidCount,big.mark=',')
paper.diseaseOntologyCancerCount <- prettyNum(paper.diseaseOntologyCancerCount,big.mark=',')


roleCounts <- plyr::count(cancermineCollated[,c('role'),drop=F])
rolePlot <- barchart(freq ~ role, 
                     ylab="Citation #",
                     ylim=c(0,1.1*max(roleCounts$freq)),
                         roleCounts, 
                         scales=list(x=list(rot=45)),
                         col="black")
rolePlot <- arrangeGrob(rolePlot,top='(a)')

topCount <- 20

geneCounts <- plyr::count(cancermineCollated[,c('gene_normalized'),drop=F])
geneCounts <- geneCounts[order(geneCounts$freq,decreasing=T),]
geneCounts$gene_normalized <- factor(as.character(geneCounts$gene_normalized), levels=as.character(geneCounts$gene_normalized))
topCancersPlot <- barchart(freq ~ gene_normalized, 
                           ylab="Citation #",
                           ylim=c(0,1.1*max(geneCounts$freq)),
         geneCounts[1:topCount,], 
         scales=list(x=list(rot=45)),
         col="black")
topCancersPlot <- arrangeGrob(topCancersPlot,top='(b)')

cancerCounts <- plyr::count(cancermineCollated[,c('cancer_normalized'),drop=F])
cancerCounts <- cancerCounts[order(cancerCounts$freq,decreasing=T),]
cancerCounts$cancer_normalized <- factor(as.character(cancerCounts$cancer_normalized), levels=as.character(cancerCounts$cancer_normalized))
topGenesPlot <- barchart(freq ~ cancer_normalized, 
                         ylab="Citation #",
                         ylim=c(0,1.1*max(cancerCounts$freq)),
         cancerCounts[1:topCount,], 
         scales=list(x=list(rot=45)),
         col="black")
topGenesPlot <- arrangeGrob(topGenesPlot,top='(c)')




hasFullText <- unique(cancermineSentences[cancermineSentences$section == 'article','journal_extra_short'])

journalCounts <- plyr::count(cancermineSentences[,c('journal_extra_short'),drop=F])
journalCounts <- journalCounts[order(journalCounts$freq,decreasing=T),]
journalCounts$journal_extra_short <- factor(as.character(journalCounts$journal_extra_short), levels=as.character(journalCounts$journal_extra_short))
#journalCounts$inPMCOA <- journalCounts$journal_extra_short %in% hasFullText
#journalCounts$inPMCOA_txt <- 'Journal in PMCOA'
#journalCounts$inPMCOA_txt[!journalCounts$inPMCOA] <- 'Journal NOT in PMCOA'

journalAndSectionCounts <- plyr::count(cancermineSentences[,c('journal_extra_short','section'),drop=F])
journalAndSectionCounts <- journalAndSectionCounts[journalAndSectionCounts$section %in% c('title','article','abstract'),]
journalAndSectionCounts$section <- factor(journalAndSectionCounts$section, labels=c('Abstract','Article','Title'))
journalAndSectionCounts$journal_extra_short <- factor(as.character(journalAndSectionCounts$journal_extra_short), levels=as.character(journalCounts$journal_extra_short))


topJournalAndSectionCounts <- journalAndSectionCounts[journalAndSectionCounts$journal_extra_short %in% journalCounts$journal_extra_short[1:topCount],]

myColours <- c(brewer.pal(4,"Set2"),"#000000")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

topJournalsPlot <- barchart(freq ~ journal_extra_short, 
         ylab="Citation #",
         ylim=c(0,1.1*max(journalCounts$freq)),
         topJournalAndSectionCounts, 
         scales=list(x=list(rot=45)),
         groups=section,
         par.settings = my.settings,
         auto.key=list(space="top", columns=2, 
                       points=FALSE, rectangles=TRUE),
         stack=T)
topJournalsPlot <- arrangeGrob(topJournalsPlot,top='(d)')
grid.arrange(topJournalsPlot)

citationCounts <- data.frame(group=c('1','2-4','5-9','10-19','20+'),
                             count=c(sum(cancermineCollated$citation_count==1),
                                     sum(cancermineCollated$citation_count>=2 & cancermineCollated$citation_count<=4),
                                     sum(cancermineCollated$citation_count>=5 & cancermineCollated$citation_count<=9),
                                     sum(cancermineCollated$citation_count>=10 & cancermineCollated$citation_count<=19),
                                     sum(cancermineCollated$citation_count>=20)
                                     ))
citationCounts$group <- factor(citationCounts$group,levels=citationCounts$group)
citationPlot <- barchart(count ~ group, 
                         citationCounts, 
                         ylim=c(0,1.1*max(citationCounts$count)),
                         xlab="Number of citations",
                         ylab="# of cancer gene roles",
                         col="black")
citationPlot <- arrangeGrob(citationPlot,top='(e)')


paper.bcrabl_oncogene_cml <- cancermineCollated[cancermineCollated$gene_normalized=='BCR|ABL1' & cancermineCollated$cancer_normalized=='chronic myeloid leukemia' & cancermineCollated$role=='Oncogene','citation_count']
paper.apc_tumorsuppressor_colorectal <- cancermineCollated[cancermineCollated$gene_normalized=='APC' & cancermineCollated$cancer_normalized=='colorectal cancer' & cancermineCollated$role=='Tumor_Suppressor','citation_count']

paper.erbb2_oncogene_breast <- cancermineCollated[cancermineCollated$gene_normalized=='ERBB2' & cancermineCollated$cancer_normalized=='breast cancer' & cancermineCollated$role=='Oncogene','citation_count']
paper.erbb2_oncogene_breast_firstYear <- min(cancermineSentences[cancermineSentences$gene_normalized=='ERBB2' & cancermineSentences$cancer_normalized=='breast cancer' & cancermineSentences$role=='Oncogene','year'])

paper.kras_driver_nsclc <- cancermineCollated[cancermineCollated$gene_normalized=='KRAS' & cancermineCollated$cancer_normalized=='non-small cell lung carcinoma' & cancermineCollated$role=='Driver','citation_count']
paper.kras_driver_nsclc_firstYear <- min(cancermineSentences[cancermineSentences$gene_normalized=='KRAS' & cancermineSentences$cancer_normalized=='non-small cell lung carcinoma' & cancermineSentences$role=='Driver','year'])

paper.triplesWithSingleCitations <- sum(cancermineCollated$citation_count==1)
paper.triplesWithMultipleCitations <- sum(cancermineCollated$citation_count>1)
paper.triplesWithAnyCitations <- nrow(cancermineCollated)
paper.percSingleCitations <- round(100*paper.triplesWithSingleCitations/paper.triplesWithAnyCitations,1)

paper.triplesWithSingleCitations <- prettyNum(paper.triplesWithSingleCitations,big.mark=',')
paper.triplesWithAnyCitations <- prettyNum(paper.triplesWithAnyCitations,big.mark=',')


fig_textminedOverview <- arrangeGrob(rolePlot,topCancersPlot,topGenesPlot,topJournalsPlot, citationPlot, ncol=2)

grid.arrange(fig_textminedOverview)
