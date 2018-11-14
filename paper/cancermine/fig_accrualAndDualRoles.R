
source('cancermine/dependencies.R')

collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermineCollated <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


sentencesFilename <- 'cancermine/cancermine_sentences.tsv'
sentencesFilename <- normalizePath(sentencesFilename)
cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


cancermineDM <- cancermineSentences[cancermineSentences$year!='None' & cancermineSentences$month!='None',]
cancermineDM$year <- as.integer(as.character(cancermineDM$year))
cancermineDM$month <- as.integer(as.character(cancermineDM$month))
cancermineDM$yearMonth <- cancermineDM$year + cancermineDM$month/12
cancermineDM <- cancermineDM[order(cancermineDM$yearMonth),]

head(cancermineDM,1)

novelSentences <- cancermineDM[!duplicated(cancermineDM$matching_id),c('matching_id','yearMonth')]

collatedWithFirstPub <- merge(cancermineCollated,novelSentences,by="matching_id",all.x = TRUE)

accrualPlot <- xyplot(citation_count ~ yearMonth, 
                      collatedWithFirstPub,
                      col="black",
                      xlab="First publication date",
                      ylab="Citations"
                      #scales=(x=list(log=10))
)
accrualPlot = arrangeGrob(accrualPlot,top='(a)')

grid.arrange(accrualPlot)


















myColours <- brewer.pal(3,"Dark2")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

#selectedTriples <- c('Oncogene BCL2 follicular lymphoma','Oncogene ERBB2 breast cancer','Oncogene MYC plasmacytoma')
#selectedTriples <- c('Tumor_Suppressor KLF6 prostate cancer','Oncogene ERBB2 breast cancer','Oncogene MYC plasmacytoma')

selectedTriples <- c('Oncogene NRAS breast cancer','Oncogene ERBB2 breast cancer','Tumor_Suppressor RUNX3 breast cancer')
cancermineSentences$combined <- paste(cancermineSentences$role,cancermineSentences$gene_normalized,cancermineSentences$cancer_normalized,sep=' ')
specialTripleCounts <- plyr::count(cancermineSentences[,c('combined','year')])
specialTripleCounts <- specialTripleCounts[specialTripleCounts$combined %in% selectedTriples,]
specialTripleCounts <- specialTripleCounts[specialTripleCounts$year<2018,]
selectedTriplesOverTimePlot <- xyplot(freq ~ year, 
                                      par.settings = my.settings,
                                      specialTripleCounts,
                                      type="l",
                                      groups=combined,
                                      xlab="Publication Year",
                                      ylab="Citations per year",
                                      auto.key=list(space="top", columns=1, 
                                                    points=FALSE, rectangles=TRUE),
                                      lwd=3)
selectedTriplesOverTimePlot = arrangeGrob(selectedTriplesOverTimePlot,top='(b)')
grid.arrange(selectedTriplesOverTimePlot)

foundAfter2000 <- cancermineSentences[cancermineSentences$year>=2000,'matching_id']
disappearedOnes <- cancermineCollated[!(cancermineCollated$matching_id %in% foundAfter2000),]
disappearedOnes <- disappearedOnes[order(disappearedOnes$citation_count,decreasing=T),]

forgottenRole <- head(disappearedOnes,1)
paper.forgottenRole_role <- tolower(as.character(forgottenRole$role))
paper.forgottenRole_gene <- as.character(forgottenRole$gene_normalized)
paper.forgottenRole_cancer <- as.character(forgottenRole$cancer_normalized)
paper.forgottenRole_before2010 <- forgottenRole$citation_count


paper.forgottenRole_count <- nrow(disappearedOnes)


cancermineSentences$decade = round(cancermineSentences$year,-1)
exploringTriples <- plyr::count(cancermineSentences[,c('role','gene_normalized','cancer_normalized','decade')])
exploringTriples <- exploringTriples[exploringTriples$cancer_normalized=='breast cancer',]
exploringTriples$roleGene <- paste(exploringTriples$role,exploringTriples$gene_normalized,sep='_')
#selected <- exploringTriples[exploringTriples$decade==1990 & exploringTriples$freq > 10,'roleGene']
#exploringTriples <- exploringTriples[exploringTriples$roleGene %in% selected,]
xyplot( log10(freq) ~ decade, exploringTriples, group=roleGene, type="l")

selected <- exploringTriples[exploringTriples$decade==1990 & exploringTriples$freq > 10,'roleGene']
exploringTriples2 <- exploringTriples[exploringTriples$roleGene %in% selected,]
#head(exploringTriples2)
exploringTriples2

tripTable <- dcast(exploringTriples, roleGene ~ decade, value.var="freq", fill=0)
colnames(tripTable) <- paste("x",colnames(tripTable),sep="")
tripTable[tripTable$x2010 > 20 & tripTable$x1990==0,]

paper.dualRole_citationsRequired <- 4
paper.dualRole_percRequired <- 90


#oncogeneCounts <- plyr::count(cancermine_byCitations[cancermine_byCitations$role=='Oncogene',c('gene_normalized','cancer_normalized')])
#tsCounts <- plyr::count(cancermine_byCitations[cancermine_byCitations$role=='Tumor_Suppressor',c('gene_normalized','cancer_normalized')])

oncogeneCounts <- cancermineCollated[cancermineCollated$role=='Oncogene',c('gene_normalized','cancer_normalized','citation_count')]
tsCounts <- cancermineCollated[cancermineCollated$role=='Tumor_Suppressor',c('gene_normalized','cancer_normalized','citation_count')]

oncogeneTSCounts <- merge(oncogeneCounts,tsCounts,by=c('gene_normalized','cancer_normalized'),all=T)
colnames(oncogeneTSCounts) <- c("gene_normalized","cancer_normalized",'Oncogene','Tumor_Suppressor')
oncogeneTSCounts[is.na(oncogeneTSCounts$Oncogene),'Oncogene'] <- 0
oncogeneTSCounts[is.na(oncogeneTSCounts$Tumor_Suppressor),'Tumor_Suppressor'] <- 0

oncogeneTSCounts$total <- oncogeneTSCounts$Oncogene + oncogeneTSCounts$Tumor_Suppressor
oncogeneTSCounts <- oncogeneTSCounts[oncogeneTSCounts$total >= paper.dualRole_citationsRequired,]
oncogeneTSCounts$Oncogene_perc <- oncogeneTSCounts$Oncogene / oncogeneTSCounts$total
oncogeneTSCounts$role <- 'unclear'
oncogeneTSCounts[oncogeneTSCounts$Oncogene_perc<(1-paper.dualRole_percRequired/100),'role'] <- 'Tumor_Suppressor'
oncogeneTSCounts[oncogeneTSCounts$Oncogene_perc>(paper.dualRole_percRequired/100),'role'] <- 'Oncogene'
#oncogeneTSCounts[oncogeneTSCounts$gene_normalized=='NOTCH1',]

oncogeneTSCounts$cancer_normalized_with_citations <- paste(oncogeneTSCounts$cancer_normalized, oncogeneTSCounts$total)

justOncogenes <- oncogeneTSCounts[oncogeneTSCounts$role=='Oncogene',]
justOncogenes$cancer_normalized_with_citations <- paste(justOncogenes$cancer_normalized, ' (',justOncogenes$Oncogene,')',sep='')

justTumorSuppressors <- oncogeneTSCounts[oncogeneTSCounts$role=='Tumor_Suppressor',]
justTumorSuppressors$cancer_normalized_with_citations <- paste(justTumorSuppressors$cancer_normalized, ' (',justTumorSuppressors$Tumor_Suppressor,')',sep='')

genesThatAreOncogenes <- plyr::count(oncogeneTSCounts[oncogeneTSCounts$role=='Oncogene','gene_normalized',drop=F])
genesThatAreTumorSuppressors <- plyr::count(oncogeneTSCounts[oncogeneTSCounts$role=='Tumor_Suppressor','gene_normalized',drop=F])
genesThatAreBoth <- merge(genesThatAreOncogenes,genesThatAreTumorSuppressors,by="gene_normalized")
colnames(genesThatAreBoth) <- c('gene','oncogenic role #', 'ts role #')
genesThatAreBoth <- genesThatAreBoth[grep("miR",genesThatAreBoth$gene,invert=T),]
genesThatAreBoth <- genesThatAreBoth[grep(";",genesThatAreBoth$gene,invert=T),]



genesThatAreBothWithCancerTypes <- data.frame(gene=genesThatAreBoth$gene)
rownames(genesThatAreBothWithCancerTypes) <- genesThatAreBothWithCancerTypes$gene
genesThatAreBothWithCancerTypes$oncogenic <- ''
genesThatAreBothWithCancerTypes$tumorSuppressive <- ''
for (gene in genesThatAreBothWithCancerTypes$gene) {
  genesThatAreBothWithCancerTypes[gene,'oncogenic'] <- paste(sort(as.character(justOncogenes[justOncogenes$gene_normalized==gene,'cancer_normalized_with_citations'])), sep=', ', collapse='\n')
  genesThatAreBothWithCancerTypes[gene,'tumorSuppressive'] <- paste(sort(as.character(justTumorSuppressors[justTumorSuppressors$gene_normalized==gene,'cancer_normalized_with_citations'])), sep=', ', collapse='\n')
}

oncogeneAndTSGenesPlot <- tableGrob(genesThatAreBothWithCancerTypes,rows=NULL, cols=c('Gene','Oncogenic in','Tumor Suppressive in'))
#oncogeneAndTSGenesPlot <- tableGrob(genesThatAreBothWithCancerTypes)
oncogeneAndTSGenesPlot <- arrangeGrob(oncogeneAndTSGenesPlot,top='(c)')
grid.arrange(oncogeneAndTSGenesPlot)


fig_accrualAndDualRoles <- arrangeGrob(accrualPlot,selectedTriplesOverTimePlot,oncogeneAndTSGenesPlot,ncol=1)

#fig_left <- arrangeGrob(accrualPlot,selectedTriplesOverTimePlot,ncol=1)
#fig_accrualAndDualRoles <- arrangeGrob(fig_left,oncogeneAndTSGenesPlot,ncol=1)
#grid.arrange(fig_all)

grid.arrange(fig_accrualAndDualRoles)
