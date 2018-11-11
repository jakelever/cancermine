## @knitr fig5

source('cancermine/dependencies.R')
source('cancermine/plotHeatmapWithDendro.R')

collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermineCollated <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


sentencesFilename <- 'cancermine/cancermine_sentences.tsv'
sentencesFilename <- normalizePath(sentencesFilename)
cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
cancermineSentences$journalShort <- strtrim(cancermineSentences$journal,21)


paper.clusteringTopCancerCount <- 30
cancerCounts <- plyr::count(cancermineSentences[,c('cancer_normalized'),drop=F])
cancerCounts <- cancerCounts[order(cancerCounts$freq,decreasing=T),]
topCancers <- cancerCounts[1:paper.clusteringTopCancerCount,'cancer_normalized']

#topCancers <- 
geneCounts <- plyr::count(cancermineSentences[,c('gene_normalized'),drop=F])
geneCounts <- geneCounts[order(geneCounts$freq,decreasing=T),]
topGenes <- geneCounts[1:paper.clusteringTopCancerCount,'gene_normalized']

topCancermineCounts <- cancermineCollated
topCancermineCounts <- topCancermineCounts[topCancermineCounts$cancer_normalized %in% topCancers,]
topCancermineCounts <- topCancermineCounts[topCancermineCounts$gene_normalized %in% topGenes,]
topCancermineCounts$gene_and_role <- factor(paste(topCancermineCounts$gene_normalized, topCancermineCounts$role, sep=':'))
topCancermineCounts$cancer_normalized <- factor(as.character(topCancermineCounts$cancer_normalized))

selectedCancers <- c('acute lymphocytic leukemia','glioblastoma multiforme')
selectedCancermineCounts <- cancermineCollated
selectedCancermineCounts <- selectedCancermineCounts[selectedCancermineCounts$cancer_normalized %in% selectedCancers,]
selectedCancermineCounts$gene_and_role <- factor(paste(selectedCancermineCounts$gene_normalized, selectedCancermineCounts$role, sep=':'))
selectedCancermineCounts$cancer_normalized <- factor(as.character(selectedCancermineCounts$cancer_normalized))


cancerCountsToHeatmapData <- function(myCancerCounts) {
  sm <- sparseMatrix(i=as.integer(myCancerCounts$cancer_normalized), 
                     j=as.integer(myCancerCounts$gene_and_role),
                     x=log10(myCancerCounts$citation_count),
                     #x=cancermineCounts$citation_count,
                     dims=c(nlevels(myCancerCounts$cancer_normalized),
                            nlevels(myCancerCounts$gene_and_role)))
  dense <- as.matrix(sm)
  colnames(dense) <- levels(myCancerCounts$gene_and_role)
  rownames(dense) <- levels(myCancerCounts$cancer_normalized)
  
  filteredDense <- dense
  filteredDense <- apply(filteredDense, 1, rescale)
  dim(filteredDense)
  rowMaxs <- apply(filteredDense,1,max)
  filteredDense <- filteredDense[rowMaxs>0.2,]
  dim(filteredDense)
  
  return(t(filteredDense))
}

heatmapCols <- brewer.pal(9,'YlOrRd')

topHeatmapData <- cancerCountsToHeatmapData(topCancermineCounts)
selectedHeatmapData <- cancerCountsToHeatmapData(selectedCancermineCounts)


fig_profiles <- plotHeatmapWithDendro(topHeatmapData)
fig_profiles <- arrangeGrob(fig_profiles,top='(a)')
grid.arrange(fig_profiles)

tcgaProfiles <- read.table('cancermine/cancermineProfilesTCGA.tsv',sep='\t',header=T,row.names=1)
tcgaProfiles.melted <- melt(as.matrix(tcgaProfiles))
colnames(tcgaProfiles.melted) <- c('profile','tcga','perc')

profileOrdering = c('breast cancer','colorectal cancer','hepatocellular carcinoma','lung cancer', 'malignant glioma', 'prostate cancer', 'stomach cancer', 'none')
tcgaOrdering = c('BRCA', 'COAD', 'LIHC', 'LUAD', 'LGG', 'PRAD', 'STAD')

cancerProfiles <- grep("none",unique(tcgaProfiles.melted$profile),value=T,invert=T)
#tcgaProfiles.melted$profile <- factor(tcgaProfiles.melted$profile,levels=c("none",rev(cancerProfiles)))
#tcgaProfiles.melted$tcga <- factor(tcgaProfiles.melted$tcga,levels=rev(unique(tcgaProfiles.melted$tcga)))

tcgaProfiles.melted$profile <- factor(tcgaProfiles.melted$profile,levels=rev(profileOrdering))
tcgaProfiles.melted$tcga <- factor(tcgaProfiles.melted$tcga,levels=tcgaOrdering)

## plot data
fig_tcga <- ggplot(tcgaProfiles.melted, aes(tcga,profile)) +
  theme_bw() +
  geom_tile(aes(fill = perc)) + 
  geom_text(aes(label = round(perc, 1))) +
  #scale_fill_gradient(low = "white", high = "red") +
  scale_fill_gradient("", high = "#6f6fce", low = "#e3f1f2") +
  #theme(axis.text.x = element_text(angle = 0, hjust = .5)) +
  labs(y = "CancerMine Profile", x = "TCGA Project", fill = "% top match") + 
  theme(axis.line=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
fig_tcga = arrangeGrob(fig_tcga,top='(b)')
grid.arrange(fig_tcga)

fig_profilesAndTCGA <- arrangeGrob(fig_profiles,ggplot(),fig_tcga,heights=c(2,0.1,1))

grid.arrange(fig_profilesAndTCGA)

paper.tcga_LGG_glioma_perc <- tcgaProfiles['malignant glioma','LGG']
