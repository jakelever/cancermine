## @knitr fig5

source('cancermine/dependencies.R')
source('cancermine/plotHeatmapWithDendro.R')

collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermine <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")


selectDataForHeatmap <- function(cancermine,selectedCancers,roleCount) {
  myCancerCounts <- cancermine[cancermine$cancer_normalized %in% selectedCancers,]
  myCancerCounts$gene_and_role <- factor(paste(myCancerCounts$gene_normalized, myCancerCounts$role, sep=':'))
  
  selectedRoles <- myCancerCounts[,c('gene_and_role','citation_count')]
  selectedRoles <- selectedRoles[order(selectedRoles$citation_count,decreasing=T),]
  selectedRoles <- selectedRoles[!duplicated(selectedRoles$gene_and_role),]
  selectedRoles <- selectedRoles[1:roleCount,]
  
  myCancerCounts <- cancermine[cancermine$cancer_normalized %in% selectedCancers,]
  myCancerCounts$gene_and_role <- factor(paste(myCancerCounts$gene_normalized, myCancerCounts$role, sep=':'))
  myCancerCounts <- myCancerCounts[myCancerCounts$gene_and_role %in% selectedRoles$gene_and_role,]
  myCancerCounts$cancer_normalized <- factor(myCancerCounts$cancer_normalized)
  myCancerCounts$gene_and_role <- factor(myCancerCounts$gene_and_role)
  
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


paper.clusteringTopCancerCount <- 30
cancerCounts <- aggregate(cancermine$citation_count, by=list(cancer_normalized=cancermine$cancer_normalized), FUN=sum)
colnames(cancerCounts) <- c('cancer_normalized','total_citation_count')
cancerCounts <- cancerCounts[order(cancerCounts$total_citation_count,decreasing=T),]
topCancers <- cancerCounts[1:paper.clusteringTopCancerCount,'cancer_normalized']

topHeatmapData <- selectDataForHeatmap(cancermine,topCancers,40)

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
