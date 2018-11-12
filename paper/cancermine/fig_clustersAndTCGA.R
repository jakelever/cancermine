## @knitr fig5

source('cancermine/dependencies.R')
source('cancermine/plotHeatmapWithDendro.R')

collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermine <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

generateProfiles <- function(cancermine) {
  cancermineProfiles <- cancermine
  cancermineProfiles$gene_and_role <- factor(paste(cancermineProfiles$gene_normalized, str_replace(cancermineProfiles$role,"_"," "), sep=':'))
  cancermineProfiles$log_citation_count <- log10(cancermineProfiles$citation_count+1)
  cancermineProfiles <- cancermineProfiles[,c('cancer_normalized','gene_and_role','log_citation_count')]
  
  maxCitationsByCancer <- aggregate(cancermineProfiles$log_citation_count,by=list(cancer_normalized=cancermineProfiles$cancer_normalized),FUN=max)
  colnames(maxCitationsByCancer) <- c('cancer_normalized','max_log_citations_for_cancer')
  
  cancermineProfiles <- inner_join(cancermineProfiles,maxCitationsByCancer,by="cancer_normalized")
  cancermineProfiles$importance_score <- cancermineProfiles$log_citation_count / cancermineProfiles$max_log_citations_for_cancer
  cancermineProfiles <- cancermineProfiles[,c('cancer_normalized','gene_and_role','importance_score')]
  
  return(cancermineProfiles)
}

selectDataForHeatmap <- function(cancermineProfiles,selectedCancers,roleCount) {
  
  selectedGeneRoles <- cancermineProfiles[cancermineProfiles$cancer_normalized %in% selectedCancers,c('gene_and_role','importance_score')]
  selectedGeneRoles <- selectedGeneRoles[order(selectedGeneRoles$importance_score,decreasing=T),]
  selectedGeneRoles <- selectedGeneRoles[!duplicated(selectedGeneRoles$gene_and_role),'gene_and_role']
  selectedGeneRoles <- selectedGeneRoles[1:roleCount]
  
  selectedImportanceScores <- cancermineProfiles[cancermineProfiles$cancer_normalized %in% selectedCancers,]
  selectedImportanceScores <- selectedImportanceScores[selectedImportanceScores$gene_and_role %in% selectedGeneRoles,]
  
  dense <- acast(selectedImportanceScores, 
                 cancer_normalized ~ gene_and_role,
                 value.var="importance_score",
                 fill=0)
  
  return(dense)
}

cancermineProfiles <- generateProfiles(cancermine)

paper.clusteringTopCancerCount <- 30
cancerCounts <- aggregate(cancermine$citation_count, by=list(cancer_normalized=cancermine$cancer_normalized), FUN=sum)
colnames(cancerCounts) <- c('cancer_normalized','total_citation_count')
cancerCounts <- cancerCounts[order(cancerCounts$total_citation_count,decreasing=T),]
topCancers <- cancerCounts[1:paper.clusteringTopCancerCount,'cancer_normalized']

topHeatmapData <- selectDataForHeatmap(cancermineProfiles,topCancers,40)

moreCancers <- cancerCounts[1:50,'cancer_normalized']
allData <- selectDataForHeatmap(cancermineProfiles,moreCancers,10000)
profilePCA <- prcomp(allData)
twoPCs <- as.data.frame(profilePCA$x[,1:2])
twoPCs$cancer_normalized <- row.names(twoPCs)

ggplot(twoPCs[1:50,], aes(x= PC1, y = PC2)) + 
  geom_point(color = "blue", size = 3) + 
  geom_label_repel(aes(label = cancer_normalized),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50')

distMatrix <- as.matrix(dist(allData))
heatmap.2(distMatrix,trace="none")

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
