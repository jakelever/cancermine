## @knitr fig5

source('cancermine/dependencies.R')

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
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

heatmapCols <- brewer.pal(9,'YlOrRd')
#library(gridGraphics)
#library(grid)

#melt_topHeatmapData <- melt(topHeatmapData)
#colnames(melt_topHeatmapData) <- c('cancer_normalized','gene_normalized','value')


topHeatmapData <- cancerCountsToHeatmapData(topCancermineCounts)
selectedHeatmapData <- cancerCountsToHeatmapData(selectedCancermineCounts)

plotHeatmap <- function() {
  #pdf(file="sideeffects2.pdf")
  heatmap.2(topHeatmapData, trace='none',dendrogram='row', col=heatmapCols, margins=c(15,15))
  #topHeatmapPlot <- grab_grob()
  #topHeatmapPlot <- arrangeGrob(topHeatmapPlot,top="(a)")
  
  #heatmap.2(selectedHeatmapData, trace='none',dendrogram='row', col=heatmapCols, margins=c(15,15), cexRow=1)
  #selectedHeatmapPlot <- grab_grob()
  #selectedHeatmapPlot <- arrangeGrob(selectedHeatmapPlot,top="(b)")
  
  #fig_clusters <- arrangeGrob(topHeatmapPlot,selectedHeatmapPlot,ncol=1)
  #dev.off()
  
  #grid.arrange(fig_clusters)
}

plotHeatmap()

mat <- topHeatmapData
sample_names <- colnames(mat)

# Obtain the dendrogram
dend <- as.dendrogram(hclust(dist(mat)))
dend_data <- dendro_data(dend)

#sample_cluster <- hclust(dist(t(mat)))
#sample_cluster$labels[c(sample_cluster$order)]
tmpHeatmap <- heatmap.2(topHeatmapData, trace='none',dendrogram='row', col=heatmapCols, margins=c(15,15))


segment_data <- with(
  segment(dend_data), 
  data.frame(x = y, y = x, xend = yend, yend = xend))

gene_pos_table <- with(
  dend_data$labels, 
  data.frame(y_center = x, gene = as.character(label), height = 1))

sample_pos_table <- data.frame(sample = sample_names[tmpHeatmap$colInd], x_center = 1:length(sample_names), width = 1)
#sample_pos_table <- data.frame(sample = sample_names[tmpHeatmap$colInd], x_center = 1:, width = 1)

heatmap_data <- mat %>% 
  reshape2::melt(value.name = "expr", varnames = c("gene", "sample")) %>%
  left_join(gene_pos_table, by="gene") %>%
  left_join(sample_pos_table, by="sample")

gene_axis_limits <- with(
  gene_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1


plt_hmap <- ggplot(heatmap_data, 
                   aes(x = x_center, y = y_center, fill = expr, 
                       height = height, width = width)) + 
  geom_tile() +
  scale_fill_gradient("", high = "#b50000", low = "#fffccc") +
  scale_x_continuous(breaks = sample_pos_table$x_center, 
                     labels = sample_pos_table$sample, 
                     expand = c(0, 0)) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = gene_pos_table[, "y_center"], 
                     #labels = rep("", nrow(gene_pos_table)),
                     labels = gene_pos_table$gene,
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Cancer Gene Role", y = "Cancer Type") +
  theme_bw() +
  theme(axis.text.x = element_text(size = rel(1), hjust = 1, vjust = 0.5, angle = 90), 
        # margin: top, right, bottom, and left
        plot.margin = unit(c(.2, -.5, 0.2, 0.2), "cm"), 
        axis.line=element_blank(),
        axis.ticks=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank()) +
  theme(legend.position="left")

# Dendrogram plot
plt_dendr <- ggplot(segment_data) + 
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
  scale_x_continuous(expand = c(0, 0.5), 
                     labels=NULL) + 
  scale_y_continuous(breaks = gene_pos_table$y_center, 
                     labels = NULL,#gene_pos_table$gene, 
                     limits = gene_axis_limits, 
                     expand = c(0, 0)) + 
  labs(x = "Distance", y = "", colour = "", size = "") +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank(),
        plot.margin = unit(c(.2, 0.2, 0.2, 0), "cm"))


fig_profiles <- plot_grid(plt_hmap, plt_dendr, align = 'h', rel_widths = c(6, 1))
fig_profiles = arrangeGrob(fig_profiles,top='(a)')
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
