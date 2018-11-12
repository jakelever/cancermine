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



plotHeatmapWithDendro <- function(mat) {
  sample_names <- colnames(mat)
  
  # Obtain the dendrogram
  dend <- as.dendrogram(hclust(dist(mat)))
  dend_data <- dendro_data(dend)
  
  heatmapCols <- brewer.pal(9,'YlOrRd')
  tmpHeatmap <- heatmap.2(mat, trace='none',dendrogram='row', col=heatmapCols, margins=c(15,15))
  
  segment_data <- with(
    segment(dend_data), 
    data.frame(x = y, y = x, xend = yend, yend = xend))
  
  gene_pos_table <- with(
    dend_data$labels, 
    data.frame(y_center = x, gene = as.character(label), height = 1))
  
  sample_pos_table <- data.frame(sample = sample_names[tmpHeatmap$colInd], x_center = 1:length(sample_names), width = 1)
  
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
  
  
  heatmapAndDendro <- plot_grid(plt_hmap, plt_dendr, align = 'h', rel_widths = c(6, 1))
  return(heatmapAndDendro)
}

