
source('cancermine/dependencies.R')

myColours <- c(brewer.pal(4,"Set2"),"#000000")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)


collatedFilename <- 'cancermine/cancermine_collated.tsv'
collatedFilename <- normalizePath(collatedFilename)
cancermineCollated <- read.table(collatedFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

cancermineCollated$role <- factor(str_replace(as.character(cancermineCollated$role),"_"," "))

topCount <- 10

geneCounts <- plyr::count(cancermineCollated[,c('gene_normalized'),drop=F])
geneCounts <- geneCounts[order(geneCounts$freq,decreasing=T),]
geneCounts$gene_normalized <- factor(as.character(geneCounts$gene_normalized), levels=as.character(geneCounts$gene_normalized))

geneRoleCounts <- plyr::count(cancermineCollated[,c('gene_normalized','role'),drop=F])
geneRoleCounts <- suppressWarnings(inner_join(geneRoleCounts,geneCounts,by="gene_normalized"))
colnames(geneRoleCounts) <- c('gene_normalized','role','role_citation_count','gene_citation_count')
geneRoleCounts <- geneRoleCounts[order(geneRoleCounts$gene_citation_count,decreasing=T),]
geneRoleCounts$gene_normalized <- factor(as.character(geneRoleCounts$gene_normalized), levels=unique(as.character(geneRoleCounts$gene_normalized)))
topGenes <- geneCounts[1:topCount,'gene_normalized']

topGeneRolesPlot <- barchart(role_citation_count ~ gene_normalized, 
         ylab="Associations",
         ylim=c(0,1.1*max(geneCounts$freq)),
         geneRoleCounts[geneRoleCounts$gene_normalized %in% topGenes,], 
         scales=list(x=list(rot=45)),
         groups=role,
         par.settings = my.settings,
         auto.key=F,
         stack=T)
topGeneRolesPlot <- arrangeGrob(topGeneRolesPlot,top='(a)')


cancerCounts <- plyr::count(cancermineCollated[,c('cancer_normalized'),drop=F])
cancerCounts <- cancerCounts[order(cancerCounts$freq,decreasing=T),]
cancerCounts$cancer_normalized <- factor(as.character(cancerCounts$cancer_normalized), levels=as.character(cancerCounts$cancer_normalized))

cancerRoleCounts <- plyr::count(cancermineCollated[,c('cancer_normalized','role'),drop=F])
cancerRoleCounts <- suppressWarnings(inner_join(cancerRoleCounts,cancerCounts,by="cancer_normalized"))
colnames(cancerRoleCounts) <- c('cancer_normalized','role','role_citation_count','cancer_citation_count')
cancerRoleCounts <- cancerRoleCounts[order(cancerRoleCounts$cancer_citation_count,decreasing=T),]
cancerRoleCounts$cancer_normalized <- factor(as.character(cancerRoleCounts$cancer_normalized), levels=unique(as.character(cancerRoleCounts$cancer_normalized)))
topCancers <- cancerCounts[1:topCount,'cancer_normalized']

topCancerRolesPlot <- barchart(role_citation_count ~ cancer_normalized, 
                             ylab="Associations",
                             ylim=c(0,1.1*max(cancerCounts$freq)),
                             cancerRoleCounts[cancerRoleCounts$cancer_normalized %in% topCancers,], 
                             scales=list(x=list(rot=45)),
                             groups=role,
                             par.settings = my.settings,
                             auto.key=list(space="bottom", columns=1, 
                                           points=FALSE, rectangles=TRUE),
                             stack=T)
topCancerRolesPlot <- arrangeGrob(topCancerRolesPlot,top='(b)')


img <- readPNG("clustering.png")
fig_profiles <- rasterGrob(img, interpolate=TRUE)
fig_profiles <- arrangeGrob(fig_profiles,top='(c)')
grid.arrange(fig_profiles)

fig_correspondence <- arrangeGrob(
  arrangeGrob(
    topGeneRolesPlot,
    topCancerRolesPlot,
    ncol=1,
    heights=c(0.4,0.6)
    ),
  fig_profiles,
  nrow=1, 
  widths=c(0.33,0.67))
grid.arrange(fig_correspondence)
