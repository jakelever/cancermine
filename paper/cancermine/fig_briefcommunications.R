
source('cancermine/dependencies.R')

Driver <- read.table('cancermine/prCurves/Driver.txt',header=T)
Driver$role <- 'Driver'

Oncogene <- read.table('cancermine/prCurves/Oncogene.txt',header=T)
Oncogene$role <- 'Oncogene'

Tumor_Suppressor <- read.table('cancermine/prCurves/Tumor_Suppressor.txt',header=T)
Tumor_Suppressor$role <- 'Tumor_Suppressor'

data <- rbind(Driver,Oncogene,Tumor_Suppressor)
data <- data[order(data$precision,decreasing=T),]
data <- data[order(data$recall),]

prCurvePlot <- xyplot(precision ~ recall | role, 
                      xlab="Recall", ylab="Precision",
                      #xlim=c(0,1),ylim=c(0,1),
                      data, 
                      lwd=2,
                      type="l")

myColours <- c(brewer.pal(3,"Dark2"),"#000000")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)


selectedThresholds <- data.frame(role=c('Driver','Oncogene','Tumor_Suppressor'),
                                 threshold=c(0.80,0.76,0.92))

data <- data[order(data$threshold),]
write.table(data,'fig1a_data.tsv',sep='\t',quote=FALSE,row.names=FALSE)

thresholdPlot <- xyplot(precision + recall ~ threshold | role, 
                        xlab="Threshold", ylab="Precision / Recall",
                        #auto.key=T,
                        par.settings = my.settings,
                        auto.key=list(space="top", columns=2, 
                                      points=FALSE, rectangles=TRUE),
                        lwd=2,
                        panel = function(x, ...) {
                          role <- dimnames(trellis.last.object())[["role"]][packet.number()]
                          panel.xyplot(x, ...)
                          panel.abline(v=selectedThresholds$threshold[selectedThresholds$role==role])
                        }, 
                        data, type="l")

prAndThresholdPlot <- arrangeGrob(prCurvePlot,thresholdPlot,top='(a)')






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

write.table(geneCounts,'fig1b_genes.tsv',sep='\t',quote=FALSE,row.names=FALSE)

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


cancerCounts <- plyr::count(cancermineCollated[,c('cancer_normalized'),drop=F])
cancerCounts <- cancerCounts[order(cancerCounts$freq,decreasing=T),]
cancerCounts$cancer_normalized <- factor(as.character(cancerCounts$cancer_normalized), levels=as.character(cancerCounts$cancer_normalized))

cancerRoleCounts <- plyr::count(cancermineCollated[,c('cancer_normalized','role'),drop=F])
cancerRoleCounts <- suppressWarnings(inner_join(cancerRoleCounts,cancerCounts,by="cancer_normalized"))
colnames(cancerRoleCounts) <- c('cancer_normalized','role','role_citation_count','cancer_citation_count')
cancerRoleCounts <- cancerRoleCounts[order(cancerRoleCounts$cancer_citation_count,decreasing=T),]
cancerRoleCounts$cancer_normalized <- factor(as.character(cancerRoleCounts$cancer_normalized), levels=unique(as.character(cancerRoleCounts$cancer_normalized)))
topCancers <- cancerCounts[1:topCount,'cancer_normalized']

write.table(cancerCounts,'fig1b_cancers.tsv',sep='\t',quote=FALSE,row.names=FALSE)


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

dataOverviewPlot <- arrangeGrob(topGeneRolesPlot,topCancerRolesPlot,top='(b)')

fig_briefcommunications1 <- arrangeGrob(prAndThresholdPlot,dataOverviewPlot,ncol=2)
grid.arrange(fig_briefcommunications1)



fig_comparisonsWithHeader <- arrangeGrob(fig_comparisons,top='(a)')

img <- readPNG("clustering.png")
fig_profiles <- rasterGrob(img, interpolate=TRUE)
fig_profiles <- arrangeGrob(fig_profiles,top='(b)')
grid.arrange(fig_profiles)

fig_briefcommunications2 <- arrangeGrob(
  fig_comparisonsWithHeader,
  fig_profiles,
  nrow=2, 
  heights=c(0.4,0.6))
grid.arrange(fig_briefcommunications2)
