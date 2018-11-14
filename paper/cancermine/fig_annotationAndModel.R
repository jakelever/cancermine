## @knitr fig_annotationAndModel

source('cancermine/dependencies.R')

roleCounts <- read.table('cancermine/annotation.roleCounts.tsv', sep='\t', header=T)
rolePlot <- barchart(count ~ role, ylab="Annotations", xlab="Role", roleCounts, col="black")
rolePlot <- arrangeGrob(rolePlot, top='(a)')

interannotatorAgreement <- read.table('cancermine/annotation.interannotator.tsv', sep='\t', header=T, row.names=1)
interannotatorAgreement <- round(interannotatorAgreement,3)
grid.newpage()
interannotatorPlot <- tableGrob(interannotatorAgreement,
                                rows=c('Annotator 1','Annotator 2','Annotator 3'),
                                cols=c('Annotator 1','Annotator 2','Annotator 3'))
interannotatorPlot <- arrangeGrob(interannotatorPlot, top='(b)')
grid.arrange(interannotatorPlot)
#plot(tableGrob(interannotatorAgreement))

#ggplot(as.data.frame(table(df)), aes(x=gender, y = Freq, fill=fraud)) + geom_bar(stat="identity")



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
prCurvePlot <- arrangeGrob(prCurvePlot,top='(c)')

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
thresholdPlot <- arrangeGrob(thresholdPlot,top='(d)')

paper.lenient_driver_precision <- data[data$threshold==0.5 & data$role=='Driver','precision']
paper.lenient_driver_recall <- data[data$threshold==0.5 & data$role=='Driver','recall']
paper.lenient_oncogene_precision <- data[data$threshold==0.5 & data$role=='Oncogene','precision']
paper.lenient_oncogene_recall <- data[data$threshold==0.5 & data$role=='Oncogene','recall']
paper.lenient_ts_precision <- data[data$threshold==0.5 & data$role=='Tumor_Suppressor','precision']
paper.lenient_ts_recall <- data[data$threshold==0.5 & data$role=='Tumor_Suppressor','recall']

paper.lenient_avg_precision <- round(mean(c(paper.lenient_driver_precision,paper.lenient_oncogene_precision,paper.lenient_ts_precision)),1)
paper.lenient_avg_recall <- round(mean(c(paper.lenient_driver_recall,paper.lenient_oncogene_recall,paper.lenient_ts_recall)),1)

#data <- data[order(data$precision),]
#data[data$precision>0.85,]

fig_annotationAndModel <- arrangeGrob(rolePlot,interannotatorPlot,prCurvePlot,thresholdPlot, ncol=2)

grid.arrange(fig_annotationAndModel)

