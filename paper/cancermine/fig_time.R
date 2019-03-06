
source('cancermine/dependencies.R')

sentencesFilename <- 'cancermine/cancermine_sentences.tsv'
sentencesFilename <- normalizePath(sentencesFilename)
cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")

cancermineSentences <- cancermineSentences[order(cancermineSentences$year),]
cancermineSentences$isNovel <- !duplicated(cancermineSentences[,c('role','gene_hugo_id','cancer_id')])



novelYearCounts <- plyr::count(cancermineSentences[cancermineSentences$isNovel,c('year'),drop=F])
novelYearCounts$novel <- 'Novel'
notnovelYearCounts <- plyr::count(cancermineSentences[!cancermineSentences$isNovel,c('year'),drop=F])
notnovelYearCounts$novel <- 'Not Novel'

novelAndNotNovelYearCounts <- rbind(novelYearCounts,notnovelYearCounts)
novelAndNotNovelYearCounts$yearTxt <- as.character(novelAndNotNovelYearCounts$year)
novelAndNotNovelYearCounts <- novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year<2019,]

novelNotNovel2018 <- novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2018,]
novel2018 <- novelNotNovel2018[novelNotNovel2018$novel=='Novel','freq']
notnovel2018 <- novelNotNovel2018[novelNotNovel2018$novel=='Not Novel','freq']

paper.notnovel2018Perc = round(100*notnovel2018/(novel2018+notnovel2018))
paper.novel2018Perc = round(100*novel2018/(novel2018+notnovel2018))

paper.rolesIn2018 <- sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2018,'freq'])
paper.avgRolesPerMonth <- round(sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2018,'freq'])/12)
paper.avgNovelRolesPerMonth <- round(sum(novel2018/12))
paper.rolesIn2018 <- prettyNum(paper.rolesIn2018,big.mark=',')

myColours <- brewer.pal(7,"Spectral")
myColours <- c(myColours[1],myColours[length(myColours)])
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

labels <- sort(unique(novelAndNotNovelYearCounts$year))
labels[labels%%5!=0] <- ''
ratePlot <- barchart(freq ~ yearTxt, 
                     novelAndNotNovelYearCounts, 
                     groups=novel,
                     ylab="Mentions per year",
                     par.settings = my.settings,
                     auto.key=list(space="top", columns=2, 
                                   points=FALSE, rectangles=TRUE),
                     scales=list(x=list(rot=45,labels=labels)),
                     stack=T,
                     horizontal=F)
ratePlot = arrangeGrob(ratePlot,top='(a)')
grid.arrange(ratePlot)






title <- plyr::count(cancermineSentences[cancermineSentences$section=='title',c('year'),drop=F])
title$section <- 'Title'
abstract <- plyr::count(cancermineSentences[cancermineSentences$section=='abstract',c('year'),drop=F])
abstract$section <- 'Abstract'
article <- plyr::count(cancermineSentences[cancermineSentences$section=='article',c('year'),drop=F])
article$section <- 'Article'

sectionData <- rbind(title,abstract,article)
sectionData$yearTxt <- as.character(sectionData$year)
sectionData <- sectionData[sectionData$year<2019,]

myColours <- brewer.pal(3,"Spectral")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

labels <- sort(unique(sectionData$year))
labels[labels%%5!=0] <- ''
sectionPlot <- barchart(freq ~ yearTxt, 
                        sectionData, 
                     groups=section,
                     ylab="Mentions per year",
                     par.settings = my.settings,
                     auto.key=list(space="top", columns=1, 
                                   points=FALSE, rectangles=TRUE),
                     scales=list(x=list(rot=45,labels=labels)),
                     stack=T,
                     horizontal=F)
sectionPlotWithoutHeader <- sectionPlot
sectionPlot = arrangeGrob(sectionPlot,top='(c)')
grid.arrange(sectionPlot)

paper.novelGeneRolesPerMonth2018 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$year==2018,]) / 12)
paper.novelDriverPerMonth2018 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Driver' & cancermineSentences$year==2018,]) / 12)
paper.novelOncogenePerMonth2018 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Oncogene' & cancermineSentences$year==2018,]) / 12)
paper.novelTumorSuppressorPerMonth2018 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Tumor_Suppressor' & cancermineSentences$year==2018,]) / 12)

driver <- plyr::count(cancermineSentences[cancermineSentences$role=='Driver',c('year'),drop=F])
driver$role <- 'Driver'
oncogene <- plyr::count(cancermineSentences[cancermineSentences$role=='Oncogene',c('year'),drop=F])
oncogene$role <- 'Oncogene'
tumor_suppressor <- plyr::count(cancermineSentences[cancermineSentences$role=='Tumor_Suppressor',c('year'),drop=F])
tumor_suppressor$role <- 'Tumor Suppressor'

roleData <- rbind(driver,oncogene,tumor_suppressor)
roleData$yearTxt <- as.character(roleData$year)
roleData <- roleData[roleData$year<2019,]

paper.driver2018 <- roleData[roleData$year==2018 & roleData$role=='Driver','freq']
paper.oncogene2018 <- roleData[roleData$year==2018 & roleData$role=='Oncogene','freq']
paper.tumorsuppressor2018 <- roleData[roleData$year==2018 & roleData$role=='Tumor Suppressor','freq']


#paper.driverPerMonth <- round(paper.driver2018/12)
#paper.oncogenePerMonth <- round(paper.oncogene2018/12)
#paper.tumorsuppressorPerMonth <- round(paper.tumorsuppressor2018/12)

paper.driver2018 <- prettyNum(paper.driver2018,big.mark=',')
paper.oncogene2018 <- prettyNum(paper.oncogene2018,big.mark=',')
paper.tumorsuppressor2018 <- prettyNum(paper.tumorsuppressor2018,big.mark=',')

myColours <- brewer.pal(3,"Spectral")
my.settings <- list(
  superpose.polygon=list(col=myColours),
  #strip.background=list(col=myColours),
  superpose.line=list(col=myColours),
  strip.border=list(col="black")
)

labels <- sort(unique(roleData$year))
labels[labels%%5!=0] <- ''
rolePlot <- barchart(freq ~ yearTxt, 
                     roleData, 
                     groups=role,
                     ylab="Mentions per year",
                     par.settings = my.settings,
                     auto.key=list(space="top", columns=1, 
                                   points=FALSE, rectangles=TRUE),
                     scales=list(x=list(rot=45,labels=labels)),
                     stack=T,
                     horizontal=F)
rolePlot = arrangeGrob(rolePlot,top='(b)')
grid.arrange(rolePlot)

unique(as.character(cancermineSentences$subsection))
cancermineSentences$subsection <- as.character(cancermineSentences$subsection)
cancermineSentences[cancermineSentences$subsection=='conclusions','subsection'] <- 'conclusion'
cancermineSentences[cancermineSentences$subsection=='materials and methods','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='background','subsection'] <- 'introduction'
cancermineSentences[cancermineSentences$subsection=='materials and methods','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='study design','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='ethics statement','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='statistical analysis','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='data analysis','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='abbreviations','subsection'] <- 'introduction'
cancermineSentences[cancermineSentences$subsection=='summary','subsection'] <- 'introduction'
cancermineSentences[cancermineSentences$subsection=='supporting information','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='material and methods','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='method','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='statistical analyses','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='additional information','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='materials','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='conflicts of interest','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='data collection','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='author contributions','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='case report','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='patients and method','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='patients and methods','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='limitations','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='results and discussion results','subsection'] <- 'results'
cancermineSentences[cancermineSentences$subsection=='results and discussion','subsection'] <- 'results'
cancermineSentences[cancermineSentences$subsection=='supplementary material','subsection'] <- 'methods'
cancermineSentences[cancermineSentences$subsection=='None','subsection'] <- 'Unable to identify'
#cancermineSentences[cancermineSentences$subsection!='abbreviations',]

cancermineSentences$subsection <- factor(as.character(cancermineSentences$subsection),levels=c('introduction','methods','results','discussion','conclusion','Unable to identify'))

subsectionCounts <- plyr::count(cancermineSentences[cancermineSentences$section=='article',c('subsection'),drop=F])
subsectionPlot <- barchart(freq ~ subsection, 
                           subsectionCounts,
                           ylab="Mentions",
                           scales=list(x=list(rot=45)),
                           col="black")
subsectionPlotWithoutHeader <- subsectionPlot
subsectionPlot <- arrangeGrob(subsectionPlot,top='(d)')


subsectionCountsWithNovelty <- plyr::count(cancermineSentences[cancermineSentences$section=='article',c('subsection','isNovel'),drop=F])
subsectionCountsWithNovelty$isNovelTxt <- revalue(as.character(subsectionCountsWithNovelty$isNovel),c('FALSE'='Not Novel','TRUE'='Novel'))
subsectionPlotWithNovelty <- barchart(freq ~ subsection, 
         subsectionCountsWithNovelty,
         groups=isNovelTxt,
         auto.key=list(space="top", columns=2, 
                       points=FALSE, rectangles=TRUE),
         xlab="Section of paper",
         ylab="Mentions",
         ylim=c(0,1.1*max(subsectionCountsWithNovelty$freq)),
         scales=list(x=list(rot=45)))








genesAndYears <- cancermineSentences[,c('gene_entrez_id','year')]
genesAndYears <- genesAndYears[!duplicated(genesAndYears),]

uniqueGenesByYear <- plyr::count(genesAndYears[,c('year'),drop=F])
uniqueGenesByYear$yearTxt <- as.character(uniqueGenesByYear$year)
uniqueGenesByYear <- uniqueGenesByYear[uniqueGenesByYear$year<2019,]

uniqueGenesPlot <- barchart(freq ~ yearTxt, 
                     uniqueGenesByYear,
                     col='black',
                     ylab="Mentions per year",
                     scales=list(x=list(rot=45,labels=labels)),
                     horizontal=F)
uniqueGenesPlot = arrangeGrob(uniqueGenesPlot,top='(e)')
grid.arrange(uniqueGenesPlot)


cancersAndYears <- cancermineSentences[,c('cancer_id','year')]
cancersAndYears <- cancersAndYears[!duplicated(cancersAndYears),]

uniqueCancersByYear <- plyr::count(cancersAndYears[,c('year'),drop=F])
uniqueCancersByYear$yearTxt <- as.character(uniqueCancersByYear$year)
uniqueCancersByYear <- uniqueCancersByYear[uniqueCancersByYear$year<2019,]

uniqueCancersPlot <- barchart(freq ~ yearTxt, 
                            uniqueCancersByYear,
                            col='black',
                            ylab="Mentions per year",
                            scales=list(x=list(rot=45,labels=labels)),
                            horizontal=F)
uniqueCancersPlot = arrangeGrob(uniqueCancersPlot,top='(f)')
grid.arrange(uniqueCancersPlot)




fig_time <- arrangeGrob(ratePlot,rolePlot,sectionPlot,subsectionPlot,uniqueGenesPlot,uniqueCancersPlot,ncol=2)

grid.arrange(fig_time)


sentencesWithNovelty <- cancermineSentences[cancermineSentences$section=='article',]
sentencesWithNovelty <- sentencesWithNovelty[order(sentencesWithNovelty$year,sentencesWithNovelty$month,sentencesWithNovelty$day),]
sentencesWithNovelty$isNovel <- !duplicated(sentencesWithNovelty[,c('role','gene_entrez_id','cancer_id')])

noveltyWithSubsectionTable <- table(sentencesWithNovelty$isNovel,sentencesWithNovelty$subsection)
noveltyWithSubsectionTest <- chisq.test(noveltyWithSubsectionTable)

paper.noveltyWithSubsection_pvalue <- signif(noveltyWithSubsectionTest$p.value,2)
