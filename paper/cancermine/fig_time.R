
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
novelAndNotNovelYearCounts <- novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year<2018,]

novelNotNovel2017 <- novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2017,]
novel2017 <- novelNotNovel2017[novelNotNovel2017$novel=='Novel','freq']
notnovel2017 <- novelNotNovel2017[novelNotNovel2017$novel=='Not Novel','freq']

paper.notnovel2017Perc = round(100*notnovel2017/(novel2017+notnovel2017))
paper.novel2017Perc = round(100*novel2017/(novel2017+notnovel2017))

paper.rolesIn2017 <- sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2017,'freq'])
paper.avgRolesPerMonth <- round(sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2017,'freq'])/12)
paper.avgNovelRolesPerMonth <- round(sum(novel2017/12))
paper.rolesIn2017 <- prettyNum(paper.rolesIn2017,big.mark=',')

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
sectionData <- sectionData[sectionData$year<2018,]

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
                     par.settings = my.settings,
                     auto.key=list(space="top", columns=1, 
                                   points=FALSE, rectangles=TRUE),
                     scales=list(x=list(rot=45,labels=labels)),
                     stack=T,
                     horizontal=F)
sectionPlot = arrangeGrob(sectionPlot,top='(c)')
grid.arrange(sectionPlot)


paper.novelDriverPerMonth2017 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Driver' & cancermineSentences$year==2017,]) / 12)
paper.novelOncogenePerMonth2017 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Oncogene' & cancermineSentences$year==2017,]) / 12)
paper.novelTumorSuppressorPerMonth2017 <- round(nrow(cancermineSentences[cancermineSentences$isNovel & cancermineSentences$role=='Tumor_Suppressor' & cancermineSentences$year==2017,]) / 12)

driver <- plyr::count(cancermineSentences[cancermineSentences$role=='Driver',c('year'),drop=F])
driver$role <- 'Driver'
oncogene <- plyr::count(cancermineSentences[cancermineSentences$role=='Oncogene',c('year'),drop=F])
oncogene$role <- 'Oncogene'
tumor_suppressor <- plyr::count(cancermineSentences[cancermineSentences$role=='Tumor_Suppressor',c('year'),drop=F])
tumor_suppressor$role <- 'Tumor Suppressor'

roleData <- rbind(driver,oncogene,tumor_suppressor)
roleData$yearTxt <- as.character(roleData$year)
roleData <- roleData[roleData$year<2018,]

paper.driver2017 <- roleData[roleData$year==2017 & roleData$role=='Driver','freq']
paper.oncogene2017 <- roleData[roleData$year==2017 & roleData$role=='Oncogene','freq']
paper.tumorsuppressor2017 <- roleData[roleData$year==2017 & roleData$role=='Tumor Suppressor','freq']


#paper.driverPerMonth <- round(paper.driver2017/12)
#paper.oncogenePerMonth <- round(paper.oncogene2017/12)
#paper.tumorsuppressorPerMonth <- round(paper.tumorsuppressor2017/12)

paper.driver2017 <- prettyNum(paper.driver2017,big.mark=',')
paper.oncogene2017 <- prettyNum(paper.oncogene2017,big.mark=',')
paper.tumorsuppressor2017 <- prettyNum(paper.tumorsuppressor2017,big.mark=',')

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
                           scales=list(x=list(rot=45)),
                           col="black")
subsectionPlot <- arrangeGrob(subsectionPlot,top='(d)')








genesAndYears <- cancermineSentences[,c('gene_entrez_id','year')]
genesAndYears <- genesAndYears[!duplicated(genesAndYears),]

uniqueGenesByYear <- plyr::count(genesAndYears[,c('year'),drop=F])
uniqueGenesByYear$yearTxt <- as.character(uniqueGenesByYear$year)
uniqueGenesByYear <- uniqueGenesByYear[uniqueGenesByYear$year<2018,]

uniqueGenesPlot <- barchart(freq ~ yearTxt, 
                     uniqueGenesByYear,
                     col='black',
                     scales=list(x=list(rot=45,labels=labels)),
                     horizontal=F)
uniqueGenesPlot = arrangeGrob(uniqueGenesPlot,top='(e)')
grid.arrange(uniqueGenesPlot)


cancersAndYears <- cancermineSentences[,c('cancer_id','year')]
cancersAndYears <- cancersAndYears[!duplicated(cancersAndYears),]

uniqueCancersByYear <- plyr::count(cancersAndYears[,c('year'),drop=F])
uniqueCancersByYear$yearTxt <- as.character(uniqueCancersByYear$year)
uniqueCancersByYear <- uniqueCancersByYear[uniqueCancersByYear$year<2018,]

uniqueCancersPlot <- barchart(freq ~ yearTxt, 
                            uniqueCancersByYear,
                            col='black',
                            scales=list(x=list(rot=45,labels=labels)),
                            horizontal=F)
uniqueCancersPlot = arrangeGrob(uniqueCancersPlot,top='(f)')
grid.arrange(uniqueCancersPlot)




fig_time <- arrangeGrob(ratePlot,rolePlot,sectionPlot,subsectionPlot,uniqueGenesPlot,uniqueCancersPlot,ncol=2)

grid.arrange(fig_time)
