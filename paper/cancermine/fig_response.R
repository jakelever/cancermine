source('cancermine/dependencies.R')

sentencesFilename <- 'cancermine/cancermine_unfiltered.tsv'
sentencesFilename <- normalizePath(sentencesFilename)

cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
cancermineSentences <- cancermineSentences[!is.na(cancermineSentences$gene_entrez_id),]


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

#paper.notnovel2018Perc = round(100*notnovel2018/(novel2018+notnovel2018))
#paper.novel2018Perc = round(100*novel2018/(novel2018+notnovel2018))

#paper.rolesIn2018 <- sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2018,'freq'])
#paper.avgRolesPerMonth <- round(sum(novelAndNotNovelYearCounts[novelAndNotNovelYearCounts$year==2018,'freq'])/12)
paper.avgNovelRolesPerMonthWithoutThreshold <- round(sum(novel2018/12))
#paper.rolesIn2018 <- prettyNum(paper.rolesIn2018,big.mark=',')

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
ratePlotLowThreshold <- barchart(freq ~ yearTxt, 
                     novelAndNotNovelYearCounts, 
                     groups=novel,
                     ylab="Mentions per year",
                     par.settings = my.settings,
                     auto.key=list(space="top", columns=2, 
                                   points=FALSE, rectangles=TRUE),
                     scales=list(x=list(rot=45,labels=labels)),
                     stack=T,
                     horizontal=F)
ratePlotLowThreshold = arrangeGrob(ratePlotLowThreshold,top='(a) Low Thresholds')


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
ratePlot = arrangeGrob(ratePlot,top='(b) High Thresholds')




rateResponsePlot <- arrangeGrob(ratePlotLowThreshold,ratePlot)
grid.arrange(rateResponsePlot)


unfilteredFilename <- 'cancermine/cancermine_unfiltered.tsv'
unfilteredFilename <- normalizePath(unfilteredFilename)

cancermineUnfiltered <- read.table(unfilteredFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
cancermineUnfiltered <- cancermineUnfiltered[!is.na(cancermineUnfiltered$gene_entrez_id),]
cancermineUnfiltered <- cancermineUnfiltered[cancermineUnfiltered$year<2019,]



cancermineUnfiltered$predictprobGroup <- paste(as.character(floor(cancermineUnfiltered$predictprob*10)/10),'-',as.character(floor(cancermineUnfiltered$predictprob*10)/10+0.1))

predictProbGroupCounts <- plyr::count(cancermineUnfiltered[,c('predictprobGroup','year'),drop=F])
predictProbGroupCounts <- predictProbGroupCounts[order(predictProbGroupCounts$year),]
predictProbGroupCounts$yearTxt <- as.character(predictProbGroupCounts$year)

labels <- sort(unique(predictProbGroupCounts$year))
labels[labels%%5!=0] <- ''

barchart(freq ~ yearTxt, 
         predictProbGroupCounts, 
         groups=predictprobGroup,
         ylab="Mentions per year",
         auto.key=list(space="top", columns=2, 
                       points=FALSE, rectangles=TRUE),
         scales=list(x=list(rot=45,labels=labels)),
         stack=T,
         horizontal=F)


lotsComparisonsPlot <- arrangeGrob(arrangeGrob(fig_comparisons_lowthreshold,top="(a) Low Thresholds"),
                                   arrangeGrob(fig_comparisons,top="(b) High Thresholds"))
grid.arrange(lotsComparisonsPlot)






sentencesFilename <- 'cancermine/cancermine_sentences.tsv'
sentencesFilename <- normalizePath(sentencesFilename)

cancermineSentences <- read.table(sentencesFilename,header=T,sep='\t',quote='',comment.char='',encoding="UTF-8")
cancermineSentences <- cancermineSentences[!is.na(cancermineSentences$gene_entrez_id),]


cancermineSentences <- cancermineSentences[order(cancermineSentences$year),]
cancermineSentences$isNovel <- !duplicated(cancermineSentences[,c('role','gene_hugo_id','cancer_id')])

#search <- cancermineSentences[cancermineSentences$isNovel==TRUE & cancermineSentences$subsection=='results',c('year','formatted_sentence')]

search <- cancermineSentences[cancermineSentences$isNovel==TRUE & cancermineSentences$subsection=='introduction',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=F,sep='\t',row.names=F,col.names=T,quote=F)

search <- cancermineSentences[cancermineSentences$isNovel==FALSE & cancermineSentences$subsection=='introduction',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=T,sep='\t',row.names=F,col.names=F,quote=F)

search <- cancermineSentences[cancermineSentences$isNovel==TRUE & cancermineSentences$subsection=='results',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=T,sep='\t',row.names=F,col.names=F,quote=F)

search <- cancermineSentences[cancermineSentences$isNovel==FALSE & cancermineSentences$subsection=='results',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=T,sep='\t',row.names=F,col.names=F,quote=F)

search <- cancermineSentences[cancermineSentences$isNovel==TRUE & cancermineSentences$subsection=='discussion',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=T,sep='\t',row.names=F,col.names=F,quote=F)

search <- cancermineSentences[cancermineSentences$isNovel==FALSE & cancermineSentences$subsection=='discussion',c('isNovel','subsection','role','gene_normalized','cancer_normalized','pmid','formatted_sentence')]
write.table(tail(search,5),file='novelExamples.txt',append=T,sep='\t',row.names=F,col.names=F,quote=F)

