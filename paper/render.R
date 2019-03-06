library(bookdown)
library(tinytex)
library(data.table)
library(stringr)

if (file.exists('_main.Rmd'))
  file.remove('_main.Rmd')

figureSearch <- str_match_all(paste(readLines('supplementary.Rmd'),collapse='\n'),'supfig\\w+')
figureNames <- unique(figureSearch[[1]][,1])

tableSearch <- str_match_all(paste(readLines('supplementary.Rmd'),collapse='\n'),'suptab\\w+')
tableNames <- unique(tableSearch[[1]][,1])

mainText <- paste(readLines('brief_communication.Rmd'),collapse='\n')
if (length(figureNames) > 0) {
  for (i in 1:length(figureNames)) {
    pattern <- sprintf("\\@ref(fig:%s)",figureNames[i])
    mainText <- sub(pattern,i,mainText,fixed=T)
  }
}
if (length(tableNames) > 0) {
  for (i in 1:length(tableNames)) {
    pattern <- sprintf("\\@ref(tab:%s)",tableNames[i])
    mainText <- sub(pattern,i,mainText,fixed=T)
  }
}
cat(mainText,file='brief_communication.tmp.Rmd')

bookdown::render_book("brief_communication.tmp.Rmd", "bookdown::word_document2", preview=T)
file.rename('_book/_main.docx','_book/brief_communication.docx')
bookdown::render_book("methods.Rmd", "bookdown::word_document2", preview=T)
file.rename('_book/_main.docx','_book/methods.docx')
bookdown::render_book("supplementary.Rmd", "bookdown::word_document2", preview=T)
file.rename('_book/_main.docx','_book/supplementary.docx')

if (file.exists('brief_communication.tmp.Rmd'))
  file.remove('brief_communication.tmp.Rmd')
