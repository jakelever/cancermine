library(bookdown)
library(tinytex)
library(data.table)
library(stringr)

figureSearch <- str_match_all(paste(readLines('supplementary.Rmd'),collapse='\n'),'ref:(\\w+)')
figureNames <- unique(figureSearch[[1]][,2])

mainText <- paste(readLines('correspondence.Rmd'),collapse='\n')
for (i in 1:length(figureNames)) {
  pattern <- sprintf("\\@ref(fig:%s)",figureNames[i])
  mainText <- sub(pattern,i,mainText,fixed=T)
}
cat(mainText,file='correspondence.tmp.Rmd')

bookdown::render_book("correspondence.tmp.Rmd", "bookdown::word_document2", preview=T)
file.rename('_book/_main.docx','_book/correspondence.docx')
bookdown::render_book("supplementary.Rmd", "bookdown::word_document2", preview=T)
file.rename('_book/_main.docx','_book/supplementary.docx')

