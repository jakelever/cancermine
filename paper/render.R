library(bookdown)
library(tinytex)
library(data.table)

bookdown::render_book("paper.Rmd", "bookdown::pdf_book", clean = T)
bookdown::render_book("paper.Rmd", "bookdown::word_document2")

tinytex::pdflatex('_book/thesis.tex')
#tinytex::pdflatex('thesis.tex')
