library(data.table)
library(DBI)

# Convert the sentences file into an SQLite database indexed by matching_id so that the Shiny app
# can look up sentences when needed instead of holding them all in memory

sentencesFilename <- 'cancermine_sentences.tsv'
sentencesDBFilename <- 'cancermine_sentences.sqlite'
columns <- c('matching_id','pmid','role','journal_short','year','section','subsection','formatted_sentence')

sentences <- fread(sentencesFilename,sep='\t',header=T,quote="",select=columns)

tempFilename <- paste(sentencesDBFilename,'.tmp',sep='')
if (file.exists(tempFilename))
  file.remove(tempFilename)

db <- dbConnect(RSQLite::SQLite(), tempFilename)
dbWriteTable(db, 'sentences', sentences)
invisible(dbExecute(db, 'CREATE INDEX idx_matching_id ON sentences (matching_id)'))
dbDisconnect(db)

stopifnot(file.rename(tempFilename, sentencesDBFilename))

print(paste("Success. Wrote", nrow(sentences), "sentences to", sentencesDBFilename))
