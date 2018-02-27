library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(stringr)

# Weird hack as R sometimes "forgets" its working directory
wd <- setwd(".")
setwd(wd)

# Make an empty Google analytics file (for dev version - not for production)
if (!file.exists('google-analytics.js'))
{
  file.create('google-analytics.js')
}


cancermineFilename <- 'cancermine_latest.tsv'
cancermineFilename <- normalizePath(cancermineFilename)
fileInfo <- file.info(cancermineFilename)
modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]

cancermine <- read.table(cancermineFilename,header=T,sep='\t',quote='',comment.char='')
cancermine <- as.data.table(cancermine)

# Remove citation with missing PMID (for whatever reason)
cancermine <- cancermine[cancermine$pmid!='None',]
cancermine$pmid <- as.integer(cancermine$pmid)

# Remove the entity location columns and unique the rows
nonLocationColumns <- grep("(start|end)",colnames(cancermine),invert=T)
cancermine <- cancermine[,nonLocationColumns,with=FALSE]
cancermine <- cancermine[!duplicated(cancermine),]

# Fill in some details for the citation table
cancermine$pmidLink <- paste("<a target=\"_blank\" href='https://www.ncbi.nlm.nih.gov/pubmed/", cancermine$pmid, "'>", cancermine$pmid, "</a>", sep='')
cancermine$journalShort <- strtrim(cancermine$journal,51)
cancermine[str_length(cancermine$journalShort)==51,'journalShort'] <- paste(strtrim(cancermine[str_length(cancermine$journalShort)==51,]$journalShort,50),'...',sep='')
cancermine$journalShort <- factor(cancermine$journalShort)


# Do some ordering
cancermine <- cancermine[order(cancermine$relationtype),]
cancermine <- cancermine[order(cancermine$cancer_normalized),]
cancermine <- cancermine[order(cancermine$gene_normalized),]
cancermine <- cancermine[order(cancermine$year,decreasing=T),]

cancermineCounts <- plyr::count(cancermine[,c('relationtype','gene_normalized','cancer_normalized')])
cancermineCounts <- as.data.table(cancermineCounts)
cancermineCounts <- cancermineCounts[order(cancermineCounts$relationtype),]
cancermineCounts <- cancermineCounts[order(cancermineCounts$cancer_normalized),]
cancermineCounts <- cancermineCounts[order(cancermineCounts$gene_normalized),]

cancermine$preparedText <- paste(cancermine$sentence, " <a href='https://www.ncbi.nlm.nih.gov/pubmed/", cancermine$pmid, "'>PMID:", cancermine$pmid, "</a>", sep='')

genecounts <- plyr::count(cancermine$gene_normalized)
genecounts <- genecounts[order(genecounts$freq),]

cancercounts <- plyr::count(cancermine$cancer_normalized)
cancercounts <- cancercounts[order(cancercounts$freq),]

geneNames <- sort(unique(as.character(cancermine$gene_normalized)))
cancerNames <- sort(unique(as.character(cancermine$cancer_normalized)))

colors <- brewer.pal(3,'Set2')
color_Driver <- colors[1]
color_TumorSuppressor <- colors[2]
color_Oncogene <- colors[3]

citationTableExplanation <- "<br /><br /><br /><b>Citation Table:</b><br />Select row in table above to see associated citations and sentences<br /><br />"

ui <- function(req) {
  fluidPage(
  tags$head(includeHTML("google-analytics.js")),
  titlePanel("",windowTitle="CancerMine"),
  helpText(includeHTML("subheader.html")),
  tabsetPanel(type = "tabs",
              tabPanel("By Gene", 
                       sidebarPanel(
                         selectizeInput("gene_input", "Gene", geneNames, selected = 'EGFR', multiple = FALSE, options = list(maxOptions = 2*length(geneNames))),
                         plotlyOutput("gene_overview")
                         ),
                       mainPanel(
                         plotlyOutput("gene_barchart"),
                         DT::dataTableOutput("gene_table"),
                         HTML(citationTableExplanation),
                         DT::dataTableOutput("gene_citations"),
                         helpText(paste("Last updated on",modifiedDate))
                         )
              ),
              tabPanel("By Cancer", 
                       sidebarPanel(
                         selectizeInput("cancer_input", "Cancer", cancerNames, selected = 'colorectal cancer', multiple = FALSE, options = list(maxOptions = 2*length(cancerNames))),
                         plotlyOutput("cancer_overview")
                       ),
                       mainPanel(
                         plotlyOutput("cancer_barchart"),
                         DT::dataTableOutput("cancer_table"),
                         HTML(citationTableExplanation),
                         DT::dataTableOutput("cancer_citations"),
                         helpText(paste("Last updated on",modifiedDate))
                       )
              ),
              tabPanel("Help", helpText(includeHTML("help.html"))),
              tabPanel("About", helpText(includeHTML("about.html")))
  ),
  helpText(includeHTML("cc0.html"))
  
)
}

input <- data.frame(gene_input='EGFR', stringsAsFactors=FALSE)

server <- function(input, output, session) {
  geneData <- reactive({
    table <- cancermineCounts[cancermineCounts$gene_normalized==input$gene_input,c('relationtype','gene_normalized','cancer_normalized','freq')]
    if (nrow(table)>0) {
      ordering <- aggregate(table$freq,by=list(table$cancer_normalized),FUN=sum)
      colnames(ordering) <- c('cancer_normalized','totalfreq')
      table <- dplyr::inner_join(table,ordering,by='cancer_normalized')
      table <- table[order(table$totalfreq,decreasing=TRUE),]
      table <- table[,c('relationtype','gene_normalized','cancer_normalized','freq')]
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$gene_table <- DT::renderDataTable({
    DT::datatable(geneData(),
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Role','Gene', 'Cancer', 'Sentence count'),
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  geneTableProxy = dataTableProxy('gene_table')
  
  
  output$gene_overview <- renderPlotly({
    table <- geneData()
    if (nrow(table) > 0) {
      relationcounts <- aggregate(table$freq,by=list(table$relationtype),FUN=sum)
      colnames(relationcounts) <- c('relationtype','freq')
      
      completecounts <- data.frame(relationtype=levels(relationcounts$relationtype),freq=0)
      rownames(completecounts) <- completecounts$relationtype
      completecounts[as.character(relationcounts$relationtype),'freq'] <- relationcounts$freq
      
      completecounts$color <- "#000000"
      completecounts[completecounts$relationtype=='Driver','color'] <- color_Driver
      completecounts[completecounts$relationtype=='Oncogene','color'] <- color_Oncogene
      completecounts[completecounts$relationtype=='Tumor_Suppressor','color'] <- color_TumorSuppressor
      
      completecounts <- completecounts[completecounts$freq>0,]
      
      p <- plot_ly(completecounts, x=~relationtype, y=~freq, source='gene_barchart', type = 'bar', marker=list(color=completecounts$color)) %>%
        layout(yaxis = list(title = 'Total Citation Count'), margin = list(b = 200), xaxis=list(title = "", tickangle = 90))%>% 
        config(displayModeBar = F)
      
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  output$gene_barchart <- renderPlotly({
    table <- geneData()
    if (nrow(table) > 0) {
      unmelted <- dcast(table, cancer_normalized~relationtype, sum, drop=F, value.var="freq")
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$cancer_normalized <- factor(as.character(unmelted$cancer_normalized), unique(unmelted$cancer_normalized))
      
      p <- plot_ly(unmelted, x=~cancer_normalized, y=~Tumor_Suppressor, source='gene_barchart', type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
        add_trace(y = ~Oncogene, name = 'Oncogene', marker=list(color=color_Oncogene)) %>%
        add_trace(y = ~Driver, name = 'Driver', marker=list(color=color_Driver))%>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
      
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  observe({
    s <- event_data("plotly_click", source = "gene_barchart")
    if (length(s) > 0) {
      table <- geneData()
      relationtypeOptions <- c('Tumor_Suppressor','Oncogene','Driver')
      rowNumber <- rownames(table[table$cancer_normalized==s$x & table$relationtype==relationtypeOptions[s$curveNumber+1],])
      geneTableProxy %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  output$gene_citations <- DT::renderDataTable({
    if(length(input$gene_table_rows_selected)>0) {
      table <- geneData()
      row <- table[input$gene_table_rows_selected,]
      entries <- cancermine[cancermine$relationtype==row$relationtype & cancermine$cancer_normalized==row$cancer_normalized & cancermine$gene_normalized==row$gene_normalized,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=5))
      colnames(entries) <- c('pmidLink','journalShort','year','section','sentence')
    }
    DT::datatable(entries[,c('pmidLink','journalShort','year','section','sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'Section', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20))
  })
  
  
  
  
  
  
  
  
  cancerData <- reactive({
    table <- cancermineCounts[cancermineCounts$cancer_normalized==input$cancer_input,c('relationtype','gene_normalized','cancer_normalized','freq')]
    if (nrow(table)>0) {
      ordering <- aggregate(table$freq,by=list(table$gene_normalized),FUN=sum)
      colnames(ordering) <- c('gene_normalized','totalfreq')
      table <- dplyr::inner_join(table,ordering,by='gene_normalized')
      table <- table[order(table$totalfreq,decreasing=TRUE),]
      table <- table[,c('relationtype','gene_normalized','cancer_normalized','freq')]
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$cancer_table <- DT::renderDataTable({
    DT::datatable(cancerData(),
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Role','Gene', 'Cancer', 'Sentence count'),
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  cancerTableProxy = dataTableProxy('cancer_table')
  
  
  output$cancer_overview <- renderPlotly({
    table <- cancerData()
    if (nrow(table) > 0) {
      relationcounts <- aggregate(table$freq,by=list(table$relationtype),FUN=sum)
      colnames(relationcounts) <- c('relationtype','freq')
      
      completecounts <- data.frame(relationtype=levels(relationcounts$relationtype),freq=0)
      rownames(completecounts) <- completecounts$relationtype
      completecounts[as.character(relationcounts$relationtype),'freq'] <- relationcounts$freq
      
      completecounts$color <- "#000000"
      completecounts[completecounts$relationtype=='Driver','color'] <- color_Driver
      completecounts[completecounts$relationtype=='Oncogene','color'] <- color_Oncogene
      completecounts[completecounts$relationtype=='Tumor_Suppressor','color'] <- color_TumorSuppressor
      
      completecounts <- completecounts[completecounts$freq>0,]
      
      p <- plot_ly(completecounts, x=~relationtype, y=~freq, source='gene_barchart', type = 'bar', marker=list(color=completecounts$color)) %>%
        layout(yaxis = list(title = 'Total Citation Count'), margin = list(b = 200), xaxis=list(title = "", tickangle = 90))%>% 
        config(displayModeBar = F)
      
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
  })
  
  output$cancer_barchart <- renderPlotly({
    table <- cancerData()
    if (nrow(table) > 0) {
      unmelted <- dcast(table, gene_normalized~relationtype, sum, drop=F, value.var="freq")
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$gene_normalized <- factor(as.character(unmelted$gene_normalized), unique(unmelted$gene_normalized))
      
      p <- plot_ly(unmelted, x=~gene_normalized, y=~Tumor_Suppressor, source='cancer_barchart', type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
        add_trace(y = ~Oncogene, name = 'Oncogene', marker=list(color=color_Oncogene)) %>%
        add_trace(y = ~Driver, name = 'Driver', marker=list(color=color_Driver))%>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
    } else {
      p <- plotly_empty(type='pie')%>% 
        config(displayModeBar = F)
    }
    p$elementId <- NULL
    p
    
  })
  
  observe({
    s <- event_data("plotly_click", source = "cancer_barchart")
    if (length(s) > 0) {
      table <- cancerData()
      relationtypeOptions <- c('Tumor_Suppressor','Oncogene','Driver')
      rowNumber <- rownames(table[table$gene_normalized==s$x & table$relationtype==relationtypeOptions[s$curveNumber+1],])
      cancerTableProxy %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  
  output$cancer_text <- renderPrint({
    if(length(input$cancer_table_rows_selected)==0) {
      cat("Select a row from the table to see the associated sentences and citations")
    } else {
      table <- cancerData()
      row <- table[input$cancer_table_rows_selected,]
      entries <- cancermine[cancermine$relationtype==row$relationtype & cancermine$cancer_normalized==row$cancer_normalized & cancermine$gene_normalized==row$gene_normalized,]
      
      cat(entries$preparedText,sep='<br /><br />')
    }
  })
  
  
  output$cancer_citations <- DT::renderDataTable({
    if(length(input$cancer_table_rows_selected)>0) {
      table <- cancerData()
      row <- table[input$cancer_table_rows_selected,]
      entries <- cancermine[cancermine$relationtype==row$relationtype & cancermine$cancer_normalized==row$cancer_normalized & cancermine$gene_normalized==row$gene_normalized,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=5))
      colnames(entries) <- c('pmidLink','journalShort','year','section','sentence')
    }
    DT::datatable(entries[,c('pmidLink','journalShort','year','section','sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'Section', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20))
  })
  
  
  
  
  
  
  
  
  
  
  
}

shinyApp(ui, server)
