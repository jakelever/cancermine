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

collatedFilename <- 'cancermine_collated.tsv'
sentencesFilename <- 'cancermine_sentences.tsv'

collatedFilename <- normalizePath(collatedFilename)
sentencesFilename <- normalizePath(sentencesFilename)
fileInfo <- file.info(collatedFilename)
modifiedDate <- strsplit(as.character(fileInfo$mtime), ' ')[[1]][1]


collated <- fread(collatedFilename,sep='\t',header=T,stringsAsFactors=TRUE)
sentences <- fread(sentencesFilename,sep='\t',header=T,stringsAsFactors=TRUE)

# Fill in some details for the citation table
sentences$pubmed_link <- paste("<a target=\"_blank\" href='https://www.ncbi.nlm.nih.gov/pubmed/", sentences$pmid, "'>", sentences$pmid, "</a>", sep='')


geneNames <- sort(unique(as.character(collated$gene_normalized)))
cancerNames <- sort(unique(as.character(collated$cancer_normalized)))

colors <- brewer.pal(3,'Set2')
color_Driver <- colors[1]
color_TumorSuppressor <- colors[2]
color_Oncogene <- colors[3]

citationTableExplanation <- "<br /><br /><br /><b>Citation Table:</b><br />Select row in table above to see associated citations and sentences<br /><br />"

ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}")
    ),
    titlePanel("",windowTitle="CancerMine"),
    helpText(includeHTML("subheader.html")),
    tabsetPanel(type = "tabs",
                tabPanel("By Gene", 
                         sidebarPanel(
                           selectizeInput("gene_input", "Gene", geneNames, selected = 'EGFR', multiple = FALSE, options = list(maxOptions = 2*length(geneNames))),
                           plotlyOutput("gene_overview"),
                           width=3
                         ),
                         mainPanel(
                           plotlyOutput("gene_barchart"),
                           downloadButton("gene_download_collated_all", label = "Download All", class='rightAlign'),
                           downloadButton("gene_download_collated_shown", label = "Download Shown", class='rightAlign'),
                           DT::dataTableOutput("gene_table"),
                           
                           HTML(citationTableExplanation),
                           downloadButton("gene_download_sentences_all", label = "Download All Sentences", class='rightAlign'),
                           downloadButton("gene_download_sentences_above", label = "Download Sentences for this Gene", class='rightAlign'),
                           downloadButton("gene_download_sentences_shown", label = "Download Shown Sentences", class='rightAlign'),
                           DT::dataTableOutput("gene_citations"),
                           helpText(paste("Last updated on",modifiedDate)),
                           helpText(includeHTML("cc0.html"))
                         )
                ),
                tabPanel("By Cancer", 
                         sidebarPanel(
                           selectizeInput("cancer_input", "Cancer", cancerNames, selected = 'colorectal cancer', multiple = FALSE, options = list(maxOptions = 2*length(cancerNames))),
                           plotlyOutput("cancer_overview"),
                           width=3
                         ),
                         mainPanel(
                           plotlyOutput("cancer_barchart"),
                           downloadButton("cancer_download_collated_all", label = "Download All", class='rightAlign'),
                           downloadButton("cancer_download_collated_shown", label = "Download Shown", class='rightAlign'),
                           DT::dataTableOutput("cancer_table"),
                           
                           HTML(citationTableExplanation),
                           downloadButton("cancer_download_sentences_all", label = "Download All Sentences", class='rightAlign'),
                           downloadButton("cancer_download_sentences_above", label = "Download Sentences for this Cancer", class='rightAlign'),
                           downloadButton("cancer_download_sentences_shown", label = "Download Shown Sentences", class='rightAlign'),
                           DT::dataTableOutput("cancer_citations"),
                           helpText(paste("Last updated on",modifiedDate)),
                           helpText(includeHTML("cc0.html"))
                         )
                ),
                tabPanel("Help", helpText(includeHTML("help.html"))),
                tabPanel("About", helpText(includeHTML("about.html")))
    )
    
  )
}

input <- data.frame(gene_input='EGFR', cancer_input='prostate cancer', stringsAsFactors=FALSE)

server <- function(input, output, session) {
  geneData <- reactive({
    table <- collated[collated$gene_normalized==input$gene_input,]
    if (nrow(table)>0) {
      ordering <- aggregate(table$citation_count,by=list(table$cancer_normalized),FUN=sum)
      colnames(ordering) <- c('cancer_normalized','citation_count_for_cancer')
      table <- dplyr::inner_join(table,ordering,by='cancer_normalized')
      table <- table[order(table$citation_count_for_cancer,decreasing=TRUE),]
      table <- table[,colnames(table) != 'citation_count_for_cancer']
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$gene_table <- DT::renderDataTable({
    table <- geneData()
    DT::datatable(table[,c('role','gene_normalized','cancer_normalized','citation_count')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Role','Gene', 'Cancer', 'Citation #'),
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  geneTableProxy = dataTableProxy('gene_table')
  
  
  output$gene_overview <- renderPlotly({
    table <- geneData()
    if (nrow(table) > 0) {
      relationcounts <- aggregate(table$citation_count,by=list(table$role),FUN=sum)
      colnames(relationcounts) <- c('role','citation_count_for_role')
      
      completecounts <- data.frame(role=levels(relationcounts$role),citation_count_for_role=0)
      rownames(completecounts) <- completecounts$role
      completecounts[as.character(relationcounts$role),'citation_count_for_role'] <- relationcounts$citation_count_for_role
      
      completecounts$color <- "#000000"
      completecounts[completecounts$role=='Driver','color'] <- color_Driver
      completecounts[completecounts$role=='Oncogene','color'] <- color_Oncogene
      completecounts[completecounts$role=='Tumor_Suppressor','color'] <- color_TumorSuppressor
      
      p <- plot_ly(completecounts, x=~role, y=~citation_count_for_role, source='gene_barchart', type = 'bar', marker=list(color=completecounts$color)) %>%
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
      unmelted <- dcast(table, cancer_normalized~role, sum, drop=F, value.var="citation_count")
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
      roleOptions <- c('Tumor_Suppressor','Oncogene','Driver')
      rowNumber <- rownames(table[table$cancer_normalized==s$x & table$role==roleOptions[s$curveNumber+1],])
      geneTableProxy %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  output$gene_citations <- DT::renderDataTable({
    if(length(input$gene_table_rows_selected)>0) {
      table <- geneData()
      row <- table[input$gene_table_rows_selected,]
      entries <- sentences[sentences$matching_id==row$matching_id,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    DT::datatable(entries[,c('pubmed_link','journal_short','year','section','subsection','formatted_sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  
  output$gene_download_collated_all <- downloadHandler(
    filename = function() {
      return("cancermine_collated.tsv")
    },
    content = function(file) {
      outdata <- collated
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_shown <- downloadHandler(
    filename = function() {
      return("cancermine_collated_subset.tsv")
    },
    content = function(file) {
      outdata <- geneData()
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  output$gene_download_sentences_all <- downloadHandler(
    filename = function() {
      return("cancermine_sentences.tsv")
    },
    content = function(file) {
      write.table(sentences, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_sentences_above <- downloadHandler(
    filename = function() {
      return("cancermine_sentences_gene.tsv")
    },
    content = function(file) {
      table <- geneData()
      entries <- sentences[sentences$matching_id %in% table$matching_id,]
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_sentences_shown <- downloadHandler(
    filename = function() {
      return("cancermine_sentences_subset.tsv")
    },
    content = function(file) {
      if(length(input$gene_table_rows_selected)>0) {
        table <- geneData()
        row <- table[input$gene_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
  
  
  
  
  
  
  cancerData <- reactive({
    table <- collated[collated$cancer_normalized==input$cancer_input,]
    if (nrow(table)>0) {
      ordering <- aggregate(table$citation_count,by=list(table$gene_normalized),FUN=sum)
      colnames(ordering) <- c('gene_normalized','citation_count_for_gene')
      table <- dplyr::inner_join(table,ordering,by='gene_normalized')
      table <- table[order(table$citation_count_for_gene,decreasing=TRUE),]
      table <- table[,colnames(table) != 'citation_count_for_gene']
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$cancer_table <- DT::renderDataTable({
    table <- cancerData()
    DT::datatable(table[,c('role','gene_normalized','cancer_normalized','citation_count')],
                  selection = 'single',
                  rownames = FALSE,
                  colnames=c('Role','Gene', 'Cancer', 'Citation #'),
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  cancerTableProxy = dataTableProxy('cancer_table')
  
  
  output$cancer_overview <- renderPlotly({
    table <- cancerData()
    if (nrow(table) > 0) {
      relationcounts <- aggregate(table$citation_count,by=list(table$role),FUN=sum)
      colnames(relationcounts) <- c('role','citation_count_for_role')
      
      completecounts <- data.frame(role=levels(relationcounts$role),citation_count=0)
      rownames(completecounts) <- completecounts$role
      completecounts[as.character(relationcounts$role),'citation_count_for_role'] <- relationcounts$citation_count_for_role
      
      completecounts$color <- "#000000"
      completecounts[completecounts$role=='Driver','color'] <- color_Driver
      completecounts[completecounts$role=='Oncogene','color'] <- color_Oncogene
      completecounts[completecounts$role=='Tumor_Suppressor','color'] <- color_TumorSuppressor
      
      p <- plot_ly(completecounts, x=~role, y=~citation_count_for_role, source='gene_barchart', type = 'bar', marker=list(color=completecounts$color)) %>%
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
      unmelted <- dcast(table, gene_normalized~role, sum, drop=F, value.var="citation_count")
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
      roleOptions <- c('Tumor_Suppressor','Oncogene','Driver')
      rowNumber <- rownames(table[table$gene_normalized==s$x & table$role==roleOptions[s$curveNumber+1],])
      cancerTableProxy %>% selectRows(as.numeric(rowNumber))
    }
  })
  
  
  output$cancer_citations <- DT::renderDataTable({
    if(length(input$cancer_table_rows_selected)>0) {
      table <- cancerData()
      row <- table[input$cancer_table_rows_selected,]
      entries <- sentences[sentences$matching_id==row$matching_id,]
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    DT::datatable(entries[,c('pubmed_link','journal_short','year','section','subsection','formatted_sentence'),],
                  selection = 'none',
                  rownames = FALSE,
                  colnames=c('PMID','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                  escape = FALSE,
                  options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
  })
  
  
  
  output$cancer_download_collated_all <- downloadHandler(
    filename = function() {
      return("cancermine_collated.tsv")
    },
    content = function(file) {
      outdata <- collated
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$cancer_download_collated_shown <- downloadHandler(
    filename = function() {
      return("cancermine_collated_subset.tsv")
    },
    content = function(file) {
      outdata <- cancerData()
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  output$cancer_download_sentences_all <- downloadHandler(
    filename = function() {
      return("cancermine_sentences.tsv")
    },
    content = function(file) {
      write.table(sentences, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$cancer_download_sentences_above <- downloadHandler(
    filename = function() {
      return("cancermine_sentences_cancer.tsv")
    },
    content = function(file) {
      table <- cancerData()
      entries <- sentences[sentences$matching_id %in% table$matching_id,]
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$cancer_download_sentences_shown <- downloadHandler(
    filename = function() {
      return("cancermine_sentences_subset.tsv")
    },
    content = function(file) {
      if(length(input$cancer_table_rows_selected)>0) {
        table <- cancerData()
        row <- table[input$cancer_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
}

shinyApp(ui, server)
