library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)
library(data.table)
library(stringr)
library(ggplot2)
library(heatmaply)

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

citationTableExplanation <- "<br /><br /><br /><b><a name='citationtable'>Citation Table:</a></b><br />Select row in table above to see associated citations and sentences<br /><br />"

########################################################
# Aggregate Data for Collapsed Version (no gene roles) #
########################################################

collated_noroles <- collated[, list(sum(citation_count),paste(matching_id,collapse=',')), by=list(cancer_id,cancer_normalized,gene_hugo_id,gene_entrez_id,gene_normalized)]

names(collated_noroles)[names(collated_noroles) == "V1"] = "citation_count"
names(collated_noroles)[names(collated_noroles) == "V2"] = "matching_ids"

collated_noroles <- collated_noroles[order(collated_noroles$cancer_normalized),]
collated_noroles <- collated_noroles[order(collated_noroles$gene_normalized),]
collated_noroles <- collated_noroles[order(collated_noroles$citation_count,decreasing=T),]

#################################
# Aggregate data for pie charts #
#################################

cancerpiecounts <- aggregate(collated$citation_count,by=list(collated$cancer_normalized),FUN=sum)
colnames(cancerpiecounts) <- c('label','total_citation_count')
cancerpiecounts <- cancerpiecounts[order(cancerpiecounts$total_citation_count,decreasing=T),]

total <- sum(cancerpiecounts$total_citation_count)
cutoff <- (0.9/100.0) * total
cancerpiecounts$label <- as.character(cancerpiecounts$label)
cancerpiecounts[cancerpiecounts$total_citation_count<cutoff,'label'] <- 'other'
cancerpiecounts$label <- factor(cancerpiecounts$label, levels=c('other',as.character(cancerpiecounts[cancerpiecounts$total_citation_count>=cutoff,'label'])))

cancerpiecounts <- aggregate(cancerpiecounts$total_citation_count,by=list(cancerpiecounts$label),FUN=sum)
colnames(cancerpiecounts) <- c('label','total_citation_count')


genepiecounts <- aggregate(collated$citation_count,by=list(collated$gene_normalized),FUN=sum)
colnames(genepiecounts) <- c('label','total_citation_count')
genepiecounts <- genepiecounts[order(genepiecounts$total_citation_count,decreasing=T),]

total <- sum(genepiecounts$total_citation_count)
cutoff <- (0.9/100.0) * total
genepiecounts$label <- as.character(genepiecounts$label)
genepiecounts[genepiecounts$total_citation_count<cutoff,'label'] <- 'other'
genepiecounts$label <- factor(genepiecounts$label, levels=c('other',as.character(genepiecounts[genepiecounts$total_citation_count>=cutoff,'label'])))

genepiecounts <- aggregate(genepiecounts$total_citation_count,by=list(genepiecounts$label),FUN=sum)
colnames(genepiecounts) <- c('label','total_citation_count')

genesWithCitations <- dcast(collated, gene_normalized~role, sum, drop=F, value.var="citation_count")
genesWithCitations$total <- genesWithCitations$Driver + genesWithCitations$Oncogene + genesWithCitations$Tumor_Suppressor
genesWithCitations <- genesWithCitations[genesWithCitations$total>0,]
genesWithCitations <- genesWithCitations[order(genesWithCitations$total,decreasing=T),]
genesWithCitations$gene_link <- sprintf("<a href='?gene=%s' target='_blank'>%s</a>",genesWithCitations$gene_normalized,genesWithCitations$gene_normalized)
genesWithCitations$gene_normalized_lower <- str_to_lower(genesWithCitations$gene_normalized)
genelist_legend <- '<p></p><div style="display: inline-block; height:20px; width:20px; background-color:#66C2A5">&nbsp;</div> Driver <div style="display: inline-block; height:20px; width:20px; background-color:#8DA0CB">&nbsp;</div> Oncogene <div style="display: inline-block; height:20px; width:20px; background-color:#FC8D62">&nbsp;</div> Tumor Suppressor'

############################
# Prep Cancermine Profiles #
############################
source('profile_clustering.R')
cancermineProfiles <- generateProfiles(collated)

paper.clusteringTopCancerCount <- 25
paper.clusteringTopGeneRoleCount <- 25
cancerCounts <- aggregate(collated$citation_count, by=list(cancer_normalized=collated$cancer_normalized), FUN=sum)
colnames(cancerCounts) <- c('cancer_normalized','total_citation_count')
cancerCounts <- cancerCounts[order(cancerCounts$total_citation_count,decreasing=T),]
topCancers <- cancerCounts[1:paper.clusteringTopCancerCount,'cancer_normalized']

clusteringExplanation <- '<p><h3>CancerMine Profiles</h3></p><p>This plot shows the importance of strongly associated gene roles for the selected cancer types. The brighter the square, the more important that gene role is for that cancer type. It clusters similar cancer types together. You can change which cancer types to cluster by changing the choices below. For a summary of how the gene importance is calculated, see the FAQs on the help page.</p>'

ui <- function(req) {
  fluidPage(
    tags$head(
      includeHTML("google-analytics.js"),
      includeHTML("metadata.html"),
      tags$style(".rightAlign{float:right; margin-left:5px; margin-bottom: 20px;}"),
      tags$style(HTML('#clustering_update{background-color:#b3ccff}'))
    ),
    titlePanel("",windowTitle="CancerMine"),
    helpText(includeHTML("subheader.html")),
    tabsetPanel(id="maintabs",
                type = "tabs",
                tabPanel("About", 
                         includeHTML("landing-top.html"),
                         splitLayout(cellWidths = c("50%", "50%"), 
                                     plotlyOutput("gene_summary_piechart"),
                                     plotlyOutput("cancer_summary_piechart")),
                         includeHTML("landing-bottom.html"),
                         helpText(includeHTML("cc0.html"))),
                tabPanel("By Gene", 
                         sidebarPanel(
                           selectizeInput("gene_input", "Gene", geneNames, selected = 'EGFR', multiple = FALSE, options = list(maxOptions = 2*length(geneNames))),
                           plotlyOutput("gene_overview"),
                           checkboxInput("gene_collapseroles", "Collapse roles", FALSE),
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
                           selectizeInput("cancer_input", "Cancer", cancerNames, selected = 'diffuse large B-cell lymphoma', multiple = FALSE, options = list(maxOptions = 2*length(cancerNames))),
                           plotlyOutput("cancer_overview"),
                           checkboxInput("cancer_collapseroles", "Collapse roles", FALSE),
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
                tabPanel("With Gene List", 
                         includeHTML("genelist-top.html"),
                         textAreaInput("genelist_genes", "Genes:", value = "EGFR\nERBB2\nPTEN", height=200),
                         HTML(genelist_legend),
                         downloadButton("genelist_download", label = "Download", class='rightAlign'),
                         DT::dataTableOutput("genelist_resulttable"),
                         helpText(includeHTML("cc0.html"))),
                tabPanel("Clustering",
                         sidebarPanel(
                           HTML(clusteringExplanation),
                           selectInput("clustering_cancers","Cancers to Cluster",choices=levels(collated$cancer_normalized),multiple=T,selected=topCancers),
                           actionButton("clustering_clear", "Clear selection",width="45%"),
                           actionButton("clustering_usetopcancers", "Use 25 cancer list",width="45%"),
                           HTML("<p></p>"),
                           actionButton("clustering_update", "Update Plot",width="45%"),
                           width=3
                         ),
                         mainPanel(
                           plotlyOutput("clustering_output"))
                ),
                tabPanel("Help", 
                         includeHTML("help.html"),
                         helpText(includeHTML("cc0.html")))
    )
    
  )
}

input <- data.frame(gene_input='EGFR', cancer_input='prostate cancer', stringsAsFactors=FALSE)
matching_id <- '92aff62f03f4d4025638c3c771e2ab0c'

server <- function(input, output, session) {
  
  geneTableProxy = dataTableProxy('gene_table')
  
  url_event_val <- reactiveValues(count = 0)
  observe({
    if (url_event_val$count == 0) {
      
      query <- parseQueryString(session$clientData$url_search)
      if (!is.null(query[['gene']])) {
        geneName <- query[['gene']]
        updateSelectizeInput(session, "gene_input", selected = geneName)
        updateTabsetPanel(session, 'maintabs', selected = "By Gene")
        url_event_val$count <- 1
      }
      if (!is.null(query[['cancer']])) {
        cancerName <- query[['cancer']]
        updateSelectizeInput(session, "cancer_input", selected = cancerName)
        updateTabsetPanel(session, 'maintabs', selected = "By Cancer")
        url_event_val$count <- 1
      }
      if (!is.null(query[['genelist']])) {
        updateTabsetPanel(session, 'maintabs', selected = "With Gene List")
        url_event_val$count <- 1
      }
      if (!is.null(query[['clustering']])) {
        updateTabsetPanel(session, 'maintabs', selected = "Clustering")
        url_event_val$count <- 1
      }
      if (!is.null(query[['matching_id']])) {
        matching_id <- query[['matching_id']]
        warning(matching_id)
        mask <- collated$matching_id==matching_id
        geneName <- collated$gene_normalized[mask][1]
        warning(geneName)
        updateSelectizeInput(session, "gene_input", selected = geneName)
        updateTabsetPanel(session, 'maintabs', selected = "By Gene")
        
        table <- geneData()
        warning(paste('nrow(table) -',nrow(table)))
        rowNumber <- which(table$matching_id==matching_id)
        warning(paste('rowNumber -',rowNumber))
        
        currentPageNumber <- as.integer(input$gene_table_state$start / input$gene_table_state$length) + 1
        currentlySelectedRow <- input$gene_table_rows_selected
        
        perPage <- input$gene_table_state$length
        newPageNumber <- floor(rowNumber / perPage) + 1
        
        geneTableProxy %>% selectPage(newPageNumber)
        geneTableProxy %>% selectRows(rowNumber)
        
        #browseURL("#citationtable")
        if (newPageNumber == currentPageNumber && length(currentlySelectedRow) > 0 && rowNumber == currentlySelectedRow)
          url_event_val$count <- 1
      }
    }
  })
  
  geneData <- reactive({
    table <- collated[collated$gene_normalized==input$gene_input,]
    table
  })
  
  geneData_noRoles <- reactive({
    table <- collated_noroles[collated_noroles$gene_normalized==input$gene_input,]
    table
  })
  
  output$gene_table <- DT::renderDataTable({
    if (input$gene_collapseroles) {
      table <- geneData_noRoles()
      DT::datatable(table[,c('gene_normalized','cancer_normalized','citation_count')],
                    selection = 'single',
                    rownames = FALSE,
                    colnames=c('Gene', 'Cancer', 'Citation #'),
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
    } else {
      table <- geneData()
      DT::datatable(table[,c('role','gene_normalized','cancer_normalized','citation_count')],
                    selection = 'single',
                    rownames = FALSE,
                    colnames=c('Role','Gene', 'Cancer', 'Citation #'),
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
    }
  })
  
  
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
  
  gene_bar_event_val <- reactiveValues(count = 0)
  output$gene_barchart <- renderPlotly({
    table <- geneData()
    sourceName <- paste('gene_barchart_',gene_bar_event_val$count,sep='')
    if (nrow(table) > 0) {
      unmelted <- dcast(table, cancer_normalized~role, sum, drop=F, value.var="citation_count")
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$cancer_normalized <- factor(as.character(unmelted$cancer_normalized), unique(unmelted$cancer_normalized))
      
      p <- plot_ly(unmelted, x=~cancer_normalized, y=~Tumor_Suppressor, source=sourceName, type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
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
    sourceName <- paste('gene_barchart_',gene_bar_event_val$count,sep='')
    s <- event_data("plotly_click", source = sourceName)
    if (length(s) > 0) {
      if (input$gene_collapseroles) {
        table <- geneData_noRoles()
        rowNumber <- which(table$cancer_normalized==s$x)
      } else {
        table <- geneData()
        roleOptions <- c('Tumor_Suppressor','Oncogene','Driver')
        rowNumber <- which(table$cancer_normalized==s$x & table$role==roleOptions[s$curveNumber+1])
      }
      
      currentPageNumber <- as.integer(input$gene_table_state$start / input$gene_table_state$length) + 1
      currentlySelectedRow <- input$gene_table_rows_selected
      
      perPage <- input$gene_table_state$length
      newPageNumber <- floor(rowNumber / perPage) + 1
      
      if (newPageNumber != currentPageNumber)
        geneTableProxy %>% selectPage(newPageNumber)
      
      if (length(currentlySelectedRow) == 0 || rowNumber != currentlySelectedRow)
        geneTableProxy %>% selectRows(rowNumber)
      
      gene_bar_event_val$count <- gene_bar_event_val$count + 1
    }
  })
  
  output$gene_citations <- DT::renderDataTable({
    if(length(input$gene_table_rows_selected)>0) {
      if (input$gene_collapseroles) {
        table <- geneData_noRoles()
        row <- table[input$gene_table_rows_selected,]
        matching_ids <- strsplit(row$matching_ids,',')[[1]]
        entries <- sentences[sentences$matching_id %in% matching_ids,]
      } else {
        table <- geneData()
        row <- table[input$gene_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      }
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    
    if (input$gene_collapseroles) {
      DT::datatable(entries[,c('pubmed_link','role','journal_short','year','section','subsection','formatted_sentence'),],
                    selection = 'none',
                    rownames = FALSE,
                    colnames=c('PMID','Role','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                    escape = FALSE,
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
    } else {
      DT::datatable(entries[,c('pubmed_link','journal_short','year','section','subsection','formatted_sentence'),],
                    selection = 'none',
                    rownames = FALSE,
                    colnames=c('PMID','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                    escape = FALSE,
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
    }
  })
  
  
  output$gene_download_collated_all <- downloadHandler(
    filename = function() {
      return("cancermine_collated.tsv")
    },
    content = function(file) {
      if (input$gene_collapseroles) {
        outdata <- collated_noroles
      } else {
        outdata <- collated
      }
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$gene_download_collated_shown <- downloadHandler(
    filename = function() {
      return("cancermine_collated_subset.tsv")
    },
    content = function(file) {
      if (input$gene_collapseroles) {
        outdata <- geneData_noRoles()
      } else {
        outdata <- geneData()
      }
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
      if (input$gene_collapseroles) {
        table <- geneData_noRoles()
      } else {
        table <- geneData()
      }
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
        if (input$gene_collapseroles) {
          table <- geneData_noRoles()
          row <- table[input$gene_table_rows_selected,]
          matching_ids <- strsplit(row$matching_ids,',')[[1]]
          entries <- sentences[sentences$matching_id %in% matching_ids,]
        } else {
          table <- geneData()
          row <- table[input$gene_table_rows_selected,]
          entries <- sentences[sentences$matching_id==row$matching_id,]
        }
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  gene_event_val <- reactiveValues(count = 0)
  output$gene_summary_piechart <- renderPlotly({
    sourceName <- paste('gene_summary_piechart_',gene_event_val$count,sep='')
    p <- plot_ly(genepiecounts, labels = ~label, values = ~total_citation_count, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
      layout(title = paste('Genes'),
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
      config(displayModeBar = F) %>%
      layout(showlegend = FALSE)
    
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('gene_summary_piechart_',gene_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    if (length(d) > 0)
    {
      selected <- genepiecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "gene_input", selected = selected)
        updateTabsetPanel(session, 'maintabs', selected = "By Gene")
        
        gene_event_val$count <- gene_event_val$count + 1
      }
    }
  })
  
  
  
  
  
  
  
  
  
  
  cancerData <- reactive({
    table <- collated[collated$cancer_normalized==input$cancer_input,]
    table
  })
  
  cancerData_noRoles <- reactive({
    table <- collated_noroles[collated_noroles$cancer_normalized==input$cancer_input,]
    table
  })
  
  output$cancer_table <- DT::renderDataTable({
    if (input$cancer_collapseroles) {
      table <- cancerData_noRoles()
      DT::datatable(table[,c('gene_normalized','cancer_normalized','citation_count')],
                    selection = 'single',
                    rownames = FALSE,
                    colnames=c('Gene', 'Cancer', 'Citation #'),
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
    } else {
      table <- cancerData()
      DT::datatable(table[,c('role','gene_normalized','cancer_normalized','citation_count')],
                    selection = 'single',
                    rownames = FALSE,
                    colnames=c('Role','Gene', 'Cancer', 'Citation #'),
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30), stateSave = TRUE))
    }
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
  
  cancer_bar_event_val <- reactiveValues(count = 0)
  output$cancer_barchart <- renderPlotly({
    table <- cancerData()
    sourceName <- paste('cancer_barchart_',cancer_bar_event_val$count,sep='')
    if (nrow(table) > 0) {
      unmelted <- dcast(table, gene_normalized~role, sum, drop=F, value.var="citation_count")
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$gene_normalized <- factor(as.character(unmelted$gene_normalized), unique(unmelted$gene_normalized))
      
      p <- plot_ly(unmelted, x=~gene_normalized, y=~Tumor_Suppressor, source=sourceName, type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
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
    sourceName <- paste('cancer_barchart_',cancer_bar_event_val$count,sep='')
    s <- event_data("plotly_click", source = sourceName)
    if (length(s) > 0) {
      if (input$cancer_collapseroles) {
        table <- cancerData_noRoles()
        rowNumber <- which(table$gene_normalized==s$x)
      } else {
        table <- cancerData()
        roleOptions <- c('Tumor_Suppressor','Oncogene','Driver')
        rowNumber <- which(table$gene_normalized==s$x & table$role==roleOptions[s$curveNumber+1])
      }
      
      currentPageNumber <- as.integer(input$cancer_table_state$start / input$cancer_table_state$length) + 1
      currentlySelectedRow <- input$cancer_table_rows_selected
      
      perPage <- input$cancer_table_state$length
      newPageNumber <- floor(rowNumber / perPage) + 1
      
      if (newPageNumber != currentPageNumber)
        cancerTableProxy %>% selectPage(newPageNumber)
      
      if (length(currentlySelectedRow) == 0 || rowNumber != currentlySelectedRow)
        cancerTableProxy %>% selectRows(rowNumber)
      
      cancer_bar_event_val$count <- cancer_bar_event_val$count + 1
    }
  })
  
  
  output$cancer_citations <- DT::renderDataTable({
    if(length(input$cancer_table_rows_selected)>0) {
      if (input$cancer_collapseroles) {
        table <- cancerData_noRoles()
        row <- table[input$cancer_table_rows_selected,]
        matching_ids <- strsplit(row$matching_ids,',')[[1]]
        entries <- sentences[sentences$matching_id %in% matching_ids,]
      } else {
        table <- cancerData()
        row <- table[input$cancer_table_rows_selected,]
        entries <- sentences[sentences$matching_id==row$matching_id,]
      }
    } else {
      entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
      colnames(entries) <- colnames(sentences)
    }
    
    if (input$cancer_collapseroles) {
      DT::datatable(entries[,c('pubmed_link','role','journal_short','year','section','subsection','formatted_sentence'),],
                    selection = 'none',
                    rownames = FALSE,
                    colnames=c('PMID','Role','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                    escape = FALSE,
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))
    } else {
      DT::datatable(entries[,c('pubmed_link','journal_short','year','section','subsection','formatted_sentence'),],
                    selection = 'none',
                    rownames = FALSE,
                    colnames=c('PMID','Journal','Year', 'Section', 'Subsection', 'Sentence'),
                    escape = FALSE,
                    options = list(pageLength = 20, lengthMenu = c(10, 20, 30)))}
  })
  
  cancer_event_val <- reactiveValues(count = 0)
  output$cancer_summary_piechart <- renderPlotly({
    sourceName <- paste('cancer_summary_piechart_',cancer_event_val$count,sep='')
    p <- plot_ly(cancerpiecounts, labels = ~label, values = ~total_citation_count, type = 'pie', sort=FALSE, textinfo = 'label', textposition = 'inside', insidetextfont = list(color = '#FFFFFF'),source=sourceName) %>%
      layout(title = paste('Cancers'),
             xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
             yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
      config(displayModeBar = F) %>%
      layout(showlegend = FALSE)
    
    p$elementId <- NULL
    p
  })
  
  observe({
    sourceName <- paste('cancer_summary_piechart_',cancer_event_val$count,sep='')
    d <- event_data("plotly_click",source=sourceName)
    if (length(d) > 0)
    {
      selected <- cancerpiecounts[d$pointNumber+1,'label']
      if (!is.na(selected) && selected != 'other') {
        updateSelectizeInput(session, "cancer_input", selected = selected)
        updateTabsetPanel(session, 'maintabs', selected = "By Cancer")
        
        cancer_event_val$count <- cancer_event_val$count + 1
      }
    }
  })
  
  
  
  output$cancer_download_collated_all <- downloadHandler(
    filename = function() {
      return("cancermine_collated.tsv")
    },
    content = function(file) {
      if (input$cancer_collapseroles) {
        outdata <- collated_noroles
      } else {
        outdata <- collated
      }
      write.table(outdata, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  output$cancer_download_collated_shown <- downloadHandler(
    filename = function() {
      return("cancermine_collated_subset.tsv")
    },
    content = function(file) {
      if (input$cancer_collapseroles) {
        outdata <- cancerData_noRoles()
      } else {
        outdata <- cancerData()
      }
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
      if (input$cancer_collapseroles) {
        table <- cancerData_noRoles()
      } else {
        table <- cancerData()
      }
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
        if (input$cancer_collapseroles) {
          table <- cancerData_noRoles()
          row <- table[input$cancer_table_rows_selected,]
          matching_ids <- strsplit(row$matching_ids,',')[[1]]
          entries <- sentences[sentences$matching_id %in% matching_ids,]
        } else {
          table <- cancerData()
          row <- table[input$cancer_table_rows_selected,]
          entries <- sentences[sentences$matching_id==row$matching_id,]
        }
      } else {
        entries <- data.frame(matrix(nrow=0,ncol=ncol(sentences)))
        colnames(entries) <- colnames(sentences)
      }
      
      write.table(entries, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  
  
  
  
  
  
  
  
  
  genelistData <- reactive({
    selectedGenes = strsplit(str_to_lower(input$genelist_genes),'\n')[[1]]
    data <- genesWithCitations[genesWithCitations$gene_normalized_lower %in% selectedGenes,]
    data
  })
  
  output$genelist_resulttable <- DT::renderDataTable({
    data <- genelistData()

    maxvalue <- data[1,total]
    data$Driver_perc <- round(98*data$Driver/maxvalue,1)
    data$Oncogene_perc <- round(98*data$Oncogene/maxvalue,1)
    data$Tumor_Suppressor_perc <- round(98*data$Tumor_Suppressor/maxvalue,1)
    
    data$citation_bars <- sprintf("<div><div style='height:20px; width:%f%%; background-color: #66C2A5; float:left'>&nbsp;</div><div style='height:20px; width:%f%%; background-color: #8DA0CB; float:left'>&nbsp;</div><div style='height:20px; width:%f%%; background-color: #FC8D62; float:left'>&nbsp;</div></div>",data$Driver_perc, data$Oncogene_perc, data$Tumor_Suppressor_perc)

    DT::datatable(data[,c('gene_link','citation_bars','total')],
                  colnames=c('Gene','Roles','Total Citations'),
                  selection = 'none',
                  rownames = FALSE,
                  escape = FALSE,
                  options = list(paging=F,
                                 ordering=F,
                                 columnDefs = list(list(width = '80%', targets = c(1))),
                                 searching=F))
  })
  
  output$genelist_download <- downloadHandler(
    filename = function() {
      return("cancermine_genelist_results.tsv")
    },
    content = function(file) {
      data <- genelistData()
      selectedColumns <- c('gene_normalized','Driver','Oncogene','Tumor_Suppressor','total')
      data <- data[,selectedColumns,with=F]
      colnames(data) <- c('gene','driver_citations','oncogene_citations','tumor_suppressor_citations','all_citations')
      
      write.table(data, file, row.names = FALSE, sep='\t', quote=F)
    }
  )
  
  clusterGeneSet <- eventReactive(input$clustering_update, {
    input$clustering_cancers
  })
  
  output$clustering_output <- renderPlotly({
    if (input$clustering_update == 0) {
      cancers <- topCancers
    } else {
      cancers <- clusterGeneSet()
    }
    
    if (length(cancers) >= 2) {
      clusterHeatmapData <- selectDataForHeatmap(cancermineProfiles,cancers,paper.clusteringTopGeneRoleCount)
      p <- suppressMessages(suppressWarnings(
        heatmaply(clusterHeatmapData) %>% 
          layout(height = 800) %>% 
          config(displayModeBar = F)
      ))
    } else {
      p <- plotly_empty(type='pie') %>% config(displayModeBar = F)
    }
    p$elementId <- NULL
    suppressMessages(p)
    
  })
  
  
  observeEvent(input$clustering_clear, {
    updateSelectizeInput(session, "clustering_cancers", selected = F)
  })
  observeEvent(input$clustering_usetopcancers, {
    updateSelectizeInput(session, "clustering_cancers", selected = topCancers)
  })
}

shinyApp(ui, server)

