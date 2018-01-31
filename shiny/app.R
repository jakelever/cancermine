library(plotly)
library(shiny)
library(DT)
library(plyr)
library(dplyr)
library(reshape2)
library(RColorBrewer)

cancermine <- read.table('cancermine_latest.tsv',header=T,sep='\t',quote='',comment.char='')

# Remove the entity location columns and unique the rows
nonLocationColumns <- grep("(start|end)",colnames(cancermine),invert=T,value=T)
cancermine <- cancermine[,nonLocationColumns]
cancermine <- cancermine[!duplicated(cancermine),]

# Do some ordering
cancermine <- cancermine[order(cancermine$relationtype),]
cancermine <- cancermine[order(cancermine$cancer_standardized),]
cancermine <- cancermine[order(cancermine$gene_standardized),]

cancermineCounts <- plyr::count(cancermine[,c('relationtype','gene_standardized','cancer_standardized')])
cancermineCounts <- cancermineCounts[order(cancermineCounts$relationtype),]
cancermineCounts <- cancermineCounts[order(cancermineCounts$cancer_standardized),]
cancermineCounts <- cancermineCounts[order(cancermineCounts$gene_standardized),]

cancermine$preparedText <- paste(cancermine$sentence, " <a href='https://www.ncbi.nlm.nih.gov/pubmed/", cancermine$pmid, "'>", cancermine$pmid, "</a>", sep='')

#DF<-ddply(data,.(gene_standardized,cancer_standardized),transform,prop=number/sum(number))

genecounts <- plyr::count(cancermine$gene_standardized)
genecounts <- genecounts[order(genecounts$freq),]

cancercounts <- plyr::count(cancermine$cancer_standardized)
cancercounts <- cancercounts[order(cancercounts$freq),]

geneNames <- sort(unique(as.character(cancermine$gene_standardized)))
cancerNames <- sort(unique(as.character(cancermine$cancer_standardized)))

colors <- brewer.pal(3,'Set2')
color_Driver <- colors[1]
color_TumorSuppressor <- colors[2]
color_Oncogene <- colors[3]

ui <- fluidPage(
  headerPanel("CancerMine"),
  helpText("Text mined database of drivers, oncogenes and tumor suppressors in cancer"),
  tabsetPanel(type = "tabs",
              tabPanel("By Gene", 
                       sidebarPanel(
                         selectizeInput("gene_input", "Gene", geneNames, selected = 'EGFR', multiple = FALSE, options = NULL),
                         plotlyOutput("gene_piechart"),
                         htmlOutput("gene_text")
                         #verbatimTextOutput("gene_text")
                         ),
                       mainPanel(
                         plotlyOutput("gene_barchart"),
                         DT::dataTableOutput("gene_table")
                         )
              ),
              tabPanel("By Cancer", 
                       sidebarPanel(
                         selectizeInput("cancer_input", "Cancer", cancerNames, selected = 'colorectal cancer', multiple = FALSE, options = NULL),
                         plotlyOutput("cancer_piechart"),
                         htmlOutput("cancer_text")
                         #verbatimTextOutput("cancer_text")
                       ),
                       mainPanel(
                         plotlyOutput("cancer_barchart"),
                         DT::dataTableOutput("cancer_table")
                       )
              )
  )
  
)

#USPersonalExpenditure <- data.frame("Categorie"=rownames(USPersonalExpenditure), USPersonalExpenditure)
#data <- USPersonalExpenditure[,c('Categorie', 'X1960')]


testData <- data.frame(labels=c('a','b','c'),values=c(2,5,3))
plot_ly(testData, labels = ~labels, values = ~values, type = 'pie')

server <- function(input, output, session) {
  geneData <- reactive({
    table <- cancermineCounts[cancermineCounts$gene_standardized==input$gene_input,c('relationtype','gene_standardized','cancer_standardized','freq')]
    if (nrow(table)>0) {
      ordering <- aggregate(table$freq,by=list(table$cancer_standardized),FUN=sum)
      colnames(ordering) <- c('cancer_standardized','totalfreq')
      table <- dplyr::inner_join(table,ordering,by='cancer_standardized')
      table <- table[order(table$totalfreq,decreasing=TRUE),]
      table <- table[,c('relationtype','gene_standardized','cancer_standardized','freq')]
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$gene_table <- DT::renderDataTable({
    DT::datatable(geneData(),
                  selection = 'single',
                  rownames = FALSE,
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  geneTableProxy = dataTableProxy('gene_table')
  
  
  output$gene_piechart <- renderPlotly({
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
      
      plot_ly(completecounts, labels = ~relationtype, values = ~freq, type = 'pie', marker=list(colors=completecounts$color)) %>%
        layout(title = paste('Relation types for',input$gene_input),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      plotly_empty()%>% 
        config(displayModeBar = F)
    }
    
  })
  
  output$gene_barchart <- renderPlotly({
    table <- geneData()
    if (nrow(table) > 0) {
      unmelted <- dcast(geneData(), cancer_standardized~relationtype, sum, drop=F)
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$cancer_standardized <- factor(as.character(unmelted$cancer_standardized), unique(unmelted$cancer_standardized))
      
      plot_ly(unmelted, x=~cancer_standardized, y=~Tumor_Suppressor, source='gene_barchart', type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
        add_trace(y = ~Oncogene, name = 'Oncogene', marker=list(color=color_Oncogene)) %>%
        add_trace(y = ~Driver, name = 'Driver', marker=list(color=color_Driver))%>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
    } else {
      plotly_empty()%>% 
        config(displayModeBar = F)
    }
    
  })
  observe({
    
    table <- geneData()
    s <- event_data("plotly_click", source = "gene_barchart")
    relationtypeOptions <- c('Tumor_Suppressor','Oncogene','Driver')
    rowNumber <- rownames(table[table$cancer_standardized==s$x & table$relationtype==relationtypeOptions[s$curveNumber+1],])
    geneTableProxy %>% selectRows(as.numeric(rowNumber))
  })
  
  
  output$gene_text <- renderPrint({
    if(length(input$gene_table_rows_selected)==0) {
      cat("Select a row from the table to see the associated sentences and citations")
    } else {
      table <- geneData()
      row <- table[input$gene_table_rows_selected,]
      entries <- cancermine[cancermine$relationtype==row$relationtype & cancermine$cancer_standardized==row$cancer_standardized & cancermine$gene_standardized==row$gene_standardized,]
      
      cat(entries$preparedText,sep='<br /><br />')
    }
  })
  
  
  
  
  
  
  
  
  
  
  cancerData <- reactive({
    table <- cancermineCounts[cancermineCounts$cancer_standardized==input$cancer_input,c('relationtype','gene_standardized','cancer_standardized','freq')]
    if (nrow(table)>0) {
      ordering <- aggregate(table$freq,by=list(table$gene_standardized),FUN=sum)
      colnames(ordering) <- c('gene_standardized','totalfreq')
      table <- dplyr::inner_join(table,ordering,by='gene_standardized')
      table <- table[order(table$totalfreq,decreasing=TRUE),]
      table <- table[,c('relationtype','gene_standardized','cancer_standardized','freq')]
      rownames(table) <- 1:nrow(table)
    }
    table
  })
  
  output$cancer_table <- DT::renderDataTable({
    DT::datatable(cancerData(),
                  selection = 'single',
                  rownames = FALSE,
                  options = list(lengthMenu = c(5, 30, 50), pageLength = 20))
  })
  
  cancerTableProxy = dataTableProxy('cancer_table')
  
  
  output$cancer_piechart <- renderPlotly({
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
      
      plot_ly(completecounts, labels = ~relationtype, values = ~freq, type = 'pie', marker=list(colors=completecounts$color)) %>%
        layout(title = paste('Relation types for',input$cancer_input),
               xaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
               yaxis = list(showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE))%>% 
        config(displayModeBar = F)
    } else {
      plotly_empty()%>% 
        config(displayModeBar = F)
    }
    
  })
  
  output$cancer_barchart <- renderPlotly({
    table <- cancerData()
    if (nrow(table) > 0) {
      unmelted <- dcast(table, gene_standardized~relationtype, sum, drop=F)
      unmelted$total <- unmelted$Driver + unmelted$Oncogene + unmelted$Tumor_Suppressor
      unmelted <- unmelted[unmelted$total>0,]
      unmelted <- unmelted[order(unmelted$total,decreasing=T),]
      unmelted$gene_standardized <- factor(as.character(unmelted$gene_standardized), unique(unmelted$gene_standardized))
      
      plot_ly(unmelted, x=~gene_standardized, y=~Tumor_Suppressor, source='cancer_barchart', type = 'bar', name = 'Tumor Suppressor', marker=list(color=color_TumorSuppressor)) %>%
        add_trace(y = ~Oncogene, name = 'Oncogene', marker=list(color=color_Oncogene)) %>%
        add_trace(y = ~Driver, name = 'Driver', marker=list(color=color_Driver))%>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack', margin = list(b = 200), xaxis=list(title = "", tickangle = 45))%>% 
        config(displayModeBar = F)
    } else {
      plotly_empty()%>% 
        config(displayModeBar = F)
    }
    
  })
  
  observe({
    table <- cancerData()
    s <- event_data("plotly_click", source = "cancer_barchart")
    relationtypeOptions <- c('Tumor_Suppressor','Oncogene','Driver')
    rowNumber <- rownames(table[table$gene_standardized==s$x & table$relationtype==relationtypeOptions[s$curveNumber+1],])
    cancerTableProxy %>% selectRows(as.numeric(rowNumber))
  })
  
  
  output$cancer_text <- renderPrint({
    if(length(input$cancer_table_rows_selected)==0) {
      cat("Select a row from the table to see the associated sentences and citations")
    } else {
      table <- cancerData()
      row <- table[input$cancer_table_rows_selected,]
      entries <- cancermine[cancermine$relationtype==row$relationtype & cancermine$cancer_standardized==row$cancer_standardized & cancermine$gene_standardized==row$gene_standardized,]
      
      cat(entries$preparedText,sep='<br /><br />')
    }
  })
  
  
  
  
  
  
  
  
  
  
  
  
}

shinyApp(ui, server)
