library(shiny)
library(shinyjs)
library(ggplot2)
library(ggrepel)
library(DT)
library(plyr)
library(dplyr)
library(htmltools)

load("dataset_input.RData")

server <- function(input, output, session) {
  
#UI for main panel
  output[["box"]] <- renderUI({
    if (length(input$cell_type)==1 & length(input$change_dir)>0) {
      box(width = NULL,
          title = tags$span(style = "font-size: 17px; font-weight: bold;","OPTIONAL: find your gene(s) of interest"),
          status = "primary",
          solidHeader = TRUE,
      tagList(
        p(tags$span(style = "font-size: 15px;","To highlight the response of your lncRNA gene(s) of interest to glucocorticoids, you can upload a text file (.txt) and/or enter the gene(s) manually.")),
        fileInput("file_gene","Upload .txt file",accept='.txt'),
        textInput("gene","Enter the gene name"),
        actionButton("add_gene","Add",style="width:140px"),
        actionButton("delete_gene","Remove",style="width:140px"),
        br(),
        actionButton("clear","Clear all",style="width:140px"),
        br(),
        br(),
        actionBttn("display_user_gene","Search",block=TRUE,size="sm",style="simple",color="primary"),
        br(),
        dataTableOutput("entry_info"),
        htmlOutput("error_message4"),
        htmlOutput("error_message2"),
        tags$head(tags$style(HTML("#error_message2{color: black;
                                          font-size: 14px;
                                          max-height: 145px;
                                          overflow: auto;}")))
      )
      )
    }
  })

  entry_info <- reactiveValues()
  entry_info$entry <- data.frame(stringsAsFactors = FALSE)
  entry_info$list <- c()

# handle "add", "delete", "clear all", "file upload" events
# if enter genes manually, add the new gene to the list each time
  update_entry_info <- eventReactive(input$add_gene,{
    if (input$timepoint=="6-hour") {non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_6H_combined_df),rownames(sign_filtered_CD4_6H_combined_df),rownames(sign_filtered_Endo_6H_combined_df),rownames(sign_filtered_Fibr_6H_combined_df),rownames(sign_filtered_Mono_6H_combined_df),rownames(sign_filtered_Myob_6H_combined_df),rownames(sign_filtered_Neut_6H_combined_df),rownames(sign_filtered_Ostb_6H_combined_df),rownames(sign_filtered_PrAd_6H_combined_df)))}
    else {non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_2H_combined_df),rownames(sign_filtered_CD4_2H_combined_df),rownames(sign_filtered_Endo_2H_combined_df),rownames(sign_filtered_Fibr_2H_combined_df),rownames(sign_filtered_Mono_2H_combined_df),rownames(sign_filtered_Myob_2H_combined_df),rownames(sign_filtered_Neut_2H_combined_df),rownames(sign_filtered_Ostb_2H_combined_df),rownames(sign_filtered_PrAd_2H_combined_df)))}
    
    gene <- toupper(input$gene)
    output$error_message4 <- renderText("")
    if (gene%in%rownames(datasetInput())==TRUE) {
      if (length(entry_info$entry)==0) {
        output$error_message2 <- renderText("")
        entry_info$entry[nrow(entry_info$entry)+1,1] <- gene
      }
      else if (length(entry_info$entry)>0 && gene%in%entry_info$entry[,1]==FALSE){ # a new gene, no replicate in the list
        output$error_message2 <- renderText("")
        entry_info$entry[nrow(entry_info$entry)+1,1] <- gene
        print(entry_info$entry)
      }
      else if (length(entry_info$entry)>0 && gene%in%entry_info$entry[,1]==TRUE){ # gene has already been added to the list
        output$error_message2 <- renderText("<font color=\"#FF0000\">The gene is already in the list.</font>")
        entry_info$entry <- entry_info$entry
      }
    }
    else if (gene%in%rownames(datasetInput())==FALSE && gene%in%non_na_lncrna_gene_names) {
      output$error_message2 <- renderText({"<font color=\"#FF0000\">Please change the cell type, timepoint, and/or direction of change to display this gene.</font>"})
    }
    else if (gene%in%rownames(datasetInput())==FALSE && !(gene%in%non_na_lncrna_gene_names)) {
      output$error_message2 <- renderText({"<font color=\"#FF0000\"><b>WARNING</b> LncRNA gene not recognized. Please enter a valid human lncRNA gene.</font>"})
    }
  }
  )
  
# delete gene
  remove_gene <- eventReactive(input$delete_gene, {
    print(nrow(entry_info$entry))
    gene <- toupper(input$gene)
    row_num <- which(entry_info$entry==gene)
    row_num_selected <- input$entry_info_rows_selected
    if (!is.null(row_num_selected)){
      output$error_message2 <- renderText("")
      entry_info$entry <- data.frame(entry_info$entry[-row_num_selected,],stringsAsFactors = FALSE)
    }
    if (nrow(entry_info$entry)==0) {
      entry_info$entry <- data.frame()
    }
  })
  
# clear all genes from list
  clear_gene_list <- eventReactive(input$clear, {
    entry_info$entry <- data.frame()
    output$error_message2 <- renderText("")
  })
  
# if upload a .txt file, present it in the form of DataTable
  LoadFile <- eventReactive(input$file_gene, {
    inFile <- input$file_gene
    invalid_entry_index = c()
    invalid_entry_string = c()
    invalid_entry_message = ""
    genes_to_add = data.frame(Gene=NA,stringsAsFactors = FALSE)
    
    if (input$timepoint=="6-hour") {non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_6H_combined_df),rownames(sign_filtered_CD4_6H_combined_df),rownames(sign_filtered_Endo_6H_combined_df),rownames(sign_filtered_Fibr_6H_combined_df),rownames(sign_filtered_Mono_6H_combined_df),rownames(sign_filtered_Myob_6H_combined_df),rownames(sign_filtered_Neut_6H_combined_df),rownames(sign_filtered_Ostb_6H_combined_df),rownames(sign_filtered_PrAd_6H_combined_df)))}
    else {non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_2H_combined_df),rownames(sign_filtered_CD4_2H_combined_df),rownames(sign_filtered_Endo_2H_combined_df),rownames(sign_filtered_Fibr_2H_combined_df),rownames(sign_filtered_Mono_2H_combined_df),rownames(sign_filtered_Myob_2H_combined_df),rownames(sign_filtered_Neut_2H_combined_df),rownames(sign_filtered_Ostb_2H_combined_df),rownames(sign_filtered_PrAd_2H_combined_df)))}
    
    if (is.null(inFile)) entry_info$entry <- entry_info$entry
    else{ 
      genes_from_file <- read.table(inFile$datapath,stringsAsFactors = FALSE) # store the genes from the file temporarily for check before adding to the list
      for (i in 1:nrow(genes_from_file)) genes_from_file[i,] <- toupper(genes_from_file[i,]) # make sure all genes are in upper case
      for (i in 1:nrow(genes_from_file)){ # check if there is any invalid genes in the uploaded file
        if (genes_from_file[i,1]%in%non_na_lncrna_gene_names==FALSE){
          invalid_entry_index = append(invalid_entry_index,i)
          invalid_entry_string = append(invalid_entry_string,as.character(genes_from_file[i,1]))
        }   
      }
      if (is.null(invalid_entry_index)==FALSE){ # when an invalid gene is detected
        for (i in length(invalid_entry_index):1) genes_from_file <- data.frame(genes_from_file[-invalid_entry_index[i],],stringsAsFactors = FALSE) # remove invalid genes
        
        if (length(invalid_entry_string)>0){ # if there is an invalid gene, present the error message
          invalid_gene <- invalid_entry_string
          invalid_gene[2:(length(invalid_gene)+1)] <- invalid_gene[1:length(invalid_gene)]
          invalid_gene[1] <- "<font color=\"#FF0000\"><b>WARNING</b> The following genes cannot be found:</font>"
          output$error_message2 <- renderUI({HTML(paste(invalid_gene,"<br/>"))})
        }
        else # if all genes are valid, no error message displayed
          output$error_message2 <- renderText("")
      }
      
      for (i in nrow(genes_from_file):1){ # remove all duplicated genes in the first loaded file
        if (genes_from_file[i,1]%in%genes_from_file[1:i-1,1]) { 
          genes_from_file <- data.frame(genes_from_file[-i,],stringsAsFactors = FALSE)
          print(str(genes_from_file)) }
      }
      
      if (length(entry_info$entry)>0){
        if (all(genes_from_file[,1] %in% entry_info$entry[,1])==TRUE) {
          genes_to_add <- genes_from_file
        }
        else if (all(genes_from_file[,1] %in% entry_info$entry[,1])==FALSE) {
          for (i in nrow(genes_from_file):1){ # if any new gene is a duplicate of existed gene, do not add it
            if (genes_from_file[i,1] %in% entry_info$entry[,1]==TRUE){
              output$error_message4 <- renderText("<font color=\"#FF0000\"><b>WARNING</b> Duplicated gene(s) removed</font>")
              genes_from_file <- data.frame(genes_from_file[-i,])
            }
          }
          for (i in 1:nrow(genes_from_file)){ # make the list of genes with NO duplicate
            genes_to_add[i,1] <- as.character(genes_from_file[i,1])
          }
        }
      }
      else if (length(entry_info$entry)==0) genes_to_add <- genes_from_file
      
      for (i in 1:nrow(genes_from_file)){
        if (length(entry_info$entry)>0 && as.character(genes_from_file[i,1])%in%entry_info$entry[,1]==FALSE){
          names(entry_info$entry) <- "Genes"
          genes_to_add <- data.frame(genes_to_add)
          names(genes_to_add) <- "Genes"
          entry_info$entry <- rbind(entry_info$entry,genes_to_add,stringsAsFactors=FALSE)
        }
        if (length(entry_info$entry)==0){
          entry_info$entry <- rbind.fill(entry_info$entry,genes_from_file)
        }
      }
    }
    print(entry_info$entry)
  })
  
# output the gene list in form of table
  output$entry_info <- renderDataTable({
    observeEvent(input$file_gene,LoadFile())
    observeEvent(input$add_gene,update_entry_info())
    observeEvent(input$delete_gene,remove_gene())
    observeEvent(input$clear,clear_gene_list())
    entry_info$entry
    print(entry_info$entry)
  }, colnames = "Genes", selection = 'multiple'
  )
  
# change the DataTable into a vector
  ChangeToVector <- function(){
    entry_info$list <- c()
    for (i in 1:nrow(entry_info$entry)){
      entry_info$list[i] <- as.character(entry_info$entry[i,1])
    }
    return(Reduce(c,entry_info$list))
  }
  
#Filter data based on user's parameters and return filtered data
  datasetInput <- reactive({
    if (length(input$cell_type)==0){inputData = NULL}
    else {
      if (input$timepoint=="6-hour") {
        if(input$cell_type=="B cells"){inputData = sign_filtered_Bcells_6H_combined_df}
        else if (input$cell_type=="CD4+ T cells"){inputData = sign_filtered_CD4_6H_combined_df}
        else if (input$cell_type=="Endothelial cells"){inputData = sign_filtered_Endo_6H_combined_df}
        else if (input$cell_type=="Fibroblasts"){inputData = sign_filtered_Fibr_6H_combined_df}
        else if (input$cell_type=="Monocytes"){inputData = sign_filtered_Mono_6H_combined_df}
        else if (input$cell_type=="Myoblasts"){inputData = sign_filtered_Myob_6H_combined_df}
        else if (input$cell_type=="Neutrophils"){inputData = sign_filtered_Neut_6H_combined_df}
        else if (input$cell_type=="Osteoblasts"){inputData = sign_filtered_Ostb_6H_combined_df}
        else if (input$cell_type=="Preadipocytes"){inputData = sign_filtered_PrAd_6H_combined_df}
      }
      else {
        if(input$cell_type=="B cells"){inputData = sign_filtered_Bcells_2H_combined_df}
        else if (input$cell_type=="CD4+ T cells"){inputData = sign_filtered_CD4_2H_combined_df}
        else if (input$cell_type=="Endothelial cells"){inputData = sign_filtered_Endo_2H_combined_df}
        else if (input$cell_type=="Fibroblasts"){inputData = sign_filtered_Fibr_2H_combined_df}
        else if (input$cell_type=="Monocytes"){inputData = sign_filtered_Mono_2H_combined_df}
        else if (input$cell_type=="Myoblasts"){inputData = sign_filtered_Myob_2H_combined_df}
        else if (input$cell_type=="Neutrophils"){inputData = sign_filtered_Neut_2H_combined_df}
        else if (input$cell_type=="Osteoblasts"){inputData = sign_filtered_Ostb_2H_combined_df}
        else if (input$cell_type=="Preadipocytes"){inputData = sign_filtered_PrAd_2H_combined_df}
      }
    }
    dataset <- data.frame()
    
    userData <- function(data) {
      if (length(input$cell_type)==0) {
        dataset <- NULL
      }
      else if (length(input$cell_type)==1) {
        if (length(input$change_dir)==0) {
          dataset <- NULL
        }
        if (length(input$change_dir %in% c("up","down","NS"))==3) {
          dataset <- data
        }
        if (identical(input$change_dir,"up")) {
          dataset <- subset(data, signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc)
        }
        if (identical(input$change_dir,"down")) {
          dataset <- subset(data, signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc)
        }
        if (identical(input$change_dir,"NS")) {
          dataset <- subset(data, (signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj)) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc))
        }
        if (identical(input$change_dir,c("up","down"))) {
          dataset <- subset(data, signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & (log2FoldChange > input$log2fc | log2FoldChange < -input$log2fc))
        }
        if (identical(input$change_dir,c("up","NS"))) {
          dataset <- subset(data, (signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc) | ((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj)) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)))
        }
        if (identical(input$change_dir,c("down","NS"))) {
          dataset <- subset(data, (signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc) | ((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj)) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)))
        }
      }
      return(dataset)
    }
    
    userData(data = inputData)
  })
  
#generate base manhattan plot based on filtered data (filtered based on user's parameters)
  p <- reactiveValues(data = NULL)
  
  observe({
    if (length(input$cell_type)==0) {
      plot <- NULL
    }
    else if (length(input$cell_type)==1) {
      if (length(input$change_dir)==0) {
        plot <- NULL
      }
      else {
        plot <- ggplot(datasetInput()) +
          geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                            ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                                  ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))),size=1.1) +
          scale_color_manual(values=c("red","blue","light grey","grey29"),
                            breaks=c("up","down","NS_1","NS_2")) +
          ylab(expression(atop(bold("GC response"), "signed -log10(padj)"))) +
          ggtitle(input$cell_type) +
          theme_classic() +  
          scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
          scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
          theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                axis.title.y = element_text(size = 14),
                axis.text.y=element_text(size=13),
                axis.title.x = element_text(size = 14),
                axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                legend.position = "none")
      }
    }
    p$data <- plot
  })

#generate manhattan plot with user's gene(s) of interest highlighted
  observeEvent(input$display_user_gene, {
    if (any(ChangeToVector() %in% rownames(datasetInput()))==FALSE) { #if user's gene(s) of interest are not found in filtered data, return base manhattan plot
      plot <- ggplot(datasetInput()) +
        geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                          ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                                 ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))),size=1.1) +
        scale_color_manual(values=c("red","blue","light grey","grey29"),
                           breaks=c("up","down","NS_1","NS_2")) +
        ylab(expression(atop(bold("GC response"), "signed -log10(padj)"))) +
        ggtitle(input$cell_type) +
        theme_classic() +  
        scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
        scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
        theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
              axis.title.y = element_text(size = 14),
              axis.text.y=element_text(size=13),
              axis.title.x = element_text(size = 14),
              axis.text.x=element_text(size=13,angle=45,vjust=0.55),
              legend.position = "none")
    }
    else {
      highlight <- datasetInput() %>%
        filter(rownames(datasetInput()) %in% ChangeToVector())
      plot <- ggplot(datasetInput()[!(rownames(datasetInput()) %in% ChangeToVector()),]) +
        geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                          ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                                 ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))),size=1.1) +
        scale_color_manual(values=c("red","blue","light grey","grey29"),
                           breaks=c("up","down","NS_1","NS_2")) +
        geom_point(data=highlight, aes(x=plot_coord,y=signed_negLog10_padj),color="black",fill="gold1",size=3,shape=23,stroke=2) +
        geom_label_repel(data=highlight, aes(x=plot_coord,y=signed_negLog10_padj,label=rownames(highlight)), direction = 'y', vjust=ifelse(highlight$signed_negLog10_padj>=0,1,0), point.padding = 1, box.padding = 0.5, force = 5, max.overlaps = Inf) +
        ylab(expression(atop(bold("GC response"), "signed -log10(padj)"))) +
        ggtitle(input$cell_type) +
        theme_classic() +  
        scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
        scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
        theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
              axis.title.y = element_text(size = 14),
              axis.text.y=element_text(size=13),
              axis.title.x = element_text(size = 14),
              axis.text.x=element_text(size=13,angle=45,vjust=0.55),
              legend.position = "none")
    }
    p$data <- plot
  })

#render plot
  output$plot<- renderPlot({
    if(is.null(p$data)) return()
    p$data
  })

#control when to hide or show plot from user
  observe({
    if(length(input$cell_type)==0 | length(input$change_dir)==0) { #hide plot (and disable all plot interactions) if user does not select a cell type or change in direction
      shinyjs::hide("plot")
    }
    else {
      shinyjs::show("plot")
    }
  })

#create output in the form of a data table, displaying information about the clicked/brushed genes or highlighted genes on the plot
  userGeneInfo <- eventReactive(input$display_user_gene, {
    DT <- datatable(
      data = datasetInput()[rownames(datasetInput()) %in% ChangeToVector(),],
      colnames = c('chromosome names','start','end','width','strand','log2 fold change', 'Wald statistic', 'p value','adjusted<br>p value','signed<br>-log10(padj)',''),
      options = list(scrollX=TRUE,
                     autoWidth=TRUE, 
                     columnDefs = list(list(targets=c(6,8:10), visible=TRUE, width='100'), 
                                       list(targets=11, visible=FALSE), #hide last column with coordinates to plot manhattan plots
                                       list(targets='_all', className='dt-center'))),
      caption = htmltools::tags$caption(htmltools::tags$b("Information about highlighted gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE),
      escape = FALSE
    )
  })
  
  output$user_gene_info <- renderDT(userGeneInfo())
  
  clickInfo <- eventReactive(input$plot1_click, {
    datatable(
      data = nearPoints(datasetInput(), input$plot1_click, addDist = F),
      colnames = c('chromosome names','start','end','width','strand','log2 fold change', 'Wald statistic', 'p value','adjusted<br>p value','signed<br>-log10(padj)',''),
      options = list(scrollX=TRUE,
                     autoWidth=TRUE, 
                     columnDefs = list(list(targets=c(6,8:10), visible=TRUE, width='100'),
                                       list(targets=11, visible=FALSE),
                                       list(targets='_all', className='dt-center'))),
      caption = htmltools::tags$caption(htmltools::tags$b("Information about clicked gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE),
      escape = FALSE
    )
  })
  
  output$click_info <- renderDT(clickInfo())
  
  brushInfo <- eventReactive(input$plot1_brush, {
    datatable(
      data = brushedPoints(datasetInput(), input$plot1_brush),
      colnames = c('chromosome names','start','end','width','strand','log2 fold change', 'Wald statistic', 'p value','adjusted<br>p value','signed<br>-log10(padj)',''),
      options = list(scrollX=TRUE,
                     autoWidth=TRUE, 
                     columnDefs = list(list(targets=c(6,8:10), visible=TRUE, width='100'),
                                       list(targets=11, visible=FALSE),
                                       list(targets='_all', className='dt-center'))),
      caption = htmltools::tags$caption(htmltools::tags$b("Information about brushed gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE),
      escape = FALSE
    )
  })
  
  output$brush_info <- renderDT(brushInfo()) 
  
#control datable output of clicked/brush event or when user wants to display their gene(s) of interest
  i <- reactiveValues(interact = NULL, user_gene_table = NULL) # i$interact = reactive variable controlling datable output of brush/clicked events
                                                               # i$user_gene_table = reactive variable controlling datable output when user wants to query and display their gene of interest (i.e. click "search")
  
  observeEvent(input$plot1_click, {
    i$interact <-  div(DT::DTOutput("click_info"), style="width:100%")
  })
    
  observe({
    if (length(input$cell_type)==1) { #if the user clicks on a different cell type, i$user_gene_table is set to NULL
      i$user_gene_table <- NULL
    }
  })
  
  observeEvent(input$display_user_gene, {
    i$interact <- NULL #if the user clicks "search", i$interact is set to NULL (so that whatever datatable resulting from a clicked/brush event that is currently displayed will be gone)
    if (nrow(entry_info$entry)==0 | any(ChangeToVector() %in% rownames(datasetInput()))==FALSE) { #if user's gene(s) of interest are not found in filtered data or if no gene is in query list, set i$user_gene_table to NULL
      i$user_gene_table <- NULL
    }
    else {
      i$user_gene_table <- div(DT::DTOutput("user_gene_info"), style="width:100%")
    }
  })
  
  observeEvent(input$plot1_brush, {
    i$interact <- div(DT::DTOutput("brush_info"), style="width:100%")
  })
  
  output$interaction <- renderUI({
    if(is.null(i$user_gene_table)) { 
      i$interact
    }
    else {
      tagList(i$user_gene_table,
              i$interact)
    }
  })
  
# clear brush after anything that triggers an action is clicked
  makeReactiveBinding("plot1_brush")
  observe({
    if (length(input$cell_type)==1 | length(input$timepoint)==1 | length(input$change_dir)>0 | isTRUE(input$padj) | isTRUE(input$log2fc) | isTRUE(input$display_user_gene)) {
    session$resetBrush("plot1_brush")
    }
  })
  
#Handle any additional error messages
  observe({
    if (length(input$cell_type)==0) {
      output$error_message <- renderText({"To begin, please select the cell type in which you want to view the glucocorticoid response."})
    }
    else if (length(input$cell_type)==1) {
      if (length(input$change_dir)==0) {
        output$error_message <- renderText({"Please select the change in direction(s)."})
      }
      else {
        output$error_message <- renderText("Use the cursor to click on or brush across dots on plot to retrieve information about selected genes.")
        output$error_message2 <- renderText("")
        output$error_message3 <- renderUI({})
        output$error_message4 <- renderText("")
      }
    }
  })
  
  observeEvent (input$display_user_gene, {
    output$error_message2 <- renderText("")
    output$error_message4 <- renderText("")
    if (all(ChangeToVector() %in% rownames(datasetInput()))) {
      output$error_message3 <- renderText("")
    }
    else if (all(ChangeToVector() %in% rownames(datasetInput()))==FALSE) {
      user_gene <- ChangeToVector()[!(ChangeToVector() %in% rownames(datasetInput()))]
      user_gene[2:(length(user_gene)+1)] <- user_gene[1:length(user_gene)]
      user_gene[1] <- "<font color=\"#FF0000\"><b>WARNING</b> Please change the cell type, timepoint, and/or direction of change to display the following genes:</font>"
      if (is.na(user_gene[2])) {
        output$error_message3 <- renderUI({})
      }
      else {
        output$error_message3 <- renderUI({HTML(paste(user_gene,"<br/>"))})
      }
    }
  })
}
