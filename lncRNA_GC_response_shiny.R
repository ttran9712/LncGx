library(shiny)
library(shinyWidgets)
library(ggplot2)
library(DT)
library(plyr)
library(dplyr)
library(htmltools)

load("dataset_input.RData")

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      radioButtons("cell_type", label = "Cell type",
                   choices = list("B cells","CD4+ T cells","Endothelial cells","Fibroblasts","Monocytes","Myoblasts","Neutrophils","Osteoblasts","Preadipocytes"),selected = F),
      checkboxGroupInput("change_dir","Direction of change",
                         choiceNames = c("Up-regulated","Down-regulated","NS"),
                         choiceValues = c("up","down","NS"),
                         selected = c("up","down","NS")),
      p("Drag the sliders to choose the desired adjusted p-value and log2 fold change cutoffs"),
      shinyWidgets::sliderTextInput("padj","Adjusted p-value cutoff", choices=c(0.00001,0.0001,0.001,0.01,0.1), selected=0.01, grid = T),
      sliderInput("log2fc","Log2 fold change cutoff (+/-)", min=1, max=10, value=1, step = 0.1),
      uiOutput("textInput"),
      br(),
      dataTableOutput("entry_info"),
      htmlOutput("error_message4"),
      htmlOutput("error_message2"),
      tags$head(tags$style(HTML("#error_message2{color: black;
                                          font-size: 14px;
                                          max-height: 145px;
                                          overflow: auto;}")))
    ),
    mainPanel(
      textOutput("error_message"),
      tags$head(tags$style("#error_message{color: red;
                                          font-size: 18px;}"
      )),
      uiOutput("plotUserGeneUI"),
      uiOutput("plot"),
      uiOutput("click_info_UI"),
      uiOutput("brush_info_UI"),
    )
  )
)


server <- function(input, output) {
  output$textInput <- renderUI({
    if (length(input$cell_type)==1 & length(input$change_dir)>0) {
      tagList(
        p("Use the cursor to click on or brush across dots on plot to retrieve information about selected genes. To highlight the response of genes of interest to GC, you can upload a text file (.txt) and/or enter the genes manually.", style = "font-size:11pt"),
        fileInput("file_gene","Upload txt file",accept='.txt'),
        textInput("gene","Enter the gene name"),
        actionButton("add_gene","Add",style="width:132px"),
        actionButton("delete_gene","Remove",style="width:132px"),
        actionButton("clear","Clear all",style="width:132px"),
        br(),
        br(),
        actionBttn("display_user_gene","Display my gene",block=TRUE,size="sm",style="simple",color="primary")
      )
    }
  })
  
  #remember log2fc and padj input when user wants to display gene
  slider <- reactiveValues()
  slider$log2fc_display <- c()
  slider$padj_display <- c()
  observeEvent(input$display_user_gene, {
    slider$log2fc_display <- input$log2fc
    slider$padj_display <- input$padj
  })
  
  #UI for displaying plot with gene of interest
  Display <- eventReactive(input$display_user_gene, {
    if (any(ChangeToVector() %in% rownames(datasetInput()))==FALSE) {
      tagList(
      htmlOutput("error_message3"),
      tags$head(tags$style(HTML("#error_message3{color: black;
                                          font-size: 14px;
                                          max-height: 145px;
                                          overflow: auto;}"))),
      br(),
      plotOutput("plot1",
                 click = "plot1_click",
                 brush = "plot1_brush")
      )
    }
    else {
    tagList(
      htmlOutput("error_message3"),
      tags$head(tags$style(HTML("#error_message3{color: black;
                                          font-size: 14px;
                                          max-height: 145px;
                                          overflow: auto;}"))),
      br(),
      plotOutput("plot_user_gene",
                 click = "plot1_click",
                 brush = "plot1_brush"),
      DTOutput("user_gene_info")
    )
    }
  })
  
    output$plotUserGeneUI <- renderUI({
    if (isTRUE(input$display_user_gene>0 && length(input$change_dir)==0)) {
      NULL
    }
    if (isTRUE(input$display_user_gene>0 && (input$log2fc!=slider$log2fc_display || input$padj!=slider$padj_display) && length(input$change_dir)>0)) {
      plotOutput("plot1",
                 click = "plot1_click",
                 brush = "plot1_brush")
    }
    else if (isTRUE(input$display_user_gene>0 && input$log2fc==slider$log2fc_display && input$padj==slider$padj_display && length(input$change_dir)>0)) {
      Display()
    } 
  })
  
  #UI for displaying base plot without queried gene displayed or displaying DT once plot is clicked or brushed
  output$plot <- renderUI({
    if(isTRUE(input$display_user_gene==0)) {
      plotOutput("plot1",
                 click = "plot1_click",
                 brush = "plot1_brush")
    }
  })
  
  output$click_info_UI <- renderUI({
    if(length(input$plot1_click)>0 & length(input$plot1_brush)==0) {
      DTOutput("click_info")
    }
  })
  
  output$brush_info_UI <- renderUI({
    if(length(input$plot1_brush)>0) {
      DTOutput("brush_info")
    }
  })

  #Store user's chosen dataset in memory
  datasetInput <- reactive({
    if (length(input$cell_type)==0){inputData = NULL}
    else if(input$cell_type=="B cells"){inputData = sign_filtered_Bcells_6H_combined_df}
    else if (input$cell_type=="CD4+ T cells"){inputData = sign_filtered_CD4_6H_combined_df}
    else if (input$cell_type=="Endothelial cells"){inputData = sign_filtered_Endo_6H_combined_df}
    else if (input$cell_type=="Fibroblasts"){inputData = sign_filtered_Fibr_6H_combined_df}
    else if (input$cell_type=="Monocytes"){inputData = sign_filtered_Mono_6H_combined_df}
    else if (input$cell_type=="Myoblasts"){inputData = sign_filtered_Myob_6H_combined_df}
    else if (input$cell_type=="Neutrophils"){inputData = sign_filtered_Neut_6H_combined_df}
    else if (input$cell_type=="Osteoblasts"){inputData = sign_filtered_Ostb_6H_combined_df}
    else if (input$cell_type=="Preadipocytes"){inputData = sign_filtered_PrAd_6H_combined_df}
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
  
  entry_info <- reactiveValues()
  entry_info$entry <- data.frame(stringsAsFactors = FALSE)
  entry_info$list <- c()
  
  # if enter genes manually, add the new gene to the list each time
  update_entry_info <- eventReactive(input$add_gene,{
    non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_6H_combined_df),rownames(sign_filtered_CD4_6H_combined_df),rownames(sign_filtered_Endo_6H_combined_df),rownames(sign_filtered_Fibr_6H_combined_df),rownames(sign_filtered_Mono_6H_combined_df),rownames(sign_filtered_Myob_6H_combined_df),rownames(sign_filtered_Neut_6H_combined_df),rownames(sign_filtered_Ostb_6H_combined_df),rownames(sign_filtered_PrAd_6H_combined_df)))
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
      output$error_message2 <- renderText({"<font color=\"#FF0000\">Please change the cell type and/or direction of change to display this gene.</font>"})
    }
    else if (gene%in%rownames(datasetInput())==FALSE && !(gene%in%non_na_lncrna_gene_names)) {
      output$error_message2 <- renderText({"<font color=\"#FF0000\"><b>WARNING</b> Gene symbol not recognized. Please enter a valid human gene symbol.</font>"})
    }
  }
  )
  # delete gene
  remove_gene <- eventReactive(input$delete_gene, {
    gene <- toupper(input$gene)
    row_num <- which(entry_info$entry==gene)
    row_num_selected <- input$entry_info_rows_selected
    if (!is.null(row_num_selected)){
      output$error_message2 <- renderText("")
      entry_info$entry <- data.frame(entry_info$entry[-row_num_selected,],stringsAsFactors = FALSE)
    }
    else {
      if (length(row_num)>0){
        output$error_message2 <- renderText("")
        entry_info$entry <- data.frame(entry_info$entry[-row_num,],stringsAsFactors = FALSE)
      }
      if (length(row_num)==0 & gene!="") {
        output$error_message2 <- renderText("<font color=\"#FF0000\"><b>WARNING</b> Gene symbol not recognized. Please enter a valid human gene symbol.</font>")
        entry_info$entry <- entry_info$entry
      }
    }
    })
  
  # clear all genes from list
  clear_gene_list <- eventReactive(input$clear, {
    entry_info$entry <- data.frame(entry_info$entry[-c(1:nrow(entry_info$entry)),],stringsAsFactors = FALSE)
    output$error_message2 <- renderText("")
  })
  
  # if upload a txt file, present it in the form of DataTable
  LoadFile <- eventReactive(input$file_gene, {
    inFile <- input$file_gene
    invalid_entry_index = c()
    invalid_entry_string = c()
    invalid_entry_message = ""
    genes_to_add = data.frame(Gene=NA,stringsAsFactors = FALSE)
    non_na_lncrna_gene_names <- unique(c(rownames(sign_filtered_Bcells_6H_combined_df),rownames(sign_filtered_CD4_6H_combined_df),rownames(sign_filtered_Endo_6H_combined_df),rownames(sign_filtered_Fibr_6H_combined_df),rownames(sign_filtered_Mono_6H_combined_df),rownames(sign_filtered_Myob_6H_combined_df),rownames(sign_filtered_Neut_6H_combined_df),rownames(sign_filtered_Ostb_6H_combined_df),rownames(sign_filtered_PrAd_6H_combined_df)))
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
  
  #base plot
  plotInput <- reactive({
    if (length(input$cell_type)==0) {
      plot <- NULL
    }
    else if (length(input$cell_type)==1) {
      if (length(input$change_dir)==0) {
        plot <- NULL
      }
      else {
        if (length(input$change_dir %in% c("up","down","NS"))==3) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                              ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                                     ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))),size=1.1) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +  
            scale_color_manual(values=c("red","blue","light grey","grey29"),
                               breaks=c("up","down","NS_1","NS_2")) +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,"up")) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color="up")) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values="red",
                               breaks="up") +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,"down") && input$log2fc>=1) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color="down")) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values="blue",
                               breaks="down") +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,"NS")) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2"))) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values=c("light grey","grey29"),
                               breaks=c("NS_1","NS_2")) +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,c("up","down"))) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up","down"))) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values=c("red","blue"),
                               breaks=c("up","down")) +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,c("up","NS"))) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                              ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values=c("red","light grey","grey29"),
                               breaks=c("up","NS_1","NS_2")) +
            scale_x_continuous(name = NULL, limits = c(0,chr_sizes$cumsum[24]), expand = c(0.005,0), breaks = xaxis_lab, labels = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY")) +
            scale_y_continuous(limits = c(-270,270), breaks = seq(-250,250,50)) +
            theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
                  axis.title.y = element_text(size = 14),
                  axis.text.y=element_text(size=13),
                  axis.title.x = element_text(size = 14),
                  axis.text.x=element_text(size=13,angle=45,vjust=0.55),
                  legend.position = "none")
        }
        if (identical(input$change_dir,c("down","NS"))) {
          plot <- ggplot(datasetInput()) +
            geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                              ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))) +
            ylab("GC response \n signed -log10(padj)") +
            ggtitle(input$cell_type) +
            theme_classic() +
            scale_color_manual(values=c("blue","light grey","grey29"),
                               breaks=c("down","NS_1","NS_2")) +
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
    }
    return(plot)
  })
  
  #plot displaying queried gene(s)
  plotUserGene <- eventReactive(input$display_user_gene, {
    highlight <- datasetInput() %>%
      filter(rownames(datasetInput()) %in% ChangeToVector())
    plot <- ggplot(datasetInput()[!(rownames(datasetInput()) %in% ChangeToVector()),]) +
      geom_point(aes(x=plot_coord, y=signed_negLog10_padj, color=ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange > input$log2fc,"up",
                                                                        ifelse(signed_negLog10_padj*sign(log2FoldChange) >= -log10(input$padj) & log2FoldChange < -input$log2fc,"down",
                                                                               ifelse((signed_negLog10_padj*sign(log2FoldChange) < -log10(input$padj) | (log2FoldChange > -input$log2fc & log2FoldChange < input$log2fc)) & chr_names=="chr1" | chr_names=="chr3" | chr_names=="chr5" | chr_names=="chr7" | chr_names=="chr9" | chr_names=="chr11" | chr_names=="chr13" | chr_names=="chr15" | chr_names=="chr17" | chr_names=="chr19" | chr_names=="chr21" | chr_names=="chrX","NS_1","NS_2")))),size=1.1) +
      scale_color_manual(values=c("red","blue","light grey","grey29"),
                         breaks=c("up","down","NS_1","NS_2")) +
      geom_point(data=highlight, aes(x=plot_coord,y=signed_negLog10_padj),color="black",fill="gold1",size=3,shape=23,stroke=2) +
      ylab("GC response \n signed -log10(padj)") +
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
    return(plot)
  })
  
  #output
  userGeneInfo <- eventReactive(input$display_user_gene, {
    DT <- datatable(
      data = datasetInput()[rownames(datasetInput()) %in% ChangeToVector(),],
      caption = htmltools::tags$caption(htmltools::tags$b("Information about highlighted gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE)
      )
  })
  
  output$plot1 <- renderPlot({plotInput()})
  
  output$plot_user_gene <- renderPlot({plotUserGene()})
  
  output$user_gene_info <- renderDT({userGeneInfo()})
  
  clickInfo <- eventReactive(input$plot1_click, {
    datatable(
      data = nearPoints(datasetInput(), input$plot1_click, addDist = F),
      caption = htmltools::tags$caption(htmltools::tags$b("Information about clicked gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE)
    )
  })
  
  output$click_info <- renderDT({clickInfo()})
  
  brushInfo <- eventReactive(input$plot1_brush, {
    datatable(
      data = brushedPoints(datasetInput(), input$plot1_brush),
      caption = htmltools::tags$caption(htmltools::tags$b("Information about brushed gene(s)"), style = "color:purple; caption-side: top; text-align:left;", rownames=FALSE)
    )
  })
  
  output$brush_info <- renderDT({brushInfo()})  
  
  #Error messages
  observe({
    if (length(input$cell_type)==0) {
      output$error_message <- renderText({"Please select a cell type."})
    }
    else if (length(input$cell_type)==1) {
      if (length(input$change_dir)==0) {
        output$error_message <- renderText({"Please select the change in direction(s)."})
      }
      else {
        output$error_message <- renderText("")
      }
    }
  })
  
  observeEvent (input$display_user_gene, {
    if (all(ChangeToVector() %in% rownames(datasetInput()))) {
      output$error_message3 <- renderText("")
    }
    else if (all(ChangeToVector() %in% rownames(datasetInput()))==FALSE) {
      user_gene <- ChangeToVector()[!(ChangeToVector() %in% rownames(datasetInput()))]
      user_gene[2:(length(user_gene)+1)] <- user_gene[1:length(user_gene)]
      user_gene[1] <- "<font color=\"#FF0000\"><b>WARNING</b> Please change the cell type and/or direction of change to display the following genes:</font>"
      output$error_message3 <- renderUI({HTML(paste(user_gene,"<br/>"))})
    }
  })
}


shinyApp(ui, server)