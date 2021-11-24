library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)

ui <- fluidPage(
        dashboardPage(
          dashboardHeader(title = "LncGx", titleWidth = 300),
            dashboardSidebar(disable = TRUE),
            dashboardBody(
              tags$head(tags$style(HTML('.same-row {max-width: 200px; display: table-cell; padding-right: 10px;}'))),
              fluidRow(
                column(3, 
                       box(width = NULL,
                           title = tags$span(style = "font-size: 17px; font-weight: bold;","Generate your plot:"),
                           status = "primary",
                           solidHeader = TRUE,
                           radioButtons("cell_type", tags$span(style = "font-size: 16px;","Cell type"),
                                        choiceNames = list(
                                          tags$span(style = "font-size: 15px;","B cells"),
                                          tags$span(style = "font-size: 15px;","CD4+ T cells"),
                                          tags$span(style = "font-size: 15px;","Endothelial cells"),
                                          tags$span(style = "font-size: 15px;","Fibroblasts"),
                                          tags$span(style = "font-size: 15px;","Monocytes"),
                                          tags$span(style = "font-size: 15px;","Myoblasts"),
                                          tags$span(style = "font-size: 15px;","Neutrophils"),
                                          tags$span(style = "font-size: 15px;","Osteoblasts"),
                                          tags$span(style = "font-size: 15px;","Preadipocytes")),
                                        choiceValues = c("B cells", "CD4+ T cells", "Endothelial cells", "Fibroblasts", "Monocytes", "Myoblasts","Neutrophils", "Osteoblasts", "Preadipocytes"),
                                        selected = F),
                            radioButtons("timepoint", tags$span(style = "font-size: 16px;","Timepoint"),
                                        choiceNames = list(
                                          tags$span(style = "font-size: 15px;","2-hour"),
                                          tags$span(style = "font-size: 15px;","6-hour")),
                                        choiceValues = c("2-hour", "6-hour"),
                                        selected = "6-hour"),
                            checkboxGroupInput("change_dir",tags$span(style = "font-size: 16px;","Direction of change"),
                                        choiceNames = list(
                                          tags$span(style = "font-size: 15px;","Upregulation"),
                                          tags$span(style = "font-size: 15px;","Downregulation"),
                                          tags$span(style = "font-size: 15px;","Not significant")),
                                        choiceValues = c("up","down","NS"),
                                        selected = c("up","down","NS")),
                            p(tags$span(style = "font-size: 15px;","Drag the sliders to choose the desired adjusted p-value and log2 fold change cutoffs")),
                            shinyWidgets::sliderTextInput("padj",tags$span(style = "font-size: 16px;","Adjusted p-value cutoff"), choices=c(0.00001,0.0001,0.001,0.01,0.1), selected=0.01, grid = F),
                            sliderInput("log2fc",tags$span(style = "font-size: 16px;","Log2 fold change cutoff (+/-)"), min=1, max=10, value=1, step = 0.1, ticks = FALSE)
                            ),
                       uiOutput("box")
                       ),
                column(9,
                       box(width = NULL,
                           status = "primary",
                           solidHeader = TRUE,
                           useShinyjs(),
                           textOutput("error_message"),
                           br(),
                           tags$head(tags$style("#error_message{font-size: 17px;}")),
                           plotOutput("plot",
                                      width = "100%",
                                      height = "400px",
                                      click = "plot1_click",
                                      brush = "plot1_brush"),
                           uiOutput("interaction"),
                           uiOutput("error_message3")
                           )
                        )
                        ))))
