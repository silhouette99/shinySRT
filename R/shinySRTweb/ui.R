library(shiny)
library(shinyhelper)
library(shinydashboard)
library(data.table)
library(Matrix)
library(DT)
library(scales)
library(ggiraph)
library(shinythemes)
library(bslib)
require(shinySRT)




ui <- navbarPage(
  theme = shinythemes::shinytheme('united'),
  title = 'shinySRT web',
  tabPanel(
    shinybusy::add_busy_spinner(spin = "cube-grid",
                                color = "red"),
    title = 'Analysis',
    # titlePanel((strong(HTML(
    #   "Analysis"
    # )))),
    # sidebarLayout(
    br(),
    br(),
    br(),
    br(),
    fluidRow(
      sidebarPanel(
        title = 'Input',
        width = 2,
        numericInput(
          inputId = 'resv',
          label = 'Resolusion',
          min = 0,
          max = 2,
          value = 1,
          step = 0.1
        ),
        numericInput(
          inputId = 'npcs',
          label = 'NPCs',
          min = 10,
          max = 50,
          value = 20
        ),
        numericInput(
          inputId = 'level_line',
          label = 'Unit level',
          min = 10,
          max = 50,
          value = 20
        ),
        br(),
        radioButtons(
          inputId = 'genemap',
          label = 'Gene mapping',
          choices = c(FALSE, TRUE)
        ),
        radioButtons(
          inputId = 'sp',
          label = 'species',
          choices = c('Human' = 'hg', 'Mouse' = 'mm')
        ),
        br(),
        actionButton(inputId = 'go', label = 'Go to shiny', icon("gear")),
        tags$hr(),
        br(),
        actionButton(inputId = 'res_web', label = 'Refresh', icon("refresh")),
        br()
      ),
      mainPanel(
        width = 10,
        HTML('<center><img src="shinysrt.png" width = 50%></center>'),
        tags$hr(),
        tabsetPanel(
          tabPanel(title = 'Seurat object',
                   shiny::fluidRow(
                     br(),
                     br(),
                     column(
                       width = 8,
                       fileInput(
                         inputId = 'obj_input1',
                         label = 'Choose Object1 Rds File',
                         accept = c(".Rds", '.rds')
                       ),
                       fileInput(
                         inputId = 'obj_input2',
                         label = 'Choose Object2 Rds File',
                         accept = c(".Rds", '.rds')
                       ),
                       br(),
                       actionButton('run_button_obj', 'Run Analysis', icon = icon('play'))
                       
                     ),
                     column(1),
                     column(
                       width = 3,
                       h4('This interface can generate shinySRT, requiring only the provision of relevant spatial transcriptomics data to produce a visual interface for spatial transcriptomics.', style = 'font-size:15px;color:#a2a1a8;font-family:arial'),
                       h4('Seurat object: The Seurat object should be stored in .rds format. Once uploaded, click "Run Analysis" to process the data. After click "Go to shiny" to display the standard interface of shinySRT.', style = 'font-size:15px;color:#a2a1a8;font-family:arial')
                     )
                   )),
          tabPanel(
            title = 'SpaceRanger output',
            br(),
            br(),
            fluidRow(
              column(
                width = 8,
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'mtx_10x_input1',
                      label = 'Choose matrix File1',
                      accept = c(".h5", '.mtx')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'gene_10x_input1',
                      label = 'Choose rownames File1',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'barcode_10x_input1',
                      label = 'Choose colnames File1',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  )
                ),
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'image_10x_input1',
                      label = 'Choose image File1',
                      accept = c(".png")
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'json_10x_input1',
                      label = 'Choose scalefactor File1',
                      accept = c(".json")
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'position_10x_input1',
                      label = 'Choose position File1',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  )
                ),
                br(),
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'mtx_10x_input2',
                      label = 'Choose matrix File2',
                      accept = c(".h5", '.mtx')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'gene_10x_input2',
                      label = 'Choose rownames File2',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'barcode_10x_input2',
                      label = 'Choose colnames File2',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  )
                ),
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'image_10x_input2',
                      label = 'Choose image File2',
                      accept = c(".png")
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'json_10x_input2',
                      label = 'Choose scalefactor File2',
                      accept = c(".json")
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'position_10x_input2',
                      label = 'Choose position File2',
                      accept = c(".csv", '.tsv')
                    ),
                    width = 4
                  )
                )
              ),
              column(width = 1),
              column(
                width = 3,
                h4('SpaceRanger output: Ulpoad raw data obtained from SpaceRanger. Sparse matrix need to upload the gene names (features.tsv) and barcode files (barcodes.tsv).', style = 'font-size:15px;color:#a2a1a8;font-family:arial')
              )
            ),
            br(),
            actionButton('run_button_10x', 'Run Analysis', icon = icon('play'))
          ),
          tabPanel(
            title = 'Data matrix',
            br(),
            br(),
            fluidRow(
              column(
                width = 8,
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'mtx_input1',
                      label = 'Choose matrix File1',
                      accept = c(".csv", '.tsv', '.txt', '.xlsx')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'meta_input1',
                      label = 'Choose meta File1',
                      accept = c(".csv", '.tsv', '.txt', '.xlsx')
                    ),
                    width = 4
                  ),
                  column(
                    radioButtons(
                      inputId = 'x_rev1',
                      label = 'X Reverse1',
                      choices = c(FALSE, TRUE),
                      selected = FALSE
                    ),
                    width = 4
                  )
                ),
                fluidRow(column(
                  fileInput(
                    inputId = 'image_input1',
                    label = 'Choose image File1',
                    accept = c(".png", '.jpeg', '.tiff')
                  ),
                  width = 4
                ),
                column(
                  fileInput(
                    inputId = 'position_input1',
                    label = 'Choose position File1',
                    accept = c(".csv", '.tsv', '.txt', '.xlsx')
                  ),
                  width = 4
                )),
                br(),
                fluidRow(
                  column(
                    fileInput(
                      inputId = 'mtx_input2',
                      label = 'Choose matrix File2',
                      accept = c(".csv", '.tsv', '.txt', '.xlsx')
                    ),
                    width = 4
                  ),
                  column(
                    fileInput(
                      inputId = 'meta_input2',
                      label = 'Choose meta File2',
                      accept = c(".csv", '.tsv', '.txt', '.xlsx')
                    ),
                    width = 4
                  ),
                  column(
                    radioButtons(
                      inputId = 'x_rev2',
                      label = 'X Reverse2',
                      choices = c(FALSE, TRUE),
                      selected = FALSE
                    ),
                    width = 4
                  )
                ),
                fluidRow(column(
                  fileInput(
                    inputId = 'image_input2',
                    label = 'Choose image File2',
                    accept = c(".png", '.jpeg', '.tiff')
                  ),
                  width = 4
                ),
                column(
                  fileInput(
                    inputId = 'position_input2',
                    label = 'Choose position File2',
                    accept = c(".csv", '.tsv', '.txt', '.xlsx')
                  ),
                  width = 4
                ))
              ),
              column(width = 1),
              column(
                width = 3,
                h4('Data matrix: upload the expression matrix of spatial transcriptomics data.', style = 'font-size:15px;color:#a2a1a8;font-family:arial')
              )
            ),
            br(),
            actionButton('run_button_mtx', 'Run Analysis', icon = icon('play'))
          )
        )
      )
    ),
    tags$hr(),
    textOutput(outputId = 'ex'),
    uiOutput(outputId = 'fullPage'),
    br()
    
    # )
  ),tabPanel(title = 'About',width = 9,
             HTML('<center><img src="shinysrt.png" width = 70%></center>'),
             h1(strong('shinySRT web operation')),
             br(),br(),
             h3('NPCs'),
             h4('The PCA dimensions used during the Seurat data processing steps of FindNeighbors, FindClusters, and RunUMAP.'),
             tags$hr(),br(),
             h3('Resolusion'),
             h4('The "resolution" parameter in FindClusters is used to adjust the granularity of clustering.'),
             tags$hr(),br(),
             h3('Gene mapping'),
             h4('When "Gene mapping = TRUE," it will convert between gene symbol IDs and gene Ensembl IDs.'),
             tags$hr(),br(),
             h3('species'),
             h4('The source species of the sample'),
             tags$hr(),br(),
             h3('Go to shiny'),
             h4('Open the interface generated by shinySRT.'),
             tags$hr(),br(),
             h3('Refresh'),
             h4('Refresh the interface. If you want to upload files to generate a new interface, please refresh first.'),
             br(),br())
  
)
