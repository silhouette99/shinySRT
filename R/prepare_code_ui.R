#' ui code prepared for shinyspatial
#'
#'
#'
#' @param df_select default information to provide the slice number
#'
#'
#' ui_load : load package in ui scirpts
#' ui_head : with web's header and menu
#'
#' @return data files required for shiny app
#'
#' @export
## loaded library
ui_load <- function() {
  glue::glue(
    'library(shiny)\n',
    'library(shinyhelper)\n',
    'library(shinydashboard)\n',
    'library(data.table)\n',
    'library(Matrix)\n',
    'library(DT)\n',
    'library(scales)\n',
    'meta <- readRDS(\'meta.Rds\')\n',
    'meta_group <- readRDS(\'meta_group.Rds\')\n',
    'genesets <- readRDS(\'genesets.Rds\')\n',
    'df_select <- readRDS(\'df_select.Rds\')\n\n\n'
  )
}
## ui header
ui_head <- function(title){
  glue::glue(
    '### Start server code\n',
    'ui <- dashboardPage(\n',
    '  ### project title\n',
    '  dashboardHeader(title = \"{title}\"),\n',
    '  ### side meta menu\n',
    '  dashboardSidebar(\n',
    '    ### page menu\n',
    '    sidebarMenu(\n',
    '      menuItem(\"SpotInfo vs GeneExpr\",tabName = \"dashboard\",icon = icon(\"dashboard\")),\n',
    '      menuItem(\"GeneExpr vs GeneExpr\",tabName = \"d_gene_ex\",icon = icon(\"sun\")),\n',
    '      menuItem(\"Viol / Box data chart\",tabName = \"vio\",icon = icon(\"desktop\")),\n',
    '      menuItem(\"Portion chart\", tabName = \"portion\", icon = icon(\"bars\")),\n',
    '      menuItem(\"Heat / Dot plot\", tabName = \"heatdot\", icon = icon(\"gear\")),\n',
    '      menuItem(\"Gene coexpression\",tabName = \"coexp\",icon = icon(\"code\"))),\n',
    '    ### import new spot information\n',
    '    shiny::selectizeInput(inputId = \"new_title\",label = \'group name\',choices = NULL,multiple = FALSE,options = list(create = TRUE)),\n',
    '    textAreaInput(\"impspifo\",HTML(\'Import new spots information\'),height = \"200px\",value = NULL),\n',
    '    actionButton(inputId = \'swi_spinfo\',label = \'Import\',icon = icon(\'house\')),\n',
    '    sliderInput(inputId = \'point_size\',label = \'spots size\',min = 0,ticks = FALSE,max = 1,value = 1)),\n',
    '  ### function page\n',
    '  dashboardBody(\n',
    '    tabItems(\n',
    '      tabItem(tabName = \"dashboard\",\n',
    '              h2(strong(HTML(\"SpotInfo vs GeneExpr\"))),\n',
    '              h4(\"Spot information vs gene expression on reduced dimensions\"),\n',
    '              \"In this tab, users can visualise both Spot information and gene \",\n',
    '              \"expression side-by-side on slice spatial representions.\",\n',
    '              br(),br(),\n\n'
  )
}
## ui page1
page1_ui_head <- function(){
  glue::glue(
    'fluidRow(column(12,style = \"border-bottom: 2px solid black\",\n',
    '     ### dlclick key\n',
    '     keys::useKeys(), keys::keysInput(inputId = \"keys\", keys = c(\"d\", \"e\")), br(),br(),\n',
    '     ### first row:the control panel\n',
    '     fluidRow(\n',
    '       ### spot information selected (1th row 1th column)\n',
    '       column(h4(strong(\"Spot information\")),width = 4,br(),\n',
    '              shiny::selectizeInput(inputId = \"meta1_select\",label = \'select group\',choices = meta_group$group,multiple = FALSE,\n',
    '                                    selected = df_select$meta1,options = list(create = TRUE),width = \"50%\")),\n',
    '       ### gene selected (1th row 2th column)\n',
    '       column(width = 4,h4(strong(\"Gene expression\")),br(),\n',
    '              shiny::selectInput(inputId = \"gene\",label = \"gene names:\",choices = NULL,width = \"50%\"),\n',
    '              shiny::selectInput(inputId = \"plotout\",label = \"highlight\",selected = FALSE,width = \'50%\',choices = c(FALSE, TRUE))),\n',
    '       ### add spots information (1th row 3th column)\n',
    '       column(width = 4, h4(strong(\"Select spot\")), br(),\n',
    '              fluidRow(\n',
    '                column(width = 4,\n',
    '                       fluidRow(\n',
    '                         ### new meta column names\n',
    '                         shiny::selectizeInput(inputId = \"group1\",label = \'Group name\',choices = NULL,multiple = FALSE,options = list(create = TRUE))),\n',
    '                       fluidRow(\n',
    '                         ### clusters names\n',
    '                         shiny::selectizeInput(inputId = \"cluster1\",label = \'annotation clusters\',choices = NULL,multiple = FALSE,options = list(create = TRUE))\n',
    '                       ),br(),\n',
    '                       actionButton(inputId = \'circle_select\',label = \'action\',icon = icon(\"bar-chart-o\"))),\n',
    '                column(width = 1),\n',
    '                column(width = 6,br(),\n',
    '                       fluidRow(\n',
    '                         ### save new infomation\n',
    '                         column(width = 6,\n',
    '                                shiny::actionButton(inputId = \"Import_1\",icon = icon(\"th\"),label = \"Import the meta\")),\n',
    '                         ### reset\n',
    '                         column(width = 6,\n',
    '                                shiny::actionButton(inputId = \"reset1\",icon = icon(\"refresh\"),label = \"reset all\"))),\n',
    '                       ### click function\n',
    '                       shiny::selectInput(inputId = \'p_tool\',label = \'Pen or Rubber\',choices = c(\'Pen\', \'Rubber\'),selected = \'Pen\',width = \'40%\'),\n',
    '                       br(),\n',
    '                       ### click area size\n',
    '                       fluidRow(\n',
    '                         sliderInput(inputId = \'select_size\',label = \'selected area\',min = 10,ticks = FALSE,\n',
    '                                     max = df_select$maxarea,\n',
    '                                     value = df_select$maxarea / 5)))))),br(),\n',
    '         ### second row:spatial plot\n',
    '         fluidRow(\n\n'
  )
}
ui_p1_p <- function(df_select){
  p1_code <- lapply(1:length(unique(df_select$slice)), function(x){
    if(x == length(unique(df_select$slice))){
      ends = ''
    }else{
      ends = ','
    }

    glue::glue(
      'fluidRow(column(4,style = \"border-right: 2px solid black\",\n',
      '         box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'group\'),status = \"primary\",solidHeader = TRUE,\n',
      '             shiny::plotOutput(outputId = \'bg_plot{x}\',click = \'bg_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '                    hover = hoverOpts(id = \"bg_polyg_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = df_select$image_size)),\n',
      '         fluidRow(column(6,br(),downloadButton(\"bg_spot{x}.pdf\", \"Download PDF\"),downloadButton(\"bg_spot{x}.png\", \"Download PNG\")),\n',
      '           column(6,div(style = \"display:inline-block\",numericInput(\"bg_spot{x}.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
      '                  div(style = \"display:inline-block\",numericInput( \"bg_spot{x}.w\",\"PDF / PNG width:\",width = \"138px\",min = 4,max = 20,value = 12,step = 0.5))))),\n',
      '  column(4,style = \"border-right: 2px solid black\",\n',
      '         box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'feature\'),status = \"primary\",solidHeader = TRUE,\n',
      '             shiny::plotOutput(outputId = \'feature{x}\',click = \'f_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '                    hover = hoverOpts(id = \"f_polyg_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = df_select$image_size)),\n',
      '         fluidRow(column(6,br(),downloadButton(\"ft_spot{x}.pdf\", \"Download PDF\"),downloadButton(\"ft_spot{x}.png\", \"Download PNG\")),\n',
      '           column(6,div(style = \"display:inline-block\",numericInput(\"ft_spot{x}.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
      '                  div(style = \"display:inline-block\",numericInput(\"ft_spot{x}.w\",\"PDF / PNG width:\",width = \"138px\",min = 4,max = 20,value = 12,step = 0.5))))),\n',
      '  column(4,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'highlight\'),status = \"primary\",solidHeader = TRUE,\n',
      '      container(width = 12,shiny::div(class = \"large-plot\",\n',
      '          shiny::plotOutput(outputId = \'highlight{x}\',click = \'h_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '            hover = hoverOpts(id = \"h_polyg_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = df_select$image_size),\n',
      '          shiny::plotOutput(outputId = \'highlight{x}_sm\',click = \'h_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '            hover = hoverOpts(id = \"h_polyg_{x}\",delay = 100,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = df_select$image_size),\n',
      '          shiny::tags$style(\"\\n                        .large-plot {{\\n                            position: relative;\\n                        }}\\n                        #highlight{x} {{\n                            position: absolute;\\n                        }}\\n                        #highlight{x}_sm {{\\n                            position: absolute;\\n                        }\\n\\n                      \")\n',
      '        ))))){ends}\n\n'
    )

  })
  return(p1_code)
}
p1_tab_ui <- function(){
  glue::glue(
    '),br(),br(),),\n',
    'fluidRow(column(width = 7,align = \"center\",br(),br(),box(status = \"primary\",solidHeader = TRUE,collapsible = TRUE,DTOutput(outputId = \'spotsinfor\'),width = 10)),\n',
    '  column(5,br(),br(),h4(\"spot numbers / statistics\"),radioButtons(\"splt\",\"Split continuous spot info into:\",choices = c(\"Quartile\", \"Decile\"),\n',
    '          selected = \"Decile\",inline = TRUE),shiny::sliderInput(inputId = \'stat_thrd\',label = \"Gene expression threshold\",min = 0,max = 1,value = 0),\n',
    '         box(width = 12,align = \"center\",status = \"primary\",solidHeader = TRUE,collapsible = TRUE,DTOutput(outputId = \'statis_gene\')))))),\n\n'
  )
}
## ui page2
page2_ui_head <- function(){
  glue::glue(
    'tabItem(tabName = \"d_gene_ex\",h2(strong(HTML(\"GeneExpr vs GeneExpr\"))),\n',
    '   h4(\'Gene expression vs gene expression on spatial represention\'),br(),br(),\n',
    '     fluidRow(column(12,style = \"border-bottom: 2px solid black\",br(),br(),\n',
    '       fluidRow(column(width = 4,h4(strong(\"Select spot\")),br(),fluidRow(column(12,\n',
    '                      fluidRow(column( 4,shiny::selectizeInput(inputId = \"meta2_select\",label = \'select group\',choices = meta_group$group,\n',
    '                                                               multiple = FALSE,selected = df_select$meta1,options = list(create = TRUE))),\n',
    '                               column(4,shiny::selectizeInput(inputId = \"group2\",label = \'Group name\',choices = NULL,multiple = FALSE,options = list(create = TRUE))),\n',
    '                               column(4,shiny::selectizeInput(inputId = \"cluster2\",label = \'annotation clusters\',choices = NULL,multiple = FALSE,options = list(create = TRUE)))),\n',
    '                      fluidRow(column(2,actionButton(inputId = \'circle_select2\',label = \'action\',icon = icon(\"bar-chart-o\"))),\n',
    '                               column(2,shiny::actionButton(inputId = \"reset2\",icon = icon(\"refresh\"),label = \"reset all\")),\n',
    '                               column(2,shiny::actionButton(inputId = \"Import_2\",icon = icon(\"th\"),label = \"Import\")),column(1),\n',
    '                               column(5,shiny::selectInput(inputId = \'p_tool_2\',label = \'Pen or Rubber\',choices = c(\'Pen\', \'Rubber\'),selected = \'Pen\'))),\n',
    '                      fluidRow(column(12,sliderInput(inputId = \'select_size2\',label = \'selected area\',min = 10,ticks = FALSE,max = df_select$maxarea *10,value = df_select$maxarea / 5)))))),\n',
    '             column(width = 4,h4(strong(\"Gene expression\")),br(),\n',
    '                    shiny::selectInput(inputId = \"gene_f1\",label = \"gene names:\",choices = NULL,width = \"50%\"),\n',
    '                    shiny::selectInput(inputId = \"plotout2\",label = \"highlight\",selected = FALSE,width = \'50%\',choices = c(FALSE, TRUE))),\n',
    '             column(width = 4,h4(strong(\"Gene expression\")),br(),\n',
    '                    shiny::selectInput(inputId = \"gene_f2\",label = \"gene names:\",choices = NULL,width = \"50%\"))),br(),\n',
    '                fluidRow(\n'
  )
}
ui_p2_p <- function(df_select) {
  p2_code <- lapply(1:length(unique(df_select$slice)), function(x) {
    if (x == length(unique(df_select$slice))) {
      ends = '),br(),br()))),'
    } else{
      ends = ','
    }

    glue::glue(
      'fluidRow(column(4,style = \"border-right: 2px solid black\",\n',
      ' box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'spot selected\'),status = \"primary\",solidHeader = TRUE,\n',
      '  container(width = 12,\n',
      '    shiny::div(class = \"large-plot\",shiny::plotOutput(outputId = \'bg_plot_f{x}\',click = \'bg_f_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '      hover = hoverOpts(id = \"bg_f_polyg_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = 600),\n',
      '        shiny::plotOutput(outputId = \'plot_f{x}_sm\',click = \'bg_f_click_{x}\',dblclick = \'h_dbl_click1\',\n',
      '      hover = hoverOpts(id = \"bg_f_polyg_{x}\",delay = 100,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = 600),\n',
      '       shiny::tags$style(\"\\n                        .large-plot {{\\n                            position: relative;\\n                        }}\\n                        #bg_plot_f{x} {{\\n                            position: absolute;\\n                        }}\\n                        #plot_f{x}_sm {{\\n                            position: absolute;\\n                        }}\\n\\n                      \"))))),\n',
      '  column(4,style = \"border-right: 2px solid black\",\n',
      '         box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'gene 1\'),status = \"primary\",solidHeader = TRUE,\n',
      '             shiny::plotOutput(outputId = \'gene_feature1_{x}\',click = \'gf_click_1_{x}\',dblclick = \'h_dbl_click1\',\n',
      '                               hover = hoverOpts(id = \"gf_polyg_1_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = 600)),\n',
      '         fluidRow(column(6,br(),downloadButton(\"ft_g1_spot{x}.pdf\", \"Download PDF\"),downloadButton(\"ft_g1_spot{x}.png\", \"Download PNG\")),\n',
      '                  column(6,div(style = \"display:inline-block\",numericInput(\"ft_g1_spot{x}.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
      '                         div(style = \"display:inline-block\",numericInput(\"ft_g1_spot{x}.w\",\"PDF / PNG width:\",width = \"138px\",min = 4,max = 20,value = 12,step = 0.5))))),\n',
      '  column(4,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{x}], \'gene 2\'),status = \"primary\",solidHeader = TRUE,\n',
      '               shiny::plotOutput(outputId = \'gene_feature2_{x}\',click = \'gf_click_2_{x}\',dblclick = \'h_dbl_click1\',\n',
      '                                 hover = hoverOpts(id = \"gf_polyg_2_{x}\",delay = 10,delayType = \"throttle\",clip = TRUE,nullOutside = TRUE),height = 600)),\n',
      '         fluidRow(column(6,br(),downloadButton(\"ft_g2_spot{x}.pdf\", \"Download PDF\"),downloadButton(\"ft_g2_spot{x}.png\", \"Download PNG\")),\n',
      '                  column(6,div(style = \"display:inline-block\",numericInput(\"ft_g2_spot{x}.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
      '                         div(style = \"display:inline-block\",numericInput(\"ft_g2_spot{x}.w\",\"PDF / PNG width:\",width = \"138px\",min = 4,max = 20,value = 12,step = 0.5\n',
      '                         )))))){ends}\n\n'
    )
  })
  return(p2_code)
}
## ui page3
page3_ui_head <- function(){
  glue::glue(
    'tabItem(tabName = \"vio\",h2(strong(HTML(\"Viol / Box data chart\"))),h4(\'Spot information / gene expression violin plot / box plot\'),\n',
    '   \'In this tab, users can visualise the gene expression or continuous spot information (e.g. Number of UMIs / module score) across groups of spots (e.g. libary / clusters).\',\n',
    '   br(),br(),\n',
    '   fluidRow(column(12,fluidRow(column(3,style = \"border-right: 2px solid black\",\n',
    '     box(width = 12,br(),br(),\n',
    '         fluidRow(shiny::selectizeInput(inputId = \"vio_x\",label = \'Spot information - x\',choices = meta_group[info == T, group],\n',
    '                                        selected = df_select$meta1,multiple = FALSE,options = list(create = TRUE))),br(),\n',
    '         fluidRow(shiny::selectInput(choices = NULL,multiple = FALSE,inputId = \'vio_y\',label = \'statistical info - y\')),br(),\n',
    '         fluidRow(shiny::selectizeInput(choices = c(\'No\', meta_group[info == T, group]),\n',
    '                                        selected = \'No\',inputId = \'vio_gb\',label = \'group by\',options = list(create = TRUE))),\n',
    '         shiny::radioButtons(\"plot_type\",\"Plot type:\",choices = c(\"violin\", \"boxplot\"),selected = \"violin\",inline = TRUE),\n',
    '         checkboxInput(\"point_ex\", \"Show data points\", value = FALSE),br(),\n',
    '         fluidRow(shiny::selectizeInput(inputId = \"subspot\",label = \'subset spot\',choices = meta_group$group,\n',
    '                                        selected = \'seurat_clusters\',multiple = F,options = list(create = TRUE))),br(),\n',
    '         fluidRow(column(12,actionButton(\"Tog1\", \"Toggle to subset spot\"),\n',
    '                         conditionalPanel(condition = \"input.Tog1 % 2 == 1\",\n',
    '                                          shiny::selectInput(inputId = \"meta_tog2\",label = \'select spot information to subset:\',\n',
    '                                                             choices = meta_group[info == T, group],selected = df_select$meta1),\n',
    '                                          uiOutput(\"sc1a1sub1.ui\"),\n',
    '                                          actionButton(\"sc1a1sub1all\", \"Select all groups\", class = \"btn btn-primary\"),\n',
    '                                          actionButton(\"sc1a1sub1non\", \"Deselect all groups\", class = \"btn btn-primary\")))),br(),br(),\n',
    '         fluidRow(h4(strong(\'Vlnplot/Boxplot\')),\n',
    '                  column(12,fluidRow(\n',
    '                    column(6,br(),downloadButton(\"vb_spot.pdf\", \"Download PDF\"),downloadButton(\"vb_spot.png\", \"Download PNG\")),\n',
    '                    column(6,div(style = \"display:inline-block\",\n',
    '                                 numericInput(\"vb_spot.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
    '                                                            div(style = \"display:inline-block\",numericInput(\"vb_spot.w\",\"PDF / PNG width:\",width = \"138px\",\n',
    '                                                                                                                 min = 4,max = 20,value = 12,step = 0.5)))))))),\n\n',
  )
}
p3_vg_plot <- function(df_select){
  i = 1
  all_vbg <- glue::glue('column(9,\n')
  while(i <= length(df_select$slice)){
    if(i%%2 > 0){
      if(i == length(df_select$slice)){
        vbgp <- glue::glue('fluidRow(\n',
                           '  column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                           '    status = "primary",solidHeader = TRUE,shiny::plotOutput(outputId = \'vb_bg_plot{i}\',click = \'vb_click_{i}\',height = df_select$image_size)))),\n\n')
      }else{
        vbgp <- glue::glue('fluidRow(\n',
                           '  column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                           '    status = "primary",solidHeader = TRUE,shiny::plotOutput(outputId = \'vb_bg_plot{i}\',click = \'vb_click_{i}\',height = df_select$image_size))),\n\n')
      }
      # vbgp <- paste(vbgp,collapse = '')
    }else{
      vbgp <- glue::glue('  column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                         '    status = "primary",solidHeader = TRUE,shiny::plotOutput(outputId = \'vb_bg_plot{i}\',click = \'vb_click_{i}\',height = df_select$image_size)))),\n\n')
      # vbgp <- paste(vbgp,collapse = '')
    }
    all_vbg <- paste(all_vbg,vbgp,sep = '')
    i = i+1
  }
  vb <- c(
    '  fluidRow(column(12,\n',
    '    box(width = 12,height = 880,status = \"primary\",solidHeader = TRUE,\n',
    '      shiny::plotOutput(outputId = \'vb_plot\',height = 800))))))))),\n'
  )
  vb <- paste(vb,collapse = '')
  all_vbg <- paste(all_vbg,vb,sep = '')
  return(glue::glue('{all_vbg}'))
}
## page4
page4_ui <- function(df_select){
  page4_ui_head <- c(
    'tabItem(tabName = \"portion\",h2(strong(HTML(\"Portion data chart\"))),h4(\'Proportion / spot numbers across different spot information\'),\n',
    '        \'In this tab, users can visualise the composition of spots based on one discrete spots information across another discrete spots information.\',\n',
    '        br(),br(),fluidRow(column(12,fluidRow(\n',
    '          column(3,style = \"border-right: 2px solid black\",\n',
    '                 box(width = 12,br(),br(),fluidRow(\n',
    '                   shiny::selectizeInput(inputId = \"por_x\",label = \'Spot information - x\',choices = meta_group[info == T, group],\n',
    '                                         selected = df_select$meta1,multiple = FALSE,options = list(create = TRUE))),br(),\n',
    '                   fluidRow(shiny::selectizeInput(inputId = \"por_group\",label = \'color by information - group\',\n',
    '                                                  choices = meta_group[info == T, group],selected = df_select$meta2,multiple = FALSE,options = list(create = TRUE))),\n',
    '                   br(),shiny::radioButtons(\"tab_type\",\"Plot type:\",choices = c(\"percentage\", \"number\"),selected = \"percentage\",inline = TRUE),\n',
    '                   checkboxInput(\"flipxy\", \"Flip x y\", value = FALSE),br(),\n',
    '                   fluidRow(shiny::selectizeInput(inputId = \"subspot2\",label = \'subset spot\',choices = meta_group[info == T, group],\n',
    '                                                  selected = df_select$meta1,multiple = F,options = list(create = TRUE))),br(),\n',
    '                   fluidRow(column(12,actionButton(\"Tog2\", \"Toggle to subset spot\"),\n',
    '                       conditionalPanel(condition = \"input.Tog2 % 2 == 1\",\n',
    '                                        shiny::selectInput(inputId = \"meta_tog_por2\",label = \'select spot information to subset:\',\n',
    '                                                           choices = meta_group[info == T, group],selected = df_select$meta1),\n',
    '                                        uiOutput(\"spotsub1.ui\"),\n',
    '                                        actionButton(\"spotsuball\", \"Select all groups\", class = \"btn btn-primary\"),\n',
    '                                        actionButton(\"spotsubnon\", \"Deselect all groups\", class = \"btn btn-primary\")))),\n',
    '                   fluidRow(h4(strong(\'Barplot\')),column(12,\n',
    '                              fluidRow(column(6,br(),downloadButton(\"bar_spot.pdf\", \"Download PDF\"),downloadButton(\"bar_spot.png\", \"Download PNG\")),\n',
    '                                       column(6,\n',
    '                                              div(style = \"display:inline-block\",numericInput(\"bar_spot.h\",\"PDF / PNG height:\",width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
    '                                              div(style = \"display:inline-block\",numericInput(\"bar_spot.w\",\"PDF / PNG width:\",width = \"138px\",min = 4,max = 20,value = 12,step = 0.5)))))))),\n\n'
  )

  page4_ui_head <- paste(page4_ui_head,collapse = '')

  i = 1
  all_pbg <- glue::glue('column(9,\n')
  while(i <= length(df_select$slice)){
    if(i%%2 > 0){
      if(i == length(df_select$slice)){
        pogp <- glue::glue(
          '          fluidRow(\n',
          '            column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
          '               status = \"primary\",solidHeader = TRUE,shiny::plotOutput(outputId = \'por_bg_plot{i}\',click = \'por_click_{i}\',height = df_select$image_size)))),\n\n'
        )
      }else{
        pogp <- glue::glue(
          '          fluidRow(\n',
          '            column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
          '               status = \"primary\",solidHeader = TRUE,shiny::plotOutput(outputId = \'por_bg_plot{i}\',click = \'por_click_{i}\',height = df_select$image_size))),\n\n'
        )
      }
      # pogp <- paste(pogp,collapse = '')
    }else{
      pogp <- glue::glue('            column(6,box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                         '               status = \"primary\",solidHeader = TRUE,shiny::plotOutput(outputId = \'por_bg_plot{i}\',click = \'por_click_{i}\',height = df_select$image_size)))),\n\n')
      # pogp <- paste(pogp,collapse = '')
    }
    all_pbg <- paste(all_pbg,pogp,sep = '')
    i = i+1
  }

  po <- glue::glue('            fluidRow(column(12,box(width = 12,height = 880,status = \"primary\",solidHeader = TRUE,shiny::plotOutput(outputId = \'por_plot\', height = 800))))))))),\n\n')


  all_pbg <- paste(page4_ui_head,all_pbg,po,sep = '')

  return(glue::glue('{all_pbg}'))
}
## page5
page5_ui <- function(df_select){
  page5_ui_head <- glue::glue('tabItem(tabName = \"heatdot\",h2(strong(HTML(\"Gene expression heatmap/dotplot\"))),\n',
                              '  h4(\'In this tab, users can visualise the gene expression patterns of multiple genes grouped by categorical spot information (e.g. library / cluster).\'),\n',
                              '  \'The normalised expression are averaged, log-transformed and then plotted.\',br(),br(),\n',
                              '  fluidRow(column(12,fluidRow(column(3,style = \"border-right: 2px solid black\",box(width = 12,br(),br(),\n',
                              '        fluidRow(\n',
                              '          textAreaInput(\"inpgenelist\",\n',
                              '            HTML(\"List of gene names <br />\n',
                              '                          (Max 50 genes, separated <br />\n',
                              '                           by , or ; or newline):\"),height = \"200px\",\n',
                              '            value = paste0(df_select$genes, collapse = \", \")) %>%\n',
                              '            helper(type = \"inline\",size = \"m\",fade = TRUE,\n',
                              '              title = \"List of genes to plot on bubbleplot / heatmap\",\n',
                              '              content = c(\"Input genes to plot\",\"- Maximum 50 genes (due to ploting space limitations)\",\"- Genes should be separated by comma, semicolon or newline\"))),br(),\n',
                              '        fluidRow(shiny::selectizeInput(inputId = \"expgroup\",label = \'Group by:\',choices = meta_group[info == T, group],\n',
                              '            selected = df_select$meta1,multiple = FALSE,options = list(create = TRUE))),br(),\n',
                              '        shiny::radioButtons(\"exp_plot_type\",\"Plot type:\",choices = c(\"dotplot\", \"heatmap\"),selected = \"dotplot\",inline = TRUE),\n',
                              '        checkboxInput(\"scl_exp\", \"Scale gene expression\", value = TRUE),\n',
                              '        checkboxInput(\"rclust\", \"Cluster rows (gene)\", value = TRUE),\n',
                              '        checkboxInput(\"cclust\", \"Cluster columns (sample)\", value = FALSE),br(),\n',
                              '        fluidRow(shiny::selectizeInput(\n',
                              '            inputId = \"subspot3\",label = \'subset spot\',choices = meta_group[info == T, group],\n',
                              '            selected = df_select$meta1,multiple = F,options = list(create = TRUE))),br(),\n',
                              '        fluidRow(column(\n',
                              '          12,actionButton(\"Tog3\", \"Toggle to subset spot\"),\n',
                              '          conditionalPanel(condition = \"input.Tog3 % 2 == 1\",\n',
                              '            shiny::selectInput(inputId = \"meta_tog3\",label = \'select spot information to subset:\',\n',
                              '              choices = meta_group[info == T, group],selected = df_select$meta1),\n',
                              '            uiOutput(\"geneexpsub.ui\"),\n',
                              '            actionButton(\"geneexpsuball\", \"Select all groups\", class = \"btn btn-primary\"),\n',
                              '            actionButton(\"geneexpsubnon\", \"Deselect all groups\", class = \"btn btn-primary\")))),\n',
                              '        fluidRow(h4(strong(\'Dotplot / Heatmap\')),column(12,\n',
                              '          fluidRow(column(6,br(),\n',
                              '              downloadButton(\"heat_spot.pdf\", \"Download PDF\"),\n',
                              '              downloadButton(\"heat_spot.png\", \"Download PNG\")),\n',
                              '            column(6,div(style = \"display:inline-block\",\n',
                              '                numericInput(\"heat_spot.h\",\"PDF / PNG height:\",\n',
                              '                  width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
                              '              div(style = \"display:inline-block\",\n',
                              '                numericInput(\"heat_spot.w\",\"PDF / PNG width:\",\n',
                              '                  width = \"138px\",min = 4,max = 20,value = 12,step = 0.5)))))))),\n\n')

  i = 1
  all_dhbg <- glue::glue('column(9,\n')
  while(i <= length(df_select$slice)){
    if(i%%2 > 0){
      if(i == length(df_select$slice)){
        dhgp <- glue::glue('                     fluidRow(column(6,\n',
                           '                      box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                           '                        status = \"primary\",solidHeader = TRUE,\n',
                           '                        shiny::plotOutput(outputId = \'exp_bg_plot{i}\',click = \'exp_click_{i}\',height = df_select$image_size)))),\n\n'
        )
      }else{
        dhgp <- glue::glue(
          '                    fluidRow(column(6,\n',
          '                      box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
          '                        status = \"primary\",solidHeader = TRUE,\n',
          '                        shiny::plotOutput(outputId = \'exp_bg_plot{i}\',click = \'exp_click_{i}\',height = df_select$image_size))),\n\n'
        )
      }
      # dhgp <- paste(dhgp,collapse = '')
    }else{
      dhgp <- glue::glue('                    column(6,\n',
                         '                      box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'spot selected\'),\n',
                         '                        status = \"primary\",solidHeader = TRUE,\n',
                         '                        shiny::plotOutput(outputId = \'exp_bg_plot{i}\',click = \'exp_click_{i}\',height = df_select$image_size)))),\n\n')
      # dhgp <- paste(dhgp,collapse = '')
    }
    all_dhbg <- paste(all_dhbg,dhgp,sep = '')
    i = i+1
  }

  hd <- glue::glue('                    fluidRow(column(12,box(width = 12,height = 880,status = \"primary\",solidHeader = TRUE,shiny::plotOutput(outputId = \'exp_plot\', height = 800))))))))),')

  all_dhbg <- paste(page5_ui_head,all_dhbg,hd,sep = '')

  return(all_dhbg)
}
##page 6
page6_ui <- function(df_select){
  page6_ui_head <- glue::glue(
    'tabItem(tabName = \"coexp\",h2(strong(HTML(\"Gene coexpression\"))),h4(\'Coexpression of two genes\'),\n',
    '  \'In this tab, users can visualise the coexpression of two genes\',br(),br(),\n',
    '  fluidRow(column(12,fluidRow(column(width = 4,h4(strong(\"Select spot\")),br(),\n',
    '    fluidRow(column(12,\n',
    '      fluidRow(column(4,shiny::selectizeInput(\n',
    '          inputId = \"meta3_select\",label = \'select group\',\n',
    '          choices = meta_group[info == T, group],selected = df_select$meta1,\n',
    '          multiple = FALSE,options = list(create = TRUE)),\n',
    '        shiny::radioButtons(\"col_sel\",\"color select:\",\n',
    '          choices = c(\"Red (Gene1); Blue (Gene2)\",\"Orange (Gene1); Blue (Gene2)\",\n',
    '            \"Red (Gene1); Green (Gene2)\",\"Green (Gene1); Blue (Gene2)\"),\n',
    '          selected = \"Red (Gene1); Blue (Gene2)\",inline = TRUE)))))),\n',
    '    column(4,h4(strong(\"Gene expresion\")),br(),\n',
    '      shiny::selectInput(inputId = \"coexgene1\",label = \"gene1 names:\",choices = NULL,width = \"50%\"),\n',
    '      shiny::selectInput(inputId = \"coexgene2\",label = \"gene2 names:\",choices = NULL,width = \"50%\")),\n',
    '    column(4,h4(strong(\'Expression ratio\')),\n',
    '          fluidRow(column(6,\n',
    '              shiny::sliderInput(inputId = \'thrd1\',\n',
    '                label = \"Gene1 expression threshold\",\n',
    '                min = 0, max = 1,value = 0),\n',
    '              shiny::sliderInput(inputId = \'thrd2\',\n',
    '                label = \"Gene2 expression threshold\",\n',
    '                min = 0,max = 1,value = 0)),\n',
    '            column(6,\n',
    '              shiny::selectizeInput(inputId = \"group_tab\",\n',
    '                label = \'select group\',choices = c(\'No\', meta_group[info == T, group]),\n',
    '                multiple = FALSE,selected = \'No\',\n',
    '                options = list(create = TRUE)),\n',
    '              actionButton(inputId = \'reset_coexp\',label = \'reset all\',icon = icon(\"refresh\")))))),\n',
    '    fluidRow(column(8,\n\n'
  )

  all_cops <- c()
  for(i in 1:length(df_select$slice)){
    cop <- glue::glue(
      '      fluidRow(column(6,style = \"border-right: 2px solid black\",\n',
      '          box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'group\'),\n',
      '            status = \"primary\",solidHeader = TRUE,\n',
      '            shiny::plotOutput(outputId = \'coexp_bg_plot{i}\',click = \'coexp_click_{i}\',height = df_select$image_size))),\n',
      '        column(6,style = \"border-right: 2px solid black\",\n',
      '          box(width = 12,height = df_select$boxsize,title = paste(unique(df_select$slice)[{i}], \'group\'),\n',
      '            status = \"primary\",solidHeader = TRUE,\n',
      '            shiny::plotOutput(outputId = \'coexpfeature{i}\',height = df_select$image_size)),\n',
      '          fluidRow(column(6,br(),\n',
      '              downloadButton(\"coexpfeature{i}.pdf\", \"Download PDF\"),\n',
      '              downloadButton(\"coexpfeature{i}.png\", \"Download PNG\")),\n',
      '            column(6,div(style = \"display:inline-block\",\n',
      '                numericInput(\"coexpfeature{i}.h\",\"PDF / PNG height:\",width = \"138px\",\n',
      '                  min = 4,max = 20,value = 6,step = 0.5)),\n',
      '              div(style = \"display:inline-block\",\n',
      '                numericInput(\"coexpfeature{i}.w\",\"PDF / PNG width:\",\n',
      '                  width = \"138px\",min = 4,max = 20,value = 12,step = 0.5)))))),\n\n',
    )
    all_cops <- paste(all_cops,cop)
  }

  ends_6 <- glue::glue(
    '      ),\n',
    '    column(4,box(width = 12,shiny::plotOutput(outputId = \'coexp_leg\')),br(),DTOutput(outputId = \'coexp_table\')))\n',
    '  ))\n',
    '))))\n\n'
  )

  return(glue::glue(page6_ui_head,all_cops,ends_6))
}
