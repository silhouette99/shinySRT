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
    'library(ggiraph)\n',
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
    '      menuItem(\"Heat / Dot plot\", tabName = \"heatdot\", icon = icon(\"gear\"))),\n',
    '    ### import new spot information\n',
    '    shiny::selectizeInput(inputId = \"new_title\",\n',
    '      label = \'group name\',choices = NULL,\n',
    '      multiple = FALSE,options = list(create = TRUE)),\n',
    '    textAreaInput(\"impspifo\",HTML(\'Import new spots information\'),\n',
    '      height = \"200px\",value = NULL),\n',
    '    actionButton(inputId = \'swi_spinfo\',\n',
    '      label = \'Import\',icon = icon(\'house\')),\n',
    '    sliderInput(inputId = \'point_size\',\n',
    '      label = \'spots size\',min = 0,\n',
    '      ticks = FALSE,max = 1,value = 1)),\n\n'
  )
}
## page1 head
ui_p1 <- function(df_select){
  ui_p1_head <- c('  dashboardBody(tabItems(tabItem(tabName = \"dashboard\",\n',
                  '     h2(strong(HTML(\"SpotInfo vs GeneExpr\"))),\n',
                  '     h4(\"Spot information vs gene expression on reduced dimensions\"),\n',
                  '     \"In this tab, users can visualise both Spot information and gene \",\n',
                  '     \"expression side-by-side on slice spatial representions.\",br(),br(),\n',
                  '     fluidRow(column(12,style = \"border-bottom: 2px solid black\",br(),br(),\n',
                  '      ### first row:the control panel\n',
                  '       fluidRow(column(width = 2,br(),br(),\n',
                  '         shiny::selectizeInput(inputId = \"meta1_select\",\n',
                  '           label = \'select group\',choices = meta_group$group,\n',
                  '           multiple = FALSE,selected = df_select$meta1,\n',
                  '           options = list(create = TRUE),width = \"50%\"),\n',
                  '           shiny::selectInput(inputId = \"gene\",\n',
                  '             label = \"gene names:\",\n',
                  '             choices = NULL,width = \"50%\"),\n',
                  '           shiny::actionButton(inputId = \"reset1\",\n',
                  '             icon = icon(\"refresh\"),label = \"reset all\")),\n',
                  '           ### spot information selected (1th row 1th column)\n',
                  '           column(3,br(),br(),h4(\"spot numbers / statistics\"),\n',
                  '             radioButtons(\"splt\",\"Split continuous spot info into:\",\n',
                  '               choices = c(\"Quartile\", \"Decile\"),selected = \"Decile\",inline = TRUE),\n',
                  '             shiny::sliderInput(inputId = \'stat_thrd\',label = \"Gene expression threshold\",\n',
                  '               min = 0,max = 1,value = 0)),\n',
                  '           column(7,box(title = \'Statistical table\',collapsed = F,width = 12,\n',
                  '               align = \"center\",status = \"primary\",\n',
                  '               solidHeader = TRUE,collapsible = TRUE,\n',
                  '               DTOutput(outputId = \'statis_gene\')))),br(),\n')
  
  ui_p1_head <- paste(ui_p1_head,collapse = '')
  slices <- length(unique(df_select$slice))
  if(slices <= 4){
    nr = 1
    if(slices <= 2){
      nc = 6
    }else{
      nc = 12/slices
    }
    
  }else{
    nr = ceiling(slices/4)
    nc = 3
  }
  
  ## bgplot
  bg_for <- lapply(1:nr, function(l){
    eachr <- 1:4+4*(l-1)
    eachr <- eachr[eachr <= slices]
    bg_for <- lapply(eachr, function(eachc){
      if(eachc == max(eachr)){
        m =')' 
      }else{
        m = ','
      }
      
      glue::glue(
        '           column({nc},style = \"border-right: 2px solid black\",\n',
        '            box(title = paste(unique(df_select$slice)[{eachc}], \'group\'),collapsed = F,\n',
        '              width = 12,align = \"center\",status = \"primary\",\n',
        '              solidHeader = TRUE,collapsible = TRUE,\n',
        '              ggiraph::girafeOutput(outputId = \'bg_plot{eachc}\',\n',
        '               height = df_select$image_size)),\n',
        '               fluidRow(column(6,br(),\n',
        '                 downloadButton(\"bg_spot{eachc}.pdf\", \"Download PDF\"),\n',
        '                 downloadButton(\"bg_spot{eachc}.png\", \"Download PNG\")),\n',
        '               column(6,\n',
        '                 div(style = \"display:inline-block\",\n',
        '                  numericInput(\"bg_spot{eachc}.h\", \"PDF / PNG height:\",\n',
        '                   width = \"138px\",min = 4,max = 20, value = 6,step = 0.5)),\n',
        '                 div(style = \"display:inline-block\",\n',
        '                  numericInput(\"bg_spot{eachc}.w\",\"PDF / PNG width:\",\n',
        '                   width = \"138px\", min = 4,max = 20,value = 12,step = 0.5))))){m}\n'
      )
      
    }) %>% unlist()%>%paste(collapse = '')
    # if(l == max(nr)){
    #   bg_for <- paste('fluidRow(\n',bg_for,'),',sep = '')
    # }else{
    #   bg_for <- paste('fluidRow(\n',bg_for,',',sep = '')
    # }
    bg_for <- paste('fluidRow(\n',bg_for,',',sep = '')
    
  }) %>% unlist()%>%paste(collapse = '')
  bg_for <- paste('fluidRow(',bg_for,sep = '')
  
  ## feature
  ft_for <- lapply(1:nr, function(l){
    eachr <- 1:4+4*(l-1)
    eachr <- eachr[eachr <= slices]
    ft_for <- lapply(eachr, function(eachc){
      if(eachc == max(eachr)){
        m ='))),' 
      }else{
        m = '),'
      }
      
      glue::glue(
        '             column({nc},style = \"border-right: 2px solid black\",\n',
        '              box(width = 12,collapsible = T,collapsed = F,\n',
        '               title = paste(unique(df_select$slice)[{eachc}], \'feature\'),\n',
        '               status = \"primary\",solidHeader = TRUE,\n',
        '               ggiraph::girafeOutput(outputId = \'feature{eachc}\',height = df_select$image_size)),\n',
        '               fluidRow(column(6,br(),\n',
        '                downloadButton(\"ft_spot{eachc}.pdf\", \"Download PDF\"),\n',
        '                downloadButton(\"ft_spot{eachc}.png\", \"Download PNG\")),\n',
        '               column(6,\n',
        '                  div(style = \"display:inline-block\",\n',
        '                    numericInput(\"ft_spot{eachc}.h\",\"PDF / PNG height:\",\n',
        '                     width = \"138px\", min = 4,max = 20,value = 6,step = 0.5)),\n',
        '                  div(style = \"display:inline-block\",\n',
        '                    numericInput(\"ft_spot{eachc}.w\",\"PDF / PNG width:\",\n',
        '                     width = \"138px\",min = 4,max = 20,value = 12,step = 0.5)))){m}\n\n',
      )
      
    }) %>% unlist()%>%paste(collapse = '')
    ft_for <- paste('fluidRow(\n',ft_for)
  }) %>% unlist()%>%paste(collapse = '')
  ft_for <- paste(ft_for,'         br(),br()))),\n',sep = '')
  ui_p1 <- paste(ui_p1_head,bg_for,ft_for)%>%glue::glue()
  return(ui_p1)
}


## page2 head
ui_p2 <- function(df_select) {
  ui_p2_head <-
    c(
      '      tabItem(tabName = \"d_gene_ex\",h2(strong(HTML(\"GeneExpr vs GeneExpr\"))),\n',
      '       h4(\'Gene expression vs gene expression on spatial represention\'),br(),br(),\n',
      '       fluidRow(column(12, style = \"border-bottom: 2px solid black\",br(),br(),\n',
      '        fluidRow(column(6,fluidRow(column(width = 6,br(),br(),\n',
      '           shiny::selectizeInput(inputId = \"meta2_select\",label = \'select group\',choices = meta_group$group,\n',
      '            multiple = FALSE,selected = df_select$meta1,options = list(create = TRUE),width = \"50%\"),\n',
      '           shiny::actionButton(inputId = \"reset2\",icon = icon(\"refresh\"),label = \"reset all\")),\n',
      '           column(6,br(),br(),shiny::selectizeInput(inputId = \"group_tab\",\n',
      '             label = \'Add group\',choices = c(\'No\', meta_group[info == T, group]),\n',
      '             multiple = FALSE,selected = \'No\',\n',
      '             options = list(create = TRUE)))))),br(),\n\n'
    )
  
  ui_p2_head <- paste(ui_p2_head,collapse = '')
  
  slices <- length(unique(df_select$slice))
  if (slices <= 4) {
    nr = 1
    if (slices <= 2) {
      nc = 6
    } else{
      nc = 12 / slices
    }
  } else{
    nr = ceiling(slices / 4)
    nc = 3
  }
  slices <- length(unique(df_select$slice))
  if (slices <= 4) {
    nr = 1
    if (slices <= 2) {
      nc = 6
    } else{
      nc = 12 / slices
    }
    
  } else{
    nr = ceiling(slices / 4)
    nc = 3
  }
  ## bgplot
  bg_for <- lapply(1:nr, function(l) {
    eachr <- 1:4+4*(l-1)
    eachr <- eachr[eachr <= slices]
    bg_for <- lapply(eachr, function(eachc) {
      if (eachc == max(eachr)) {
        m = '),'
      } else{
        m = ','
      }
      
      glue::glue(
        '           column({nc},style = \"border-right: 2px solid black\",\n',
        '           box(width = 12,collapsible = T,\n',
        '            title = paste(unique(df_select$slice)[{eachc}], \'spot selected\'),\n',
        '            status = \"primary\",solidHeader = TRUE,\n',
        '             ggiraph::girafeOutput(outputId = \'bg_plot_f{eachc}\',height = df_select$image_size))){m}\n\n',
      )
      
    }) %>% unlist() %>% paste(collapse = '')
    bg_for <- paste('fluidRow(\n', bg_for, sep = '')
  }) %>% unlist() %>% paste(collapse = '')
  bg_for <- paste('fluidRow(', bg_for, sep = '')
  
  coex_leg <- c(
    '          fluidRow(column(12,column(4,column(6,fluidRow(shiny::radioButtons(\"col_sel\",\"color select:\",\n',
    '              choices = c(\"Red (Gene1); Blue (Gene2)\",\"Orange (Gene1); Blue (Gene2)\",\n',
    '                          \"Red (Gene1); Green (Gene2)\",\"Green (Gene1); Blue (Gene2)\"),\n',
    '              selected = \"Red (Gene1); Blue (Gene2)\",inline = TRUE)),fluidRow(\n',
    '                shiny::selectInput(inputId = \"gene_f1\",label = \"gene1 names:\",choices = NULL,width = \"50%\"),\n',
    '                shiny::selectInput(inputId = \"gene_f2\",label = \"gene2 names:\",choices = NULL,width = \"50%\")\n',
    '              )),column(6,\n',
    '             shiny::sliderInput(inputId = \'thrd1\',label = \"Gene1 expression threshold\",min = 0, max = 1,value = 0),\n',
    '            shiny::sliderInput(inputId = \'thrd2\',label = \"Gene2 expression threshold\",min = 0,max = 1,value = 0))),\n',
    '            column(3,box(width = 10,height = 300,shiny::plotOutput(outputId = \'coexp_leg\'))),\n',
    '            column(5,box(collapsed = F,width = 12, align = \"center\",status = \"primary\",solidHeader = TRUE,collapsible = TRUE,DTOutput(outputId = \'coexp_table\'))))),\n\n'
  )
  coex_leg <- paste(coex_leg,collapse = '')
  
  co_for <- lapply(1:nr, function(l) {
    eachr <- 1:4+4*(l-1)
    eachr <- eachr[eachr <= slices]
    co_for <- lapply(eachr, function(eachc) {
      if (eachc == max(eachr)) {
        m = ')'
      } else{
        m = ','
      }
      glue::glue(
        '           column({nc},style = \"border-right: 2px solid black\",\n',
        '           box(width = 12,collapsible = T,\n',
        '             title = unique(df_select$slice)[{eachc}],\n',
        '             status = \"primary\",solidHeader = TRUE,\n',
        '             ggiraph::girafeOutput(outputId = \'coexpfeature{eachc}\',height = df_select$image_size)),\n',
        '           fluidRow(\n',
        '             column(6,br(),\n',
        '               downloadButton(\"coexpfeature{eachc}.pdf\", \"Download PDF\"),\n',
        '               downloadButton(\"coexpfeature{eachc}.png\", \"Download PNG\")),\n',
        '             column(6,\n',
        '               div(style = \"display:inline-block\",\n',
        '                numericInput(\"coexpfeature{eachc}.h\", \"PDF / PNG height:\",\n',
        '                 width = \"138px\",min = 4,\n',
        '                 max = 20,value = 6,step = 0.5)),\n',
        '               div(style = \"display:inline-block\",\n',
        '                numericInput(\"coexpfeature{eachc}.w\", \"PDF / PNG width:\",\n',
        '                 width = \"138px\",min = 4,\n',
        '                 max = 20,value = 12,step = 0.5))))){m}\n\n',
      )
      
    }) %>% unlist() %>% paste(collapse = '')
    co_for <- paste('fluidRow(\n', co_for, sep = '')
  }) %>% unlist() %>% paste(collapse = '')
  co_for <- paste(co_for, ')))),\n\n', sep = '')
  
  ui_p2 <- paste(ui_p2_head,bg_for,coex_leg,co_for,sep = '')%>%glue::glue()
  return(ui_p2)
}

## page3
ui_p3 <- function(df_select){
  ui_p3_head <- c(
    '      tabItem(tabName = \"vio\",h2(strong(HTML(\"Viol / Box data chart\"))),\n',
    '        h4(\'Spot information / gene expression violin plot / box plot\'),\n',
    '        \'In this tab, users can visualise the gene expression or continuous spot information (e.g. Number of UMIs / module score) across groups of spots (e.g. libary / clusters).\',\n',
    '        br(),br(),\n',
    '        fluidRow(column(12, \n',
    '        fluidRow(\n',
    '         column(3,style = \"border-right: 2px solid black\",\n',
    '          box(width = 12,br(),br(),\n',
    '           fluidRow(\n',
    '            shiny::selectizeInput(\n',
    '              inputId = \"vio_x\",label = \'Spot information - x\',choices = meta_group[info == T, group],\n',
    '              selected = df_select$meta1,multiple = FALSE,options = list(create = TRUE))),br(),\n',
    '           fluidRow(\n',
    '            shiny::selectInput(choices = NULL,\n',
    '                multiple = FALSE,inputId = \'vio_y\',\n',
    '                label = \'statistical info - y\')),br(),\n',
    '           fluidRow(\n',
    '            shiny::selectizeInput(choices = c(\'No\', meta_group[info == T, group]),\n',
    '              selected = \'No\',inputId = \'vio_gb\',label = \'group by\',options = list(create = TRUE))),\n',
    '            shiny::radioButtons(\"plot_type\",\"Plot type:\",\n',
    '              choices = c(\"violin\", \"boxplot\"),selected = \"violin\",\n',
    '               inline = TRUE),\n',
    '            checkboxInput(\"point_ex\", \"Show data points\", value = FALSE),br(),\n',
    '             fluidRow(\n',
    '              shiny::selectizeInput(\n',
    '                inputId = \"subspot\",\n',
    '                label = \'subset spot\',\n',
    '                choices = meta_group$group,\n',
    '                selected = \'seurat_clusters\',\n',
    '                multiple = F,options = list(create = TRUE))),br(),br(),\n',
    '            fluidRow(column(12,shiny::actionButton(inputId = \"reset3\",\n',
    '                icon = icon(\"refresh\"),label = \"reset all\"))),\n',
    '             fluidRow(h4(strong(\'Vlnplot/Boxplot\')),\n',
    '                column(12, fluidRow(\n',
    '                   column(6,br(),downloadButton(\"vb_spot.pdf\", \"Download PDF\"),\n',
    '                     downloadButton(\"vb_spot.png\", \"Download PNG\")),\n',
    '                   column(6,div(style = \"display:inline-block\",\n',
    '                     numericInput(\"vb_spot.h\",\"PDF / PNG height:\",\n',
    '                       width = \"138px\",min = 4,\n',
    '                       max = 20,value = 6,step = 0.5)),\n',
    '                   div(style = \"display:inline-block\",\n',
    '                     numericInput(\"vb_spot.w\",\"PDF / PNG width:\",\n',
    '                       width = \"138px\",min = 4,\n',
    '                       max = 20,value = 12,step = 0.5)))))))),\n',
    '         column(9,column(12,\n',
    '           box(width = 12,height = 780,status = \"primary\",\n',
    '             solidHeader = TRUE,ggiraph::girafeOutput(outputId = \'vb_plot\', height = 700))))),\n'
  )
  ui_p3_head <- paste(ui_p3_head,collapse = '')
  
  slices <- length(unique(df_select$slice))
  if(slices <= 4){
    nr = 1
    if(slices <= 2){
      nc = 6
    }else{
      nc = 12/slices
    }
    
  }else{
    nr = ceiling(slices/4)
    nc = 3
  }
  
  ## bgplot
  bg_for <- lapply(1:nr, function(l){
    eachr <- 1:4+4*(l-1)
    eachr <- eachr[eachr <= slices]
    bg_for <- lapply(eachr, function(eachc){
      if(eachc == max(eachr)){
        m =')' 
      }else{
        m = ','
      }
      
      glue::glue(
        '         column({nc},\n',
        '           box(width = 12,collapsible = T,\n',
        '            title = paste(unique(df_select$slice)[{eachc}], \'spot selected\'),\n',
        '            status = \"primary\",\n',
        '            solidHeader = TRUE,\n',
        '            ggiraph::girafeOutput(\n',
        '              outputId = \'vb_bg_plot{eachc}\',\n',
        '              height = df_select$image_size))){m}\n',
      )
      
    }) %>% unlist()%>%paste(collapse = '')
    if(l == max(nr)){
      bg_for <- paste('fluidRow(\n',bg_for,'',sep = '')
    }else{
      bg_for <- paste('fluidRow(\n',bg_for,',\n',sep = '')
    }
    
  }) %>% unlist()%>%paste(collapse = '')
  ui_p3 <- paste(ui_p3_head,'fluidRow(column(12,\n',bg_for,'))\n','))),',sep = '')
  return(ui_p3)
}

## page4
ui_p4 <- function(df_select){
  ui_p4_head <- c(
    '      tabItem(tabName = \"portion\",h2(strong(HTML(\"Portion data chart\"))),\n',
    '       h4(\'Proportion / spot numbers across different spot information\'),\n',
    '       \'In this tab, users can visualise the composition of spots based on one discrete spots information across another discrete spots information.\',br(),br(),\n',
    '       fluidRow(column(12, fluidRow(\n',
    '        column(3,style = \"border-right: 2px solid black\",\n',
    '         box(width = 12,br(),br(),\n',
    '           fluidRow(shiny::selectizeInput(\n',
    '            inputId = \"por_x\",label = \'Spot information - x\',\n',
    '            choices = meta_group[info == T, group],\n',
    '            selected = df_select$meta1,multiple = FALSE,\n',
    '            options = list(create = TRUE))),br(),\n',
    '           fluidRow(shiny::selectizeInput(\n',
    '            inputId = \"por_group\",label = \'color by information - group\',\n',
    '            choices = meta_group[info == T, group],\n',
    '            selected = df_select$meta2,multiple = FALSE,\n',
    '            options = list(create = TRUE))),br(),\n',
    '           shiny::radioButtons(\n',
    '            \"tab_type\",\"Plot type:\",\n',
    '            choices = c(\"percentage\", \"number\"),\n',
    '            selected = \"percentage\",inline = TRUE),\n',
    '           checkboxInput(\"flipxy\", \"Flip x y\", value = FALSE),br(),br(),\n',
    '           fluidRow(shiny::selectizeInput(\n',
    '            inputId = \"subspot2\",\n',
    '            label = \'subset spot\',\n',
    '            choices = meta_group[info == T, group],\n',
    '            selected = df_select$meta1,\n',
    '            multiple = F,options = list(create = TRUE))),br(),\n',
    '           fluidRow(column(12,shiny::actionButton(inputId = \"reset4\",\n',
    '            icon = icon(\"refresh\"),label = \"reset all\"))),br(),\n',
    '           fluidRow(h4(strong(\'Barplot\')), \n',
    '            column(12,fluidRow(column(6,br(),\n',
    '             downloadButton(\"bar_spot.pdf\", \"Download PDF\"),\n',
    '             downloadButton(\"bar_spot.png\", \"Download PNG\")),\n',
    '              column(6,\n',
    '                div(style = \"display:inline-block\",\n',
    '                 numericInput(\"bar_spot.h\",\n',
    '                  \"PDF / PNG height:\",\n',
    '                  width = \"138px\",\n',
    '                  min = 4,max = 20,\n',
    '                  value = 6,step = 0.5)),\n',
    '                div(style = \"display:inline-block\",\n',
    '                 numericInput(\"bar_spot.w\",\n',
    '                  \"PDF / PNG width:\",\n',
    '                  width = \"138px\",\n',
    '                  min = 4,max = 20,\n',
    '                  value = 12,step = 0.5)))))))),\n',
    '       column(9,fluidRow(column(12,box(width = 12,height = 780,\n',
    '         status = \"primary\",solidHeader = TRUE,\n',
    '         ggiraph::girafeOutput(outputId = \'por_plot\', height = 700)))))),\n'
  )
  
  ui_p4_head <- paste(ui_p4_head,collapse = '')
  slices <- length(unique(df_select$slice))
  if (slices <= 4) {
    nr = 1
    if (slices <= 2) {
      nc = 6
    } else{
      nc = 12 / slices
    }
    
  } else{
    nr = ceiling(slices / 4)
    nc = 3
  }
  ## bgplot
  bg_for <- lapply(1:nr, function(l) {
    eachr <- l * 1:4
    eachr <- eachr[eachr <= slices]
    bg_for <- lapply(eachr, function(eachc) {
      if (eachc == max(eachr)) {
        m = ')'
      } else{
        m = ','
      }
      
      glue::glue(
        '         column({nc},box(width = 12,collapsible = T,\n',
        '          title = paste(unique(df_select$slice)[{eachc}], \'spot selected\'),\n',
        '          status = \"primary\",solidHeader = TRUE,\n',
        '          ggiraph::girafeOutput(outputId = \'por_bg_plot{eachc}\',\n',
        '           height = df_select$image_size))){m}\n\n',
      )
      
    }) %>% unlist() %>% paste(collapse = '')
    if(l == max(nr)){
      bg_for <- paste('fluidRow(\n',bg_for,'\n',sep = '')
    }else{
      bg_for <- paste('fluidRow(\n',bg_for,',\n',sep = '')
    }
  }) %>% unlist() %>% paste(collapse = '')
  bg_for <- paste('fluidRow(column(12,\n', bg_for, '))))),',sep = '')       
  ui_p4 <- paste(ui_p4_head,bg_for,sep = '')%>%glue::glue()
  return(ui_p4)
}  


## page5
ui_p5 <- function(df_select){
  ui_p5_head <- c(
    '      tabItem(tabName = \"heatdot\",h2(strong(HTML(\"Gene expression heatmap/dotplot\"))),\n',
    '       h4(\'In this tab, users can visualise the gene expression patterns of multiple genes grouped by categorical spot information (e.g. library / cluster).\'),\n',
    '       \'The normalised expression are averaged, log-transformed and then plotted.\',br(),br(),\n',
    '       fluidRow(column(12, fluidRow(\n',
    '       column(3,style = \"border-right: 2px solid black\",\n',
    '       box(width = 12,br(),br(),\n',
    '        fluidRow(\n',
    '         textAreaInput(\"inpgenelist\", HTML(\"List of gene names <br />(Max 50 genes, separated <br />by , or ; or newline):\"),\n',
    '          height = \"200px\",value = paste0(df_select$genes, collapse = \", \")) %>%\n',
    '          helper(type = \"inline\",size = \"m\",fade = TRUE,\n',
    '           title = \"List of genes to plot on bubbleplot / heatmap\",\n',
    '           content = c(\"Input genes to plot\",\n',
    '            \"- Maximum 50 genes (due to ploting space limitations)\",\n',
    '            \"- Genes should be separated by comma, semicolon or newline\"))),br(),\n',
    '        fluidRow(shiny::selectizeInput(inputId = \"expgroup\",\n',
    '          label = \'Group by:\',choices = meta_group[info == T, group],\n',
    '          selected = df_select$meta1,multiple = FALSE,\n',
    '          options = list(create = TRUE))),br(),\n',
    '        shiny::radioButtons(\"exp_plot_type\",\n',
    '          \"Plot type:\",choices = c(\"dotplot\", \"heatmap\"),\n',
    '          selected = \"dotplot\",inline = TRUE),\n',
    '        checkboxInput(\"scl_exp\", \"Scale gene expression\", value = TRUE),\n',
    '        checkboxInput(\"rclust\", \"Cluster rows (gene)\", value = TRUE),\n',
    '        checkboxInput(\"cclust\", \"Cluster columns (sample)\", value = FALSE),br(),\n',
    '        fluidRow(shiny::selectizeInput(inputId = \"scgroup\",\n',
    '          label = \'Second Group by:\',choices = meta_group[info == T, group],\n',
    '          selected = NULL,multiple = FALSE,\n',
    '          options = list(create = TRUE))),br(),\n',
    '        fluidRow(shiny::selectizeInput(inputId = \"subspot3\",\n',
    '          label = \'subset spot\',choices = meta_group[info == T, group],\n',
    '          selected = df_select$meta1,\n',
    '          multiple = F,options = list(create = TRUE))),br(),\n',
    '        fluidRow(column(12,shiny::actionButton(inputId = \"reset5\",\n',
    '                 icon = icon(\"refresh\"),label = \"reset all\"))),\n',
    '        fluidRow(h4(strong(\'Dotplot / Heatmap\')), column(12,\n',
    '         fluidRow(column(6,br(),\n',
    '           downloadButton(\"heat_spot.pdf\", \"Download PDF\"),\n',
    '           downloadButton(\"heat_spot.png\", \"Download PNG\")),\n',
    '           column(6,\n',
    '            div(style = \"display:inline-block\",\n',
    '             numericInput(\"heat_spot.h\",\"PDF / PNG height:\",\n',
    '             width = \"138px\",min = 4,max = 20,value = 6,step = 0.5)),\n',
    '            div(style = \"display:inline-block\",\n',
    '             numericInput(\"heat_spot.w\",\"PDF / PNG width:\",\n',
    '             width = \"138px\", min = 4,max = 20,value = 12,step = 0.5)))))))),\n',
    '       column(9,fluidRow(column(12,\n',
    '           box(width = 12,height = 880,status = \"primary\",\n',
    '             solidHeader = TRUE,plotOutput(outputId = \'exp_plot\', height = 800)))))),\n'
  )
  
  ui_p5_head <- paste(ui_p5_head,collapse = '')
  
  slices <- length(unique(df_select$slice))
  if (slices <= 4) {
    nr = 1
    if (slices <= 2) {
      nc = 6
    } else{
      nc = 12 / slices
    }
    
  } else{
    nr = ceiling(slices / 4)
    nc = 3
  }
  ## bgplot
  bg_for <- lapply(1:nr, function(l) {
    eachr <- l * 1:4
    eachr <- eachr[eachr <= slices]
    bg_for <- lapply(eachr, function(eachc) {
      if (eachc == max(eachr)) {
        m = ')'
      } else{
        m = ','
      }
      
      glue::glue(
        '       column(6,box(width = 12,height = df_select$boxsize,collapsible = T,\n',
        '           title = paste(unique(df_select$slice)[{eachc}], \'spot selected\'),\n',
        '           status = \"primary\",solidHeader = TRUE,\n',
        '           ggiraph::girafeOutput(outputId = \'exp_bg_plot{eachc}\'))){m}\n'
      )
      
    }) %>% unlist() %>% paste(collapse = '\n')
    if(l == max(nr)){
      bg_for <- paste('fluidRow(\n',bg_for,'\n',sep = '')
    }else{
      bg_for <- paste('fluidRow(\n',bg_for,',\n',sep = '')
    }
  }) %>% unlist() %>% paste(collapse = '')
  
  bg_for <- paste(bg_for, '))))))',sep = '')       
  ui_p5 <- paste(ui_p5_head,bg_for,sep = '')%>%glue::glue()
  
  return(ui_p5)
}


