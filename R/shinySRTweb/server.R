library(shiny)
library(shinyhelper)
library(shinydashboard)
library(data.table)
library(Matrix)
library(DT)
library(scales)
library(ggiraph)
library(Seurat)
require(shinySRT)

source('source/loadsrt.R')




server <- function(input, output, session) {
  # input file size
  options(shiny.maxRequestSize = 500 * 1024 ^ 3)
  dirsmake <- generate_random_string(10)
  # shinySRT
  output$page_web <- renderUI({
    main_page()
  })
  
  loaded_obj <- shiny::reactiveValues(obj = NULL)
  # tmp_dirs <- shiny::reactiveValues(dirs = NULL)
  file_stat <- shiny::reactiveValues(obj = T)
  
  observeEvent(input$run_button_obj, {
    if (file_stat$obj) {
      if (!is.null(input$obj_input1) | !is.null(input$obj_input2)) {
        if (!is.null(input$obj_input1)) {
          obj_in1 <- readRDS(input$obj_input1$datapath)
        } else{
          obj_in1 <- NULL
        }
        if (!is.null(input$obj_input2)) {
          obj_in2 <- readRDS(input$obj_input2$datapath)
        } else{
          obj_in2 <- NULL
        }
        
        objs <-
          load_spatial_obj_web(
            obj1 = obj_in1,
            obj2 = obj_in2,
            species = input$sp,
            resolution = input$resv,
            npcs = input$npcs
          )
        # loaded_obj$obj <- objs
        # if (class(objs) == 'Seurat') {
        CreateshinySRT(
          dat = objs,
          maxlevel = input$level_line,
          gene.mapping = input$genemap,
          sp_normalize = F,
          shiny.dir = dirsmake,
          web = T
        )
        
        file_stat$obj <- F
        # }
        
      }
    }
  })
  
  
  observeEvent(input$run_button_10x, {
    if (file_stat$obj) {
      sample1 <-
        !is.null(input$mtx_10x_input1) &
        !is.null(input$image_10x_input1) &
        !is.null(input$json_10x_input1) & !is.null(input$position_10x_input1)
      sample2 <-
        !is.null(input$mtx_10x_input2) &
        !is.null(input$image_10x_input2) &
        !is.null(input$json_10x_input2) & !is.null(input$position_10x_input2)
      if (sample1 | sample2) {
        if (sample1) {
          matx1 <- input$mtx_10x_input1$datapath
          imgs1 <- input$image_10x_input1$datapath
          sclfct1 <- input$json_10x_input1$datapath
          pos1 <- input$position_10x_input1$datapath
          if (length(grep(pattern = '.mtx', matx1)) > 0) {
            if (length(input$gene_10x_input1) == 0 &
                length(input$barcode_10x_input1) == 0) {
              validate('Lack genes file and cell names!')
            } else{
              gene1 <- input$gene_10x_input1$datapath
              cell1 <- input$barcode_10x_input1$datapath
            }
          } else{
            gene1 <- NULL
            cell1 <- NULL
          }
        } else{
          matx1 <- NULL
          imgs1 <- NULL
          sclfct1 <- NULL
          pos1 <- NULL
          gene1 <- NULL
          cell1 <- NULL
        }
        
        
        if (sample2) {
          matx2 <- input$mtx_10x_input2$datapath
          imgs2 <- input$image_10x_input2$datapath
          sclfct2 <- input$json_10x_input2$datapath
          pos2 <- input$position_10x_input2$datapath
          if (length(grep(pattern = '.mtx', matx2)) > 0) {
            if (length(input$gene_10x_input2) == 0 &
                length(input$barcode_10x_input2) == 0) {
              validate('Lack genes file and cell names!')
            } else{
              gene2 <- input$gene_10x_input2$datapath
              cell2 <- input$barcode_10x_input2$datapath
            }
          } else{
            gene2 <- NULL
            cell2 <- NULL
          }
        } else{
          matx2 <- NULL
          imgs2 <- NULL
          sclfct2 <- NULL
          pos2 <- NULL
          gene2 <- NULL
          cell2 <- NULL
        }
        # loaded_obj$obj <- c(matx1,sclfct1,imgs1,pos1,
        #                     gene1,
        #                     cell1)
        objs <-
          load_spatial_10x_web(
            mat1 = matx1,
            mat_r1 = gene1 ,
            mat_c1 = cell1 ,
            image1 = imgs1 ,
            json1 = sclfct1 ,
            coordinates1 = pos1,
            mat2 = matx2,
            mat_r2 = gene2 ,
            mat_c2 = cell2 ,
            image2 = imgs2 ,
            json2 = sclfct2 ,
            coordinates2 = pos2,
            species = input$sp,
            npcs = input$npcs,
            resolution = input$resv
          )
        
        # loaded_obj$obj <- objs
        CreateshinySRT(
          dat = objs,
          maxlevel = input$level_line,
          gene.mapping = input$genemap,
          sp_normalize = F,
          shiny.dir = dirsmake,
          web = T
        )
        
        file_stat$obj <- F
      }
      
      
    }
  })
  
  
  observeEvent(input$run_button_mtx, {
    if (file_stat$obj) {
      sample1 <-
        !is.null(input$mtx_input1) &
         !is.null(input$position_input1)
      sample2 <-
        !is.null(input$mtx_input2) &
        !is.null(input$position_input2)
      if (sample1 | sample2) {
        if (sample1) {
          matx1 <- input$mtx_input1$datapath
          imgs1 <- input$image_input1$datapath
          pos1 <- input$position_input1$datapath
          if (!is.null(input$meta_input1)) {
            metas1 <- input$meta_input1$datapath
          } else{
            metas1 <- NULL
          }
        } else{
          matx1 <- NULL
          imgs1 <- NULL
          pos1 <- NULL
          metas1 <- NULL
        }
        
        
        if (sample2) {
          matx2 <- input$mtx_input2$datapath
          imgs2 <- input$image_input2$datapath
          pos2 <- input$position_input2$datapath
          if (!is.null(input$meta_input2)) {
            metas2 <- input$meta_input2$datapath
          } else{
            metas2 <- NULL
          }
        } else{
          matx2 <- NULL
          imgs2 <- NULL
          pos2 <- NULL
          metas2 <- NULL
        }
        
        objs <-
          load_spatial_mat_web(
            matx1 = matx1,
            meta1 = metas1,
            coordi1 = pos1,
            image1 = imgs1,
            x_reverse1 = input$x_rev1,
            matx2 = matx2,
            meta2 = metas2,
            coordi2 = pos2,
            image2 = imgs2,
            x_reverse2 = input$x_rev2,
            species = input$sp,
            npcs = input$npcs,
            resolution = input$resv
          )
        
        # loaded_obj$obj <- objs
        # if (class(objs) == 'Seurat') {
        CreateshinySRT(
          dat = objs,
          maxlevel = input$level_line,
          gene.mapping = input$genemap,
          sp_normalize = F,
          shiny.dir = dirsmake,
          web = T
        )
        # }
        file_stat$obj <- F
      }
    }
  })
  
  
  output[['fullPage']] <- renderUI({
    if (length(intersect(dir(dirsmake), 'server.R')) > 0) {
      appUI <- parse(file = paste0(dirsmake, '/ui.R'))
      appServer <- eval(parse(file = paste0(dirsmake, '/server.R')))
      if (!is.null(input$go) && input$go > 0) {
        # If they pressed the button once,
        # run the appServer function and evaluate the parsed appUI code
        appServer(input, output, session)
        eval(appUI)
      } else {
        #
        NULL
      }
    }
  })
  
  observeEvent(input$go, {
    output[['fullPage']] <- renderUI({
      if (length(intersect(dir(dirsmake), 'server.R')) > 0) {
        appUI <- parse(file = paste0(dirsmake, '/ui.R'))
        appServer <-
          eval(parse(file = paste0(dirsmake, '/server.R')))
        if (!is.null(input$go) && input$go > 0) {
          # If they pressed the button once,
          # run the appServer function and evaluate the parsed appUI code
          appServer(input, output, session)
          eval(appUI)
        } else {
          #
          NULL
        }
      }
    })
  })
  
  observeEvent(input$res_web, {
    session$reload()
  })
  
  
  
  session$onSessionEnded(function() {
    cat("Session Ended\n")
    # unlink(paste0(dirsmake, '/image.Rds'))
    # unlink(paste0(dirsmake, '/data.h5'))
    # unlink(paste0(dirsmake, '/df_select.Rds'))
    # unlink(paste0(dirsmake, '/genesets.Rds'))
    # unlink(paste0(dirsmake, '/meta.Rds'))
    # unlink(paste0(dirsmake, '/meta_group.Rds'))
    # unlink(isolate(dirsmake), recursive = TRUE)
    unlink(dirsmake, recursive = TRUE)
    
  })
  
  
}
