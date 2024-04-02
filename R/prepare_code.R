#' Prepare the code for shinySRT
#' 
#' 
#' Prepare the code of shinyApp for shinySRT, including ui and server. 
#' 
#'
#'
#' Generate code files required for shiny app containing only one dataset. In
#' particular, two R scripts will be generated, namely \code{server.R} and
#' \code{ui.R}. If users want to include multiple dataset in one shiny app.
#'  Note that both \code{preparedata_shinyspatial} and \code{prepare_code}
#'  functions are ran when
#' running the wrapper function \code{makespashiny}.
#'
#'
#' @param shiny.dir directory saving the shiny code and data
#' @param title the title of the given shinySRT
#'
#' @return
#'
#' @import data.table readr glue
#'
#'
#' @examples
#' prepare_code(shiny.dir = 'shinyspatial_app')
#'
#' @export
prepare_code <- function(shiny.dir = 'shinyspatial_app',title = 'spatial_example',web = F,tmpdir = '/srv/shiny-server/temp/'){
  if(web){
    filename = paste0(tmpdir,shiny.dir, "/server.R")
    df_select <- readRDS(paste0(tmpdir,shiny.dir,"/df_select.Rds"))
  }else{
    filename = paste0(shiny.dir, "/server.R")
    df_select <- readRDS(paste0(shiny.dir,"/df_select.Rds"))
  }
  
  # library
  readr::write_file(lib_server(),file = filename)
  # color
  readr::write_file(col_server(),file = filename,append = T)
  if(web){
    # load data
    readr::write_file(data_server_web(df_select = df_select,dir = shiny.dir,tmpdir = tmpdir),file = filename,append = T)
    # function
    # readr::write_file(fun_server_web(),file = filename,append = T)
    # page
  }else{
    # load data
    readr::write_file(data_server(df_select = df_select),file = filename,append = T)
    # function
    # readr::write_file(fun_server(),file = filename,append = T)
    # page
  }
  readr::write_file(fun_server(),file = filename,append = T)
  
  readr::write_file(server_heads(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p1(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p2(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p3(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p4(df_select = df_select),file = filename,append = T)
  readr::write_file(sp_server_p5(df_select = df_select),file = filename,append = T)
  
  if(length(grep(pattern = 'deconvolusion',dir(shiny.dir),value = T)) > 0){
    readr::write_file(sp_server_p6(df_select = df_select),file = filename,append = T)
  }else{
    readr::write_file('}',file = filename,append = T)
  }
  
  ##
  if(web){
    filename = paste0(tmpdir,shiny.dir, "/ui.R")
  }else{
    filename = paste0(shiny.dir, "/ui.R")
  }
  # filename = paste0(shiny.dir, "/ui.R")
  if(web){
    readr::write_file(ui_load_web(dir = shiny.dir,tmpdir = tmpdir), file = filename)
  }else{
    readr::write_file(ui_load(), file = filename)
  }
  
  readr::write_file(
    glue::glue(
      'container <- function(...) {{\n',
      '  shiny::fluidRow(shiny::column(...))\n',
      '}}\n\n\n'
    ),
    file = filename,
    append = T
  )
  
  readr::write_file(ui_head(title,df_select = df_select), file = filename, append = T)
  ## page1
  readr::write_file(ui_p1(df_select), file = filename, append = T)
  ## page2
  readr::write_file(ui_p2(df_select), file = filename, append = T)
  ## page3
  readr::write_file(ui_p3(df_select), file = filename,append = T)
  ## page4
  readr::write_file(ui_p4(df_select = df_select), file = filename,append = T)
  ## page5
  readr::write_file(ui_p5(df_select = df_select), file = filename,append = T)
  ## page6
  if(!is.null(df_select[['cell']])){
    readr::write_file(ui_p6(df_select = df_select), file = filename,append = T)
  }
  
  readr::write_file(glue::glue(')))'), file = filename,append = T)
}




