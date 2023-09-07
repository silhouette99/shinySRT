#' perpare the code for shinyspatial
#'
#'
#' #' Generate code files required for shiny app containing only one dataset. In
#' particular, two R scripts will be generated, namely \code{server.R} and
#' \code{ui.R}. If users want to include multiple dataset in one shiny app.
#'  Note that both \code{preparedata_shinyspatial} and \code{prepare_code}
#'  functions are ran when
#' running the wrapper function \code{makespashiny}.
#'
#'
#' @param shiny.dir save the shiny code and data
#'
#' @return server.R and ui.R required for shiny app
#'
#' @import data.table readr glue
#'
#'
#' @examples
#' prepare_code(shiny.dir = 'shinyspatial_app')
#'
#' @export
prepare_code <- function(shiny.dir = 'shinyspatial_app'){

  filename = paste0(shiny.dir, "/server.R")
  df_select <- readRDS(paste0(shiny.dir,"/df_select.Rds"))
  # library
  readr::write_file(lib_server(),file = filename)
  # color
  readr::write_file(col_server(),file = filename,append = T)
  # load data
  readr::write_file(data_server(),file = filename,append = T)
  # function
  readr::write_file(fun_server(),file = filename,append = T)
  # page
  readr::write_file(server_heads(df_select = df_select),file = filename,append = T)
  readr::write_file(page1_server(df_select = df_select),file = filename,append = T)
  readr::write_file(page2_server(df_select = df_select),file = filename,append = T)
  readr::write_file(page3_server(df_select = df_select),file = filename,append = T)
  readr::write_file(page4_server(df_select = df_select),file = filename,append = T)
  readr::write_file(page5_server(df_select = df_select),file = filename,append = T)
  readr::write_file(page6_server(df_select = df_select),file = filename,append = T)

  ##
  filename = paste0(shiny.dir, "/ui.R")
  readr::write_file(ui_load(), file = filename)
  readr::write_file(
    glue::glue(
      'container <- function(...) {{\n',
      '  shiny::fluidRow(shiny::column(...))\n',
      '}}\n\n\n'
    ),
    file = filename,
    append = T
  )

  readr::write_file(ui_head(), file = filename, append = T)

  ## page1
  readr::write_file(page1_ui_head(), file = filename, append = T)
  uip1 <- ui_p1_p(df_select = df_select)
  for (i in uip1) {
    readr::write_file(i, file = filename, append = T)}
  readr::write_file(p1_tab_ui(), file = filename, append = T)
  ## page2
  readr::write_file(page2_ui_head(), file = filename, append = T)
  uip2 <- ui_p2_p(df_select = df_select)
  for (i in uip2) {
    readr::write_file(i, file = filename, append = T)
  }
  ## page3
  readr::write_file(page3_ui_head(), file = filename,append = T)
  readr::write_file(p3_vg_plot(df_select = df_select), file = filename,append = T)
  ## page4
  readr::write_file(page4_ui(df_select = df_select), file = filename,append = T)
  ## page5
  readr::write_file(page5_ui(df_select = df_select), file = filename,append = T)
  ## page6
  readr::write_file(page6_ui(df_select = df_select), file = filename,append = T)
}
