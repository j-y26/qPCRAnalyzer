#' @title Run the Shiny App
#' 
#' @description This function runs the qPCRAnalyzer Shiny app.
#' 
#' @export
#' 
#' @import shiny shinyBS openxlsx dplyr
#' 
run_qPCRAnalyzer <- function() {
  shiny::runApp(system.file("shiny-scripts", package = "qPCRAnalyzer"),
    display.mode = "normal"
  )
}