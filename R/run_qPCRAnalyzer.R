# Run the Shiny App
run_qPCRAnalyzer <- function() {
  shiny::runApp(system.file("shiny-scripts", package = "qPCRAnalyzer"),
    display.mode = "normal"
  )
}