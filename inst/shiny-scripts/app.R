library(shiny)
library(shinyBS) # for interactive elements
library(dplyr)
library(qPCRAnalyzer)

# Constants to define
DT_OPTIONS <- list(
  dom = "Brtip",
  scrollX = TRUE,
  scrollCollapse = TRUE,
  pageLength = 12
)

# Define the UI for the application
ui <- fluidPage(
  # Since fileInput labels are detailed instructions, custom CSS
  # is used to make the labels not bold
  tags$style(HTML(".control-label { font-weight: normal; }")),
  tags$style(HTML("h5 { font-weight: bold;
                        font-style: italic;
                        font-size: 1.1em; }")),

  # App title
  titlePanel(tags$h1(tags$b("qPCR Analyzer"),
                     "Streamlined qPCR data analysis")),

  # === Part 1: Data Input ===

  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      h3("Data upload"),

      # File input
      fileInput(
        inputId = "file", 
        label = "Upload qCPR data file. Must be a tabular file with columns for 
                 Well, Target, Group, Sample, and Cq.",
        accept = c(".xls", ".xlsx", ".csv", ".txt", ".tsv")
      ),

      uiOutput("sidebar_analysis_header"),
      uiOutput("sidebar_analysis_ref_gene"),
      uiOutput("sidebar_analysis_run_button"),

      # Download button
      uiOutput("download_ui"),
    ),


    # Main panel for displaying outputs
    mainPanel(
      tabsetPanel(
        tabPanel(tags$b("Data Summary and Validation"), 
          # Output: Data summary
          h3("Data summary"),
          helpText("The following is a summary of the uploaded data once provided."),
          br(),
          tableOutput("data_summary"),
          

          br(),
          br(),

          # === Part 2: Well Validation and Exclusion ===
          h3("Well Validation and Exclusion"),
          helpText("- Wells that do not have valid values and those that have
                    outlier Cq values are highlighted in the table below."),
          helpText("- NA wells will be automatically excluded from the analysis."),
          helpText("- Users should manually examine Cq values for outlier wells
                    and determine if they should be excluded from the analysis."),
          helpText("- Sample-Target pairs with outlier wells are displayed at
                    the front of the table."),
          br(),
          uiOutput("na_wells_header"),
          uiOutput("na_wells"),
          uiOutput("na_wells_summary"),

          br(),
          uiOutput("outlier_wells_header"),
          uiOutput("outlier_help_text"),
          
          DT::dataTableOutput("data_table"),
          uiOutput("wells_to_exclude"),
          uiOutput("wells_to_exclude_count"),

          br(),
          helpText(tags$b("Confirm the wells to exclude before running the analysis.")),        
        ),

        # === Part 3: Data Analysis ===
        tabPanel(tags$b("Analysis Results"),
          h3("Analysis Results"),
          helpText("Select the reference gene to view the relative expression analysis results."),
          uiOutput("select_ref_gene"),
          br(),
          plotOutput("expr_plots"),
        )
      ),
    ),
  ),
)


server <- function(input, output, session) {

  # Disable main output tabs until the analysis is performed
  hideTab(
    inputId = "main_panel",
    target = "Analysis Results"
  )

  data_reactive <- reactiveVal(NULL)

  # === Part 1: Data Input (Server) ===
  # Reactive expression for the data file
  observeEvent(input$file, {
    req(input$file)
    df <- import_and_validate(input$file$datapath) %>%
      highlight_outliers()
    data_reactive(df)
  })

  # Compute a data summary table
  output$data_summary <- renderTable({
    data <- data_reactive()
      if (!is.null(data)) {
        # Create a summary table
        summary_table <- data.frame(
          Category = c("Targets", "Samples", "Groups"),
          Count = c(
            length(unique(data$Target)),  # Number of unique targets
            length(unique(data$Sample)),  # Number of unique samples
            length(unique(data$Group))    # Number of unique groups
          ),
          UniqueValues = c(
            paste(unique(data$Target), collapse = ", "),  # List of unique targets
            paste(unique(data$Sample), collapse = ", "),  # List of unique samples
            paste(unique(data$Group), collapse = ", ")    # List of unique groups
          ),
        stringsAsFactors = FALSE
        )
      
        # Return the summary table
        return(summary_table)
      }
    }, rownames = FALSE)

  # === Part 2: Well Validation and Exclusion ===

  # NA wells header
  output$na_wells_header <- renderUI({
    if (!is.null(data_reactive())) {
      tags$h5("Wells with NA Cq values:")
    }
  })

  # Display wells with NA Cq values
  output$na_wells <- renderTable({
    data <- data_reactive()
    if (!is.null(data)) {
      na_wells <- data %>%
        get_invalid_wells()
      na_wells <- na_wells[["na_wells"]]
      return(na_wells)
    }
  })
  output$na_wells_summary <- renderText({
    data <- data_reactive()
    if (!is.null(data)) {
      n_na_wells <- sum(is.na(data$Cq))
      if (n_na_wells > 0) {
        return(paste(n_na_wells, "wells with NA Cq values found."))
      } else {
        return("No wells with NA Cq values found.")
      }
    }
  })

  # Outlier well header
  output$outlier_wells_header <- renderUI({
    if (!is.null(data_reactive())) {
      tags$h5("Wells with potential outlier values:")
    }
  })

  output$outlier_help_text <- renderUI({
    if (!is.null(data_reactive())) {
      tags$p("Select wells to exclude from the analysis by clicking on the rows in the table below.")
    }
  })

  # Interactively display a data table with highlighted outlier wells
  # Organize the table to display outlier wells first
  # Keep track of selected wells for exclusion
  selected_wells <- reactiveVal(character(0))
  output$data_table <- DT::renderDataTable({
    data <- data_reactive()
    
    if (!is.null(data)) {
      # Keep the automatically determined outlier wells at the front
      data <- rbind(data[data$Outlier_Sample_no_na, ],
                    data[!data$Outlier_Sample_no_na, ]) %>%
        select(Well, Target, Group, Sample, Cq, Outlier)

      # Preselect outlier wells
      preselected_rows <- which(data$Outlier == TRUE)

      # Generate the table with interactive row selection
      DT::datatable(data,
                    selection = list(target = "row", selected = preselected_rows),  # Preselect outliers
                    options = DT_OPTIONS,
                    filter = "top", 
                    rownames = FALSE) %>%
        DT::formatStyle("Outlier",            # Highlight only the "Outlier" column
                        backgroundColor = DT::styleEqual(c(TRUE, FALSE), c("darkred", "transparent")))  # Red for TRUE, transparent for FALSE
    }
  })

  # Track selected wells for exclusion (update when the selection changes)
  observe({
    selected <- input$data_table_rows_selected  # Get the indices of selected rows

    if (is.null(selected) | length(selected) == 0) {
      selected_wells(character(0))  # Reset the selected wells
    } else {
      # Update the reactive value with the selected wells
      data <- data_reactive()
      data <- rbind(data[data$Outlier_Sample_no_na, ],
                    data[!data$Outlier_Sample_no_na, ]) %>%
        pull(Well)

      # Get the corresponding well names from selected rows
      wells_to_exclude <- data[selected]
      selected_wells(wells_to_exclude)  # Store the selected wells for exclusion
    }
  })

  # Display wells selected to be excluded
  output$wells_to_exclude <- renderUI({
    if (length(selected_wells()) > 0) {
      tagList(
        tags$b("Selected wells to exclude:"),
        tags$span(paste(selected_wells(), collapse = ", "))
      )
    } else if (!is.null(data_reactive())) {
      tags$b("No wells selected for exclusion.")
    }
  })


  # Display number of wells selected for exclusion
  output$wells_to_exclude_count <- renderText({
    if (length(selected_wells()) > 0) {
      return(paste(tags$b("Number of wells selected for exclusion:"), length(selected_wells())))
    }
  })


  # === Part 3: Data Analysis ===

  # Sidebar layouts
  output$sidebar_analysis_header <- renderUI({
    if (!is.null(data_reactive())) {
      tagList(
        br(),
        h3("Data analysis"),
        helpText("Once data is uploaded, check the summary and validate the wells before running the analysis.")
      )
    }
  })

  output$sidebar_analysis_ref_gene <- renderUI({
    if (!is.null(data_reactive())) {
      tagList(
        br(),
        helpText("Select the reference gene(s) for relative expression analysis:"),
        selectInput("ref_gene", 
                    label = "Reference Gene(s):", 
                    choices = rev(unique(data_reactive()$Target)), 
                    multiple = TRUE)
      )
    }
  })

  output$sidebar_analysis_run_button <- renderUI({
    if (!is.null(data_reactive())) {

      tagList(
        br(),
        bsButton("run_analysis", 
                 label = "Run Analysis", 
                 block = TRUE, 
                 style = "primary"),
        tags$b(tags$p("Confirm the wells to exclude before running the analysis."))
      )
    }
  })

  # Run the analysis

  result_lst <- reactiveVal(NULL)
  
  observeEvent(input$run_analysis, {
    req(data_reactive())

    # Do not proceed if no reference gene is selected
    # Provide a warning message if no reference gene is selected
    if (length(input$ref_gene) == 0) {
      showModal(modalDialog(
        title = "Warning",
        "Please select at least one reference gene for analysis.",
        easyClose = TRUE
      ))
      return(NULL)
    }

    # Perform the analysis
    analysis_result <- data_reactive() %>%
      exclude_invalid_wells(selected_wells()) %>%
      calculate_relative_expression(ref = input$ref_gene)

    # Store the analysis results
    result_lst(analysis_result)
  })

  # Conditionally display the download button only when result_lst() is available
  output$download_ui <- renderUI({
    req(result_lst())
    tagList(
      br(),
      h3("Download Analysis Results"),
      downloadButton("download_analysis", "Download qPCR Analysis Results")
    )
  })

  # Enable download button
  output$download_analysis <- downloadHandler(
    filename = function() {
      paste0(gsub("-", "", Sys.Date()), "_qPCR.xlsx")
    },
    content = function(file) {
      openxlsx::write.xlsx(result_lst(), file)
    }
  )

  # Once the analysis is performed, display the tab to view the results
  observeEvent(result_lst(), {
    showTab(
      inputId = "main_panel",
      target = "Analysis Results"
    )
  })

  # Display the select input for reference gene
  output$select_ref_gene <- renderUI({
    req(result_lst())

    tagList(
      selectInput("ref_gene_select", 
                  label = "Select reference gene to view the results:", 
                  choices = names(result_lst())[-1],
                  selected = names(result_lst())[2])
    )
  })

  # According to the selected reference gene, display the plots
  output$expr_plots <- renderPlot({
    req(result_lst())

    # Get the selected reference gene
    ref_gene <- input$ref_gene_select

    # Get the results for the selected reference gene
    results <- result_lst()[ref_gene]

    # Create a list of plots
    plots_list <- plot_expr(result_lst())[[ref_gene]]

    if (is.null(plots_list) || length(plots_list) == 0) {
        return(NULL)  # Return nothing if there are no plots
    }

    # Determine the number of columns for the plots
    n_plots <- length(plots_list)
    group_by <- ifelse(length(unique(data_reactive()$Group)) > 1, "Group", "Sample")
    if (n_plots == 1) {
      n_cols <- 1
    } else if (n_plots == 2 | group_by == "Sample") {
      n_cols <- 2
    } else {
      n_cols <- 3
    }

    # Make a combined plot
    combined_plot <- cowplot::plot_grid(plotlist = plots_list, ncol = n_cols)

    return(combined_plot)
  })








}




      





























shinyApp(ui, server)