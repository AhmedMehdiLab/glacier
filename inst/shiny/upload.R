RAND_SEED <- 444

uploadUI <- function(id) {ns <- NS(id)
  modalDialog(
    # configuration
    title = "Upload file", size = "l", easyClose = T, useShinyjs(),
    footer = tagList(modalButton("Cancel"),
                     actionButton("confirm", "Confirm")),
    
    # styling
    tags$style(str_glue("#{ns('anno')} {{height: 500px; overflow: auto;}}")),
    tags$style(str_glue("#{ns('cell')} {{height: 500px; overflow: auto;}}")),
    tags$style(str_glue("#{ns('data')} {{height: 500px; overflow: auto;}}")),
    tags$style(str_glue("#{ns('gene')} {{height: 500px; overflow: auto;}}")),
    
    # interface
    fluidRow(
      column(3, h5("Options"),
        selectInput(ns("type"), NULL, c("File type" = "", "Annotations" = "anno", "Database" = "data", "MSigDB XML" = "msig", "Seurat RDS" = "cell")),
        textInput(ns("name"), NULL, placeholder = "File name"),
        selectInput(ns("delim"), NULL, c("Commas" = ",", "Tabs" = "\t")),
        checkboxInput(ns("header"), "File has header"),
        sliderInput(ns("content"), "Content columns", 1, 10, c(1, 10)),
        numericInput(ns("info"), "Description column", 0, 10, 1),
        numericInput(ns("down"), "Downsample", 50, 1),
        helpText("If no descriptions are available, set 'Description column' to 0"), br(),
        helpText("Please wait until preview has loaded before clicking 'Confirm'")),
      column(9, h5("Preview"),
        uiOutput(ns("preview")))
    )
  )
}

uploadServer <- function(id, file, values) moduleServer(id,
  function(input, output, session) {
    ns <- session$ns
    
    observeEvent(input$type, {
      path <- file()$datapath
      name <- file_path_sans_ext(file()$name)
      output$preview <- NULL
      
      # toggle controls
      is_data <- input$type %in% c("data", "msig")
      is_delim <- input$type %in% c("anno", "data")
      toggleState("delim", is_delim)
      toggleState("header", is_delim)
      toggleState("content", is_delim)
      toggleState("info", is_delim)
      toggleState("ok", input$type != "")
      toggleState("down", input$type == "cell")
      updateTextInput(session, "name", value = name)
      
      # parse data and update controls
      if (input$type == "cell") proc <- reactive(subset(readRDS(path), downsample = input$down, seed = RAND_SEED))
      if (input$type == "msig") proc <- reactive(import_msigdb(path))
      if (is_delim) {
        temp <- reactive(glacier:::import_delim_path(path, input$delim, input$header))
        cols <- ncol(temp())
        
        updateSliderInput(session, "content", NULL, c(1, cols), 1, cols)
        updateNumericInput(session, "info", NULL, 0, 0, cols)
      }
      
      # prepare preview
      if (input$type == "anno") {
        proc <- reactive(glacier:::import_annotations_file(temp(), input$content, input$info))
        output$preview <- renderUI(tableOutput(ns("anno")))
      }
      if (input$type == "cell") output$preview <- renderUI(verbatimTextOutput(ns("cell")))
      if (input$type == "data") proc <- reactive(glacier:::import_database_file(temp(), input$content, input$info))
      if (is_data) output$preview <- renderUI(tagList(column(6, tableOutput(ns("data"))), column(6, verbatimTextOutput(ns("gene")))))
      
      # update preview
      values$name <- input$name
      values$type <- if (input$type == "msig") "data" else input$type
      values$proc <- reactive(proc())
      output$anno <- renderTable(proc() %>% head(16) %>% select(all_of(1:2), last_col(1:0)))
      output$cell <- renderText(str_c(capture.output(proc()[[]]), collapse = "\n"))
      output$data <- renderTable(proc()$gs_info %>% head(16))
      output$gene <- renderText(str_c(capture.output(str(proc()$gs_genes, comp.str = NULL, give.attr = F, give.head = T, vec.len = 1, list.len = 22)), collapse = "\n"))
    })
  }
)
