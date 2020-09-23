# constants
stats_cols <- c("Annotation", "# gene sets", "# genes", "# matches", "P-value",
                "Adjusted P-value", "Odds Ratio", "Fold Enrichment",
                "Adjusted Fold Enrichment")
option_bar <- stats_cols[2:9]
option_heat <- c("Gene Value", stats_cols[2:9])
trans_color <- c("No adjustment" = "identity", "Exponential" = "exp",
                 "Logarithmic" = "log", "Reciprocal" = "reciprocal",
                 "Square root" = "sqrt")
trans_value <- trans_color[-2]

if (requireNamespace("shiny", quietly = T)) {
  pane_input <- htmltools::tagList(
    shiny::selectInput("input.source", "Source",
                       c("Input" = "gene", "Seurat" = "cell")),
    shiny::textAreaInput("input.area", NULL, height = "calc(100vh - 468.5px)",
                         placeholder = "GENE1 0.1\nGENE2 0.2\n...",
                         resize = "vertical"),
    shiny::verbatimTextOutput("input.count", placeholder = T),
    shiny::numericInput("input.universe", "Genes in universe", 0, min = 0)
  )
  
  pane_home <- shiny::tabPanel(
    NULL, htmltools::br(), icon = shiny::icon("home", lib = "glyphicon"),
    shiny::selectInput("anno.source", "Annotations", NULL),
    shiny::selectInput(
      "anno.types", NULL,
      c("Gene Sets" = "name", "Symbols" = "syms", "Descriptions" = "info",
        "Automatic" = "auto", "Manual" = "file"), selected = "file",
      multiple = T
    ),
    shiny::textInput("anno.regex", NULL,
                     placeholder = "Filter with regular expressions"),
    shiny::verbatimTextOutput("anno.count", placeholder = T),
  
    shiny::selectInput("data.source", "Database", NULL),
    shiny::selectInput("data.category", NULL, NULL, multiple = T),
    shiny::selectInput("data.organism", NULL, NULL, multiple = T),
    shiny::verbatimTextOutput("data.count", placeholder = T),
  
    shiny::selectInput("cell.source", "Seurat data", NULL),
    shiny::selectInput("cell.initial", NULL, NULL),
    shiny::selectInput("cell.compare", NULL, NULL),
    shiny::verbatimTextOutput("cell.count", placeholder = T)
  )
  
  pane_advanced <- shiny::tabPanel(
    NULL, htmltools::br(), icon = shiny::icon("cog", lib = "glyphicon"),
    shiny::radioButtons(
      "info.source", "Description source",
      c("Annotations" = "anno", "Database" = "sets"), inline = T
    ),
    shiny::checkboxInput("name.clean",
                         "Remove underscores from gene set names", value = T),
    htmltools::hr(),
    shiny::varSelectInput("info.cols", "Information", NULL, multiple = T),
    shiny::varSelectInput("stat.cols", "Statistics", NULL, NULL, multiple = T),
    htmltools::hr(),
    shiny::fileInput("upload", NULL, buttonLabel = "Upload..."),
    shiny::downloadButton("download", "Download", class = "wide_button"),
    htmltools::tags$style(".wide_button {width: 100%}")
  )
  
  pane_bars <- shiny::tabPanel(
    NULL, htmltools::br(), icon = shiny::icon("stats", lib = "glyphicon"),
    shiny::selectInput("bars.values", "Value", option_bar,
                       selected = "Fold Enrichment"),
    shiny::selectInput("bars.values.trans", NULL, trans_value),
    shiny::selectInput("bars.colors", "Colour", option_bar,
                       selected = "Adjusted P-value"),
    shiny::selectInput("bars.colors.trans", NULL, trans_color),
    shiny::selectInput("bars.anno.order", "Sort", stats_cols),
    shiny::uiOutput("bars.anno.select"),
    shiny::checkboxInput("bars.anno.sort", "Sort plot by value")
  )
  
  pane_overlap <- shiny::tabPanel(
    NULL, htmltools::br(), icon = shiny::icon("equalizer", lib = "glyphicon"),
    shiny::selectInput("over.colors", "Colour", option_heat),
    shiny::selectInput("over.colors.trans", NULL, trans_color),
    shiny::helpText("The Odds Ratio, Fold Enrichment and Adjusted Fold \
                    Enrichment options may hide matches"),
    shiny::selectInput("over.anno.order", "Annotation order", stats_cols),
    shiny::uiOutput("over.anno.select"),
    shiny::selectInput(
      "over.gene.order", "Gene order",
      c("As entered" = "input", "Gene name" = "names", "Gene value" = "value")
    ),
    shiny::uiOutput("over.gene.select")
  )
  
  pane_cell <- shiny::tabPanel(
    NULL, htmltools::br(), icon = shiny::icon("screenshot", lib = "glyphicon"),
    shiny::selectInput("cell.plot", "Plot type", c(
      "Dot Plot" = "dot",
      "Feature Plot" = "feat",
      "Heatmap" = "heat",
      "Ridge Plot" = "ridge",
      "Violin Plot" = "violin"
    )),
    shiny::selectInput("cell.anno", "Annotation", NULL),
    shiny::selectInput("cell.genes", "Genes", NULL, multiple = T),
    shiny::checkboxInput("cell.plot.cluster", "Cluster"),
    shiny::numericInput("cell.downsample", "Downsample", 500, 1, 1000),
    htmltools::hr(),
    shiny::selectInput("cell.cluster", "Reduction", NULL)
  )
}
