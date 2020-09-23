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

view_plot <- . %>% shiny::plotOutput(height = "calc(100vh - 132.5px)")
well_text <- . %>%
  shiny::textOutput() %>%
  shiny::wellPanel(style = "height: calc(100vh - 170px); overflow-y: auto;")
view_text <- function(head, text)
  shiny::column(6, htmltools::h5(head), well_text(text))
view_comp <- function(lh, lc, rh, rc)
  shiny::fluidRow(view_text(lh, lc), view_text(rh, rc))

view_bars <- shiny::tabPanel(
  "Plot", icon = shiny::icon("stats", lib = "glyphicon"), view_plot("bars")
)
view_over <- shiny::tabPanel(
  "Overlap", icon = shiny::icon("equalizer", lib = "glyphicon"),
  view_plot("over")
)
view_cell <- shiny::tabPanel(
  "Cluster", icon = shiny::icon("screenshot", lib = "glyphicon"),
  view_plot("cell")
)
view_heat <- shiny::tabPanel(
  "Expression", icon = shiny::icon("fire", lib = "glyphicon"),
  view_plot("heat")
)
view_stat <- shiny::tabPanel(
  "Statistics", icon = shiny::icon("th", lib = "glyphicon"),
  shiny::dataTableOutput("stat")
)
view_gene <- shiny::tabPanel(
  "Gene List", icon = shiny::icon("ok-circle", lib = "glyphicon"),
  view_comp("Recognised", "gene_ok", "Not recognised", "gene_no")
)
view_sets <- shiny::tabPanel(
  "Gene Sets", icon = shiny::icon("list", lib = "glyphicon"),
  view_comp("Not in database", "uniq_anno", "Not in annotations", "uniq_sets")
)
view_vals <- shiny::tabPanel(
  "Gene Values", icon = shiny::icon("tasks", lib = "glyphicon"),
  view_comp("Assigned values", "vals_ok", "Not assigned values", "vals_no")
)
view_info <- shiny::tabPanel(
  "Information", icon = shiny::icon("info-sign", lib = "glyphicon"),
  shiny::dataTableOutput("info")
)

options <- shiny::tabsetPanel(type = "pills", pane_home, pane_advanced,
                              pane_bars, pane_overlap, pane_cell)
sidebar <- shiny::sidebarPanel(
  shiny::fluidRow(shiny::column(6, pane_input), shiny::column(6, options))
)
content <- shiny::mainPanel(
  shiny::tabsetPanel(view_bars, view_over, view_cell, view_heat, view_stat,
                     view_gene, view_sets, view_vals, view_info)
)

ui <- shiny::fluidPage(
  theme = shinythemes::shinytheme("yeti"),
  shinyjs::useShinyjs(),
  shiny::titlePanel("GLACIER: Gene List Annotation, Calculation and \
                    Illustration of Enrichment in R"),
  shiny::sidebarLayout(sidebar, content)
)
