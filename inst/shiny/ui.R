library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyWidgets)

library(stringr)

# constants
stats_cols <- c("Annotation", "# gene sets", "# genes", "# matches", "P-value",
                "Adjusted P-value", "Odds Ratio", "Fold Enrichment",
                "Adjusted Fold Enrichment")
option_bars <- stats_cols[2:9]
option_over <- c("Gene Value", stats_cols[2:9])
trans_color <- c("No adjustment" = "identity", "Exponential" = "exp",
                 "Logarithmic" = "log", "Reciprocal" = "reciprocal",
                 "Square root" = "sqrt")
trans_value <- trans_color[-2]

pane_input <- tagList(
  selectInput("input.source", "Source", c("Input" = "gene", "Seurat" = "cell")),
  textAreaInput("input.area", NULL, height = "calc(100vh - 468.5px)",
                placeholder = "GENE1 0.1\nGENE2 0.2\n...", resize = "vertical"),
  verbatimTextOutput("input.count", placeholder = T),
  numericInput("input.universe", "Genes in universe", 0, min = 0)
)

pane_home <- tabPanel(
  NULL, br(), icon = icon("home", lib = "glyphicon"),
  selectInput("anno.source", "Annotations", NULL),
  selectInput(
    "anno.types", NULL,
    c("Gene Sets" = "name", "Symbols" = "syms", "Descriptions" = "info",
      "Automatic" = "auto", "Manual" = "file"), selected = "file",
    multiple = T
  ),
  textInput("anno.regex", NULL, placeholder = "Filter with regex"),
  verbatimTextOutput("anno.count", placeholder = T),

  selectInput("data.source", "Database", NULL),
  selectInput("data.category", NULL, NULL, multiple = T),
  selectInput("data.organism", NULL, NULL, multiple = T),
  verbatimTextOutput("data.count", placeholder = T),

  selectInput("cell.source", "Seurat data", NULL),
  selectInput("cell.initial", NULL, NULL),
  selectInput("cell.compare", NULL, NULL),
  verbatimTextOutput("cell.count", placeholder = T)
)

pane_advanced <- tabPanel(
  NULL, br(), icon = icon("cog", lib = "glyphicon"),
  radioButtons("info.source", "Description source",
               c("Annotations" = "anno", "Database" = "sets"), inline = T),
  checkboxInput("name.clean", "Remove underscores from gene set names",
                value = T),
  hr(),
  varSelectInput("info.cols", "Information", NULL, multiple = T),
  varSelectInput("stat.cols", "Statistics", NULL, NULL, multiple = T),
  hr(),
  fileInput("upload", NULL, buttonLabel = "Upload..."),
  downloadButton("download", "Download", class = "wide_button"),
  tags$style(".wide_button {width: 100%}")
)

pane_bars <- tabPanel(
  NULL, br(), icon = icon("stats", lib = "glyphicon"),
  selectInput("bars.values", "Value", option_bars,
              selected = "Fold Enrichment"),
  selectInput("bars.values.trans", NULL, trans_value),
  selectInput("bars.colors", "Colour", option_bars,
              selected = "Adjusted P-value"),
  selectInput("bars.colors.trans", NULL, trans_color),
  selectInput("bars.anno.order", "Sort", stats_cols),
  uiOutput("bars.anno.select"),
  checkboxInput("bars.anno.sort", "Sort plot by value")
)

pane_overlap <- tabPanel(
  NULL, br(), icon = icon("equalizer", lib = "glyphicon"),
  selectInput("over.colors", "Colour", option_over),
  selectInput("over.colors.trans", NULL, trans_color),
  helpText("The Gene Value, Odds Ratio, Fold Enrichment and Adjusted Fold \
           Enrichment options may hide matches"),
  selectInput("over.anno.order", "Annotation order", stats_cols),
  uiOutput("over.anno.select"),
  selectInput(
    "over.gene.order", "Gene order",
    c("As entered" = "input", "Gene name" = "names", "Gene value" = "value")
  ),
  uiOutput("over.gene.select")
)

pane_cell <- tabPanel(
  NULL, br(), icon = icon("screenshot", lib = "glyphicon"),
  selectInput("cell.plot", "Plot type", c(
    "Dot Plot" = "dot",
    "Feature Plot" = "feat",
    "Heatmap" = "heat",
    "Ridge Plot" = "ridge",
    "Violin Plot" = "violin"
  )),
  selectInput("cell.anno", "Annotation", NULL),
  selectInput("cell.genes", "Genes", NULL, multiple = T),
  checkboxInput("cell.plot.cluster", "Cluster"),
  numericInput("cell.downsample", "Downsample", 500, 1, 1000),
  hr(),
  selectInput("cell.cluster", "Reduction", NULL)
)

view_plot <- function(name) plotOutput(name, height = "calc(100vh - 132.5px)")
well_text <- function(name) wellPanel(
  textOutput(name), style = "height: calc(100vh - 170px); overflow-y: auto;"
)

view_text <- function(head, text) column(6, h5(head), well_text(text))
view_comp <- function(lh, lc, rh, rc)
  fluidRow(view_text(lh, lc), view_text(rh, rc))
view_table <- function(head, table)
  column(6, h5(head), tableOutput(table), tags$style(
    str_c("#", table, " {height: calc(100vh - 170px); overflow-y: auto}")))
view_tcomp <- function(lh, lc, rh, rc)
  fluidRow(view_table(lh, lc), view_table(rh, rc))

view_bars <-
  tabPanel("Plot", icon = icon("stats", lib = "glyphicon"), view_plot("bars"))
view_over <- tabPanel(
  "Overlap", icon = icon("equalizer", lib = "glyphicon"), view_plot("over")
)
view_cell <- tabPanel(
  "Cluster", icon = icon("screenshot", lib = "glyphicon"), view_plot("cell")
)
view_heat <- tabPanel(
  "Expression", icon = icon("fire", lib = "glyphicon"), view_plot("heat")
)
view_stat <- tabPanel(
  "Statistics", icon = icon("th", lib = "glyphicon"), dataTableOutput("stat")
)
view_gene <- tabPanel(
  "Gene List", icon = icon("ok-circle", lib = "glyphicon"),
  view_comp("Recognised", "gene_ok", "Not recognised", "gene_no")
)
view_sets <- tabPanel(
  "Gene Sets", icon = icon("list", lib = "glyphicon"),
  view_comp("Not in database", "uniq_anno", "Not in annotations", "uniq_sets")
)
view_vals <- tabPanel(
  "Gene Values", icon = icon("tasks", lib = "glyphicon"),
  view_tcomp("Assigned values", "vals_ok", "Not assigned values", "vals_no")
)
view_info <- tabPanel(
  "Information", icon = icon("info-sign", lib = "glyphicon"),
  dataTableOutput("info")
)

options <- tabsetPanel(type = "pills", pane_home, pane_advanced, pane_bars,
                       pane_overlap, pane_cell)
sidebar <-
  sidebarPanel(fluidRow(column(6, pane_input), column(6, options)))
content <-
  mainPanel(tabsetPanel(view_bars, view_over, view_cell, view_heat, view_stat,
                        view_gene, view_sets, view_vals, view_info))

ui <- fluidPage(
  theme = shinytheme("yeti"),
  useShinyjs(),
  titlePanel("GLACIER: Gene List Annotation, Calculation and \
                    Illustration of Enrichment in R"),
  sidebarLayout(sidebar, content)
)
