library(shiny)
library(shinyjs)
library(shinythemes)

inputPane <- tagList(
  selectInput("input.source", "Source", c("Text Input" = "text", "Seurat" = "cell")),
  textAreaInput("input.text", NULL, resize = "vertical", height = "calc(100vh - 393.5px)", placeholder = "GENE1 0.05\nGENE2 0.1\n..."),
  verbatimTextOutput("input.count", T),
  numericInput("input.universe", "Genes in universe", 0, 0)
)

homePane <- tabPanel(
  NULL, br(), icon = icon("home", lib = "glyphicon"),
  selectInput("anno.source", "Annotations", NULL),
  selectInput("anno.types", NULL, c("Gene Sets" = "name", "Symbols" = "syms", "Descriptions" = "info", "Automatic" = "auto", "Manual" = "file"), selected = "file", multiple = T),
  textInput("anno.regex", NULL, placeholder = "Filter using regular expressions"),
  verbatimTextOutput("anno.count", T),
  
  selectInput("data.source", "Database", NULL),
  selectInput("data.categories", NULL, NULL, multiple = T),
  selectInput("data.organisms", NULL, NULL, multiple = T),
  verbatimTextOutput("data.count", T),
  
  selectInput("cell.source", "Seurat data", NULL),
  selectInput("cell.select", NULL, NULL),
  selectInput("cell.compare", NULL, NULL),
  verbatimTextOutput("cell.count", T)
)

etcPane <- tabPanel(
  NULL, br(), icon = icon("cog", lib = "glyphicon"),
  radioButtons("info.source", "Description source", c("Annotations" = "anno", "Database" = "data"), inline = T),
  checkboxInput("name.fix", "Remove underscores from gene set names", T), hr(),
  
  varSelectInput("info.columns", "Information", NULL, multiple = T),
  varSelectInput("stat.columns", "Statistics", NULL, multiple = T), hr(),
  
  fileInput("file.up", NULL, buttonLabel = "Upload..."),
  downloadButton("file.down", "Download", class = "wide"),
  tags$head(tags$style(".wide {width: 100%;}"))
)

barsPane <- tabPanel(
  NULL, br(), icon = icon("stats", lib = "glyphicon"),
  selectInput("bars.value", "Value", c("# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), selected = "Fold Enrichment"),
  selectInput("bars.value.trans", NULL, c("No adjustment" = "identity", "Logarithmic" = "log", "Reciprocal" = "reciprocal", "Square root" = "sqrt")),
  selectInput("bars.color", "Color", c("# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), selected = "Adjusted P-value"),
  selectInput("bars.color.trans", NULL, c("No adjustment" = "identity", "Exponential" = "exp", "Logarithmic" = "log", "Reciprocal" = "reciprocal", "Square root" = "sqrt")),
  selectInput("bars.anno.order", "Order", c("Annotation", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment")),
  uiOutput("bars.anno.select"),
  checkboxInput("bars.anno.sort", "Sort plot by value")
)

overPane <- tabPanel(
  NULL, br(), icon = icon("equalizer", lib = "glyphicon"),
  selectInput("over.color", "Color", c("Gene Value", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), "Adjusted P-value"),
  selectInput("over.color.trans", NULL, c("No adjustment" = "identity", "Exponential" = "exp", "Logarithmic" = "log", "Reciprocal" = "reciprocal", "Square root" = "sqrt")),
  helpText("Selecting Gene Value, Odds Ratio, Fold Enrichment or Adjusted Fold Enrichment options may hide matches"), br(),
  selectInput("over.anno.order", "Annotation order", c("Annotation", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment")),
  uiOutput("over.anno.select"),
  selectInput("over.gene.order", "Gene order", c("As entered" = "input", "Gene name" = "names", "Gene value" = "value")),
  uiOutput("over.gene.select")
)

cellPane <- tabPanel(
  NULL, br(), icon = icon("screenshot", lib = "glyphicon"),
  selectInput("cell.overview", "Overview plot", NULL),
  numericInput("cell.downsample", "Downsample", 50, 1), hr(),
  
  selectInput("cell.plot", "Expression plot", c("Dot Plot" = "dot", "Feature Plot" = "feat", "Heatmap" = "heat", "Ridge Plot" = "ridge", "Violin Plot" = "violin"), "feat"),
  selectInput("cell.anno", "Annotation", NULL),
  selectInput("cell.genes", "Genes", NULL, multiple = T),
  checkboxInput("cell.gene.match", "Restrict to genes in input", T),
  checkboxInput("cell.gene.cluster", "Cluster genes based on expression")
)

oneColumn <- function(head, body) column(6, h5(head), wellPanel(textOutput(body), class = "text-compare"))
twoColumn <- function(l_head, l_body, r_head, r_body) fluidRow(oneColumn(l_head, l_body), oneColumn(r_head, r_body))
barsView <- tabPanel("Plot", icon = icon("stats", lib = "glyphicon"), plotOutput("bars", height = "calc(100vh - 132.5px)"))
overView <- tabPanel("Overlap", icon = icon("equalizer", lib = "glyphicon"), plotOutput("over", height = "calc(100vh - 132.5px)"))
cellView <- tabPanel("Cluster", icon = icon("screenshot", lib = "glyphicon"), plotOutput("cell", height = "calc(100vh - 132.5px)"))
heatView <- tabPanel("Expression", icon = icon("fire", lib = "glyphicon"), plotOutput("heat", height = "calc(100vh - 132.5px)"))
statView <- tabPanel("Statistics", icon = icon("th", lib = "glyphicon"), dataTableOutput("stat"))
geneView <- tabPanel("Gene List", icon = icon("ok-circle", lib = "glyphicon"), twoColumn("Recognised", "gene.ok", "Not recognised", "gene.no"))
gsetView <- tabPanel("Gene Sets", icon = icon("list", lib = "glyphicon"), twoColumn("Gene sets not in annotations", "data.only", "Gene sets not in annotations", "anno.only"))
valsView <- tabPanel("Gene Values", icon = icon("tasks", lib = "glyphicon"), twoColumn("Genes with values", "vals.ok", "Genes without values", "vals.no"))
infoView <- tabPanel("Information", icon = icon("info-sign", lib = "glyphicon"), dataTableOutput("info"))

ui <- fluidPage(
  theme = shinytheme("yeti"), useShinyjs(),
  titlePanel("GLACIER: Gene List Annotation, Calculation and Illustration of Enrichment in R"),
  sidebarLayout(
    sidebarPanel(fluidRow(column(6, inputPane), column(6, tabsetPanel(type = "pills", homePane, etcPane, barsPane, overPane, cellPane)))),
    mainPanel(tabsetPanel(barsView, overView, cellView, heatView, statView, geneView, gsetView, valsView, infoView))
  ),
  tags$head(tags$style(".text-compare {height: calc(100vh - 170px); overflow-y: auto;}"))
)
