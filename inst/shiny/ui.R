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
  selectInput("anno.source", "Annotations", c("No annotations" = "")),
  selectInput("anno.types", NULL, c("Gene Sets" = "name", "Symbols" = "syms", "Descriptions" = "info", "Automatic" = "auto", "Manual" = "file"), selected = "file", multiple = T),
  textInput("anno.regex", NULL, placeholder = "Filter using regular expressions"),
  verbatimTextOutput("anno.count", T),
  
  selectInput("data.source", "Database", c("No databases" = "")),
  selectInput("data.categories", NULL, c("Filter categories" = ""), multiple = T),
  selectInput("data.organisms", NULL, c("Filter organisms" = ""), multiple = T),
  verbatimTextOutput("data.count", T),
  
  selectInput("cell.source", "Seurat data", c("No Seurat data" = "")),
  selectInput("cell.cluster", NULL, c("Inter-cluster" = "_inter_")),
  selectInput("cell.select", NULL, c("No clusters" = "")),
  selectInput("cell.compare", NULL, c("No clusters" = "")),
  verbatimTextOutput("cell.count", T)
)

etcPane <- tabPanel(
  NULL, br(), icon = icon("cog", lib = "glyphicon"),
  radioButtons("info.source", "Description source", c("Annotations" = "anno", "Database" = "data"), inline = T),
  checkboxInput("stat.sigonly", "Significant results only", T),
  checkboxInput("name.fix", "Remove underscores from gene set names", T), hr(),
  
  varSelectizeInput("info.columns", "Information", NULL, multiple = T, options = list(placeholder = "Select columns to display")),
  varSelectizeInput("stat.columns", "Statistics", NULL, multiple = T, options = list(placeholder = "Select columns to display")), hr(),
  
  p("MSigDB XML files available ", a("here", href = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/")),
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
  selectInput("cell.overview", "Overview plot", c("No Seurat data selected" = "")),
  numericInput("cell.downsample", "Downsample", 500, 3), hr(),
  
  selectInput("cell.plot", "Expression plot", c("Dot Plot" = "dot", "Feature Plot" = "feat", "Heatmap" = "heat", "Ridge Plot" = "ridge", "Violin Plot" = "violin"), "feat"),
  selectInput("cell.anno", "Annotation (adjusted P-value ≤ 0.05 only)", c("No annotations loaded" = "")),
  selectInput("cell.genes", "Genes (click to select more)", c("Select genes to display" = ""), multiple = T),
  checkboxInput("cell.gene.match", "Restrict to genes in input", T),
  checkboxInput("cell.gene.cluster", "Cluster genes based on expression")
)

twoTable <- function(left, right) fluidRow(column(6, dataTableOutput(left)), column(6, dataTableOutput(right)))
barsView <- tabPanel("Plot", icon = icon("stats", lib = "glyphicon"), plotOutput("bars", height = "calc(100vh - 132.5px)"))
overView <- tabPanel("Overlap", icon = icon("equalizer", lib = "glyphicon"), plotOutput("over", height = "calc(100vh - 132.5px)"))
cellView <- tabPanel("Cluster", icon = icon("screenshot", lib = "glyphicon"), plotOutput("cell", height = "calc(100vh - 132.5px)"))
heatView <- tabPanel("Expression", icon = icon("fire", lib = "glyphicon"), plotOutput("heat", height = "calc(100vh - 132.5px)"))
statView <- tabPanel("Statistics", icon = icon("th", lib = "glyphicon"), dataTableOutput("stat"))
contView <- tabPanel("Contents", icon = icon("book", lib = "glyphicon"), twoTable("data.gene", "cell.gene"))
geneView <- tabPanel("Recognised", icon = icon("ok", lib = "glyphicon"), twoTable("gene.ok", "gene.no"))
gsetView <- tabPanel("Mismatched", icon = icon("remove", lib = "glyphicon"), twoTable("data.only", "anno.only"))
valsView <- tabPanel("Gene Values", icon = icon("tasks", lib = "glyphicon"), twoTable("vals.ok", "vals.no"))
infoView <- tabPanel("Information", icon = icon("info-sign", lib = "glyphicon"), dataTableOutput("info"))

ui <- fluidPage(
  theme = shinytheme("yeti"), useShinyjs(),
  titlePanel("glacier: Gene List Annotation, Calculation and Illustration of Enrichment in R"),
  sidebarLayout(
    sidebarPanel(fluidRow(column(6, inputPane), column(6, tabsetPanel(type = "pills", homePane, etcPane, barsPane, overPane, cellPane)))),
    mainPanel(tabsetPanel(barsView, overView, cellView, heatView, statView, contView, geneView, gsetView, valsView, infoView))
  ),
  tags$head(tags$style(".text-compare {height: calc(100vh - 170px); overflow-y: auto;}"))
)
