library(shiny)
library(shinyjs)
library(shinythemes)

inputPane <- tagList(
  selectInput("input.source", "Source", c("Text Input" = "text", "Seurat" = "cell")),
  textAreaInput("input.text", NULL, resize = "vertical", height = "calc(100vh - 431px)", placeholder = "GENE1 0.1\nGENE2 0.05\n..."),
  verbatimTextOutput("count", TRUE),
  numericInput("input.universe", "Genes in universe", 0, 0)
)

homePane <- tabPanel(
  NULL, br(), icon = icon("home", lib = "glyphicon"),
  selectInput("anno.source", "Annotations", c("No annotations" = "")),
  selectInput("anno.types", NULL, c("Gene Sets" = "name", "Symbols" = "syms", "Descriptions" = "info", "Automatic" = "auto", "Manual" = "file"), selected = "file", multiple = TRUE),
  textInput("anno.regex", NULL, placeholder = "Filter using regular expressions"),
  
  selectInput("data.source", "Database", c("No databases" = "")),
  selectInput("data.categories", NULL, c("Filter categories" = ""), multiple = TRUE),
  selectInput("data.organisms", NULL, c("Filter organisms" = ""), multiple = TRUE),
  
  selectInput("cell.source", "Seurat", c("No Seurat data" = "")),
  selectInput("cell.cluster", NULL, c("Inter-cluster" = "_inter_")),
  selectInput("cell.group", NULL, c("No groups" = "")),
  selectInput("cell.select", NULL, c("No clusters" = "")),
  selectInput("cell.compare", NULL, c("No clusters" = "")),
  
  selectInput("expr.source", "Expression", c("No expression data" = "")),
  
  selectInput("conv", "Convert", c("(none)" = "", "Human to mouse" = "hm", "Mouse to human" = "mh"))
)

etcPane <- tabPanel(
  NULL, br(), icon = icon("cog", lib = "glyphicon"),
  radioButtons("info.source", "Description source", c("Annotations" = "anno", "Database" = "data"), inline = TRUE),
  checkboxInput("stat.sigonly", "Significant results only", TRUE),
  checkboxInput("name.fix", "Remove underscores from gene set names", TRUE), hr(),
  
  varSelectizeInput("info.columns", "Information", NULL, multiple = T, options = list(placeholder = "Select columns to display")),
  varSelectizeInput("stat.columns", "Statistics", NULL, multiple = T, options = list(placeholder = "Select columns to display")), hr(),
  
  p("MSigDB XML files available ", a("here", href = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/")),
  fileInput("file.up", NULL, buttonLabel = "Upload..."),
  downloadButton("file.down", "Download", class = "wide"),
  tags$head(tags$style(".wide {width: 100%;}"))
)

plotPane <- tabPanel(
  NULL, br(), icon = icon("stats", lib = "glyphicon"),
  selectInput("bars.value", "Plot value", c("# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), selected = "Fold Enrichment"),
  selectInput("bars.color", "Plot color", c("# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), selected = "Adjusted P-value"),
  selectInput("bars.anno.order", "Plot order", c("Annotation", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment")),
  uiOutput("bars.anno.select"),
  checkboxInput("bars.anno.sort", "Sort plot by value"), hr(),
  
  selectInput("over.color", "Overlap color", c("Gene Value", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment"), "Adjusted P-value"),
  helpText("Selecting Gene Value, Odds Ratio, Fold Enrichment or Adjusted Fold Enrichment options may hide matches"), br(),
  selectInput("over.anno.order", "Annotation order", c("Annotation", "# gene sets", "# genes", "# matches", "P-value", "Adjusted P-value", "Odds Ratio", "Fold Enrichment", "Adjusted Fold Enrichment")),
  uiOutput("over.anno.select"),
  selectInput("over.gene.order", "Gene order", c("As entered" = "input", "Gene name" = "names", "Gene value" = "value")),
  uiOutput("over.gene.select")
)

cellPane <- tabPanel(
  NULL, br(), icon = icon("screenshot", lib = "glyphicon"),
  selectInput("cell.overview", "Overview plot", c("No Seurat data selected" = "")), hr(),
  
  selectInput("cell.plot", "Expression plot", c("Dot Plot" = "dot", "Feature Plot" = "feat", "Heatmap" = "heat", "Ridge Plot" = "ridge", "Violin Plot" = "violin"), "feat"),
  selectInput("cell.anno", "Annotation", c("No annotations loaded" = "")),
  selectInput("cell.genes", "Genes (click to select more)", c("Select genes to display" = ""), multiple = TRUE),
  selectInput("cell.view.cluster", "Clusters", NULL, multiple = TRUE),
  checkboxInput("cell.gene.match", "Restrict to genes in input", TRUE),
  checkboxInput("cell.gene.cluster", "Group genes based on expression")
)

scorePane <- tabPanel(
  NULL, br(), icon = icon("scale", lib = "glyphicon"),
  selectInput("score.type", "Type", c("(none)" = "")),
  selectInput("score.method", "Method", c("Eigen expression" = "pca", "Partial least squares" = "pls", "Partial least squares discriminant analysis" = "plsda")),
  selectInput("score.plot", "Plot", c("Scatter plot" = "whiskers", "Box plot" = "box", "Violin plot" = "violin")),
  selectInput("score.anno", "Annotation", c("No annotations loaded" = "")),
  checkboxInput("score.gene.match", "Restrict to genes in input")
)

barsView <- tabPanel("Plot", icon = icon("stats", lib = "glyphicon"), plotOutput("bars", height = "calc(100vh - 132.5px)"))
overView <- tabPanel("Overlap", icon = icon("equalizer", lib = "glyphicon"), plotOutput("over", height = "calc(100vh - 132.5px)"))
cellView <- tabPanel("Cluster", icon = icon("screenshot", lib = "glyphicon"), plotOutput("cell", height = "calc(100vh - 132.5px)"))
heatView <- tabPanel("Expression", icon = icon("fire", lib = "glyphicon"), plotOutput("heat", height = "calc(100vh - 132.5px)"))
scoreView <- tabPanel("Scores", icon = icon("scale", lib = "glyphicon"), plotOutput("score", height = "calc(100vh - 432.5px)"), plotOutput("rocs", height = 300))
statView <- tabPanel("Statistics", icon = icon("th", lib = "glyphicon"), dataTableOutput("stat"))
contView <- tabPanel("Quality", icon = icon("ok", lib = "glyphicon"), fluidRow(
  column(2, dataTableOutput("cont.anno")), column(2, dataTableOutput("cont.data")),
  column(2, dataTableOutput("cont.cell")), column(3, dataTableOutput("cont.sets")), column(3, dataTableOutput("cont.input"))))
convView <- tabPanel("Convert", icon = icon("transfer", lib = "glyphicon"), dataTableOutput("conv"))
infoView <- tabPanel("Information", icon = icon("info-sign", lib = "glyphicon"), dataTableOutput("info"))

ui <- fluidPage(
  theme = shinytheme("yeti"), useShinyjs(),
  titlePanel("glacier: Gene List Annotation, Calculation and Illustration of Enrichment in R"),
  sidebarLayout(
    sidebarPanel(fluidRow(column(6, inputPane), column(6, tabsetPanel(type = "pills", homePane, etcPane, plotPane, cellPane, scorePane)))),
    mainPanel(tabsetPanel(barsView, overView, cellView, heatView, scoreView, statView, convView, contView, infoView))
  ),
  tags$head(tags$style(".text-compare {height: calc(100vh - 170px); overflow-y: auto;}"))
)
