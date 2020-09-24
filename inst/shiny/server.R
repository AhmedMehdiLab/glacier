library(shiny)
library(shinyjs)
library(shinyWidgets)

library(dplyr)
library(stringr)
library(tools)

# options
. <- NULL
bars_anno_max <- 80
over_anno_max <- 80
over_gene_max <- 100

debounce_time <- 100
datatable_opt <- list(dom = "tr", paging = F, scrollCollapse = T,
                      scrollY = "calc(100vh - 235px)")
upload_file.t <- c("(Select type)" = "", "Annotations" = "anno",
                   "Database" = "data", "MSigDB XML" = "msig",
                   "Seurat RDS" = "cell")
dim_reduction <- c(
  "Principal component analysis" = "pca",
  "Independent component analysis" = "ica",
  "t-distributed Stochastic Neighbor Embedding" = "tsne",
  "Uniform Manifold Approximation and Projection" = "umap"
)

options(shiny.maxRequestSize = 5 * 1024^3, shiny.sanitize.errors = T)
list_show <- function(list)
  capture.output(str(list, comp.str = NULL, give.attr = F, give.head = T,
                     vec.len = 1, list.len = 22))


# upload options
upload_config <- tagList(
  h5("Options"),
  selectInput("file.type", NULL, upload_file.t),
  textInput("file.name", NULL),
  selectInput("file.delim", NULL, c("Commas (CSV)" = ",", "Tabs (TSV)" = "\t")),
  checkboxInput("file.header", "Contains header"),
  sliderInput("file.cols", "Content columns", 1, 10, c(1, 10)),
  numericInput("file.info", "Description column", 0, 10, 1),
  helpText("If no descriptions are available, set 'Description column' to 0"),
  progressBar("file.progress", 0)
)

upload_dialog <- modalDialog(
  title = "Upload file", size = "l", easyClose = T,
  footer = tagList(modalButton("Cancel"), actionButton("file.ok", "Confirm")),
  fluidRow(column(3, upload_config),
           column(9, h5("Preview"), uiOutput("file.view")))
)

server <- function(input, output, session) {
  store <- reactiveValues()
  store$anno <- list()
  store$cell <- list()
  store$data <- list()

  # process upload
  observe(toggleState("file.ok", input$file.type != ""))
  observeEvent(input$file, showModal(upload_dialog))
  observeEvent(input$file.type, {
    path <- input$file$datapath
    name <- file.path_sans_ext(input$file$name)

    # activate controls
    output$file.view <- NULL
    sepr_file <- input$file.type %in% c("anno", "sets")
    sets_file <- input$file.type %in% c("msig", "sets")
    updateTextInput(session, "file.name", value = name)
    updateProgressBar(session, "file.progress", 0)

    toggleState("file.delim", sepr_file)
    toggleState("file.header", sepr_file)
    toggleState("file.cols", sepr_file)
    toggleState("file.info", sepr_file)
    updateProgressBar(session, "file.progress", 10)

    # parse data
    if (sepr_file) temp <- reactive(glacier:::import_delim_path(path, input$file.delim, input$file.header))
    if (input$file.type == "cell") store$proc <- reactive(readRDS(path))
    if (input$file.type == "msig") store$proc <- reactive(import_msigdb(path))
    updateProgressBar(session, "file.progress", 80)

    # update controls
    if (sepr_file) {
      cols <- ncol(temp())
      updateSliderInput(session, "file.cols", value = c(1, cols), min = 1, max = cols)
      updateNumericInput(session, "file.info", value = 0, min = 0, max = cols)
    }

    # prepare preview
    if (input$file.type == "anno") {
      store$proc <- reactive(temp() %>% import_annotations(input$file.cols, input$file.info))
      output$file.view <- renderUI(tagList(
        tableOutput("file.anno"),
        tags$style("#file.anno {height: 500px; overflow: auto}"))
      )
    }
    if (input$file.type == "cell") output$file.view <- renderUI(tagList(
      verbatimTextOutput("file.cell"),
      tags$style("#file.cell {height: 500px; overflow: auto}")
    ))
    if (input$file.type == "sets") store$proc <- reactive(temp() %>% import_database(input$file.cols, input$file.info))
    if (sets_file) output$file.view <- renderUI(tagList(
      column(6, tableOutput("file.sets")),
      column(6, verbatimTextOutput("file.gene")),
      tags$style("#file.sets {height: 500px; overflow: auto}"),
      tags$style("#file.gene {height: 500px; overflow: auto}")
    ))
    updateProgressBar(session, "file.progress", 100)

    # update preview
    output$file.anno <- renderTable(store$proc() %>% head(16) %>% select(all_of(1:2), last_col(1), last_col()))
    output$file.cell <- renderText(str_c(capture.output(store$proc()[[]]), collapse = "\n"))
    output$file.sets <- renderTable(store$proc()$info %>% head(16))
    output$file.gene <- renderText(store$proc()$sets %>% list_show %>% str_c(collapse = "\n"))
  })
  observeEvent(input$file.ok, {
    type <- input$file.type
    if (type == "msig") type <- "data"

    store[[type]][[input$file.name]] <- store$proc()
    store$proc <- NULL
    removeModal()
  })
# 
#   # primary controls
#   init_select <- function(ui, choices) session %>% updateselectInput(ui, choices = choices, selected = choices[1])
#   observe(init_select("anno_from", names(store$items$anno)))
#   observe(init_select("cell_from", names(store$items$cell)))
#   observe(init_select("sets_from", names(store$items$sets)))
#   observe(toggleState("gene_list", input$gene_from == "gene"))
#   observe(toggleState("cell_down", input$cell_type == "heat"))
# 
#   # load data
#   anno_raw <- reactive(store$items$anno[[req(input$anno_from)]])
#   cell_raw <- reactive(store$items$cell[[req(input$cell_from)]])
#   info_raw <- reactive(if (input$anno_info == "anno") anno_raw() %>% select("name", "info") else sets_raw()$info %>% select("name", "info"))
#   sets_raw <- reactive(store$items$sets[[req(input$sets_from)]])
#   universe <- reactive(sets_raw()$sets %>% unlist(use.names = F) %>% unique %>% length %>% max(nrow(gene_proc()$gene_vals)))
# 
#   # secondary controls
#   observe(session %>% updateselectInput("cell_vals", choices = str_c("Cluster ", levels(cell_raw()), " mean values")))
#   observe(session %>% updateNumericInput("gene_univ", value = universe(), min = universe()))
#   observe(init_select("sets_cats", levels(sets_raw()$info$category)))
#   observe(init_select("sets_orgs", levels(sets_raw()$info$organism)))
# 
#   # process data
#   anno_proc <- reactive(anno_raw() %>% process_annotations(info_raw(), input$anno_type))
#   cell_rduc <- reactive(dim_reduction[dim_reduction %in% names(cell_raw()@reductions)])
#   i_gene_list <- reactive(input$gene_list) %>% debounce(debounce_time)
#   gene_pre <- reactive(
#     if (input$gene_from == "gene") i_gene_list() %>% preproc_gene
#     else if (input$gene_from == "cell")
#       cell_raw() %>%
#       GetAssayData %>%
#       as.matrix %>%
#       rowMeans %>%
#       as.data.frame %>%
#       rownames_to_column("gene") %>%
#       rename(value = ".") %>%
#       as_tibble
#   )
#   observe({
#     values <- req(gene_pre()$value)
#     minimum <- values %>% min(na.rm = T) %>% floor
#     maximum <- values %>% max(na.rm = T) %>% ceiling
#     session %>% updateSliderInput("gene_filt", min = minimum, max = maximum, value = c(minimum, maximum))
#   })
#   gene_proc <- reactive(
#     gene_pre() %>% filter(between(value, input$gene_filt[1], input$gene_filt[2])) %>% process_gene
#   )
#   sets_proc <- reactive(sets_raw() %>% process_database(input$sets_cats, input$sets_orgs))
# 
#   # derive data
#   anno_list <- reactive(anno_proc()$annos %>% str_subset(input$anno_filt) %>% as.character)
#   anno_sets <- reactive(anno_proc()$sets_anno$name)
#   gene_list <- reactive(gene_vals()$gene)
#   gene_vals <- reactive(gene_proc()$gene_vals)
#   sets_gene <- reactive(sets_proc()$gene)
#   sets_info <- reactive(sets_proc()$info %>% select("Gene Set" = name, Information = info, Description = any_of("desc"), Category = "category", Organism = "organism"))
#   sets_list <- reactive(sets_proc()$info$name)
#   uniq_anno <- reactive(setdiff(anno_sets(), sets_list()))
#   uniq_sets <- reactive(setdiff(sets_list(), anno_sets()))
# 
#   gene_find_ok <- reactive(intersect(gene_list(), sets_gene()))
#   gene_find_no <- reactive(setdiff(gene_list(), sets_gene()))
#   gene_vals_ok <- reactive(gene_list()[!gene_proc()$no_values])
#   gene_vals_no <- reactive(gene_list()[gene_proc()$no_values])
#   heat_gene <- reactive(annotation_map(input$cell_anno, anno_proc()$sets_anno, sets_proc()$sets_gene)$gene)
# 
#   # primary information
#   output$anno_nums <- renderText(str_c(
#     length(anno_sets()), " gene sets annotated\n",
#     length(anno_list()), " unique annotations"
#   ))
#   output$cell_nums <- renderText(str_c(
#     nrow(cell_raw()), " genes\n",
#     ncol(cell_raw()), " cells\n",
#     length(levels(cell_raw())), " clusters"
#   ))
#   output$gene_nums <- renderText(str_c(
#     length(gene_list()), " unique genes\n",
#     length(gene_find_ok()), " genes recognised\n",
#     length(gene_vals_ok()), " values entered"
#   ))
#   output$sets_nums <- renderText(str_c(
#     length(sets_list()), " gene sets loaded\n",
#     length(sets_gene()), " unique genes"
#   ))
# 
#   output$gene_ok <- renderText(str_c(gene_find_ok(), collapse = ", "))
#   output$gene_no <- renderText(str_c(gene_find_no(), collapse = ", "))
#   output$uniq_anno <- renderText(str_c(uniq_anno(), collapse = ", "))
#   output$uniq_sets <- renderText(str_c(uniq_sets(), collapse = ", "))
#   output$vals_ok <- renderText(str_c(gene_vals_ok(), collapse = ", "))
#   output$vals_no <- renderText(str_c(gene_vals_no(), collapse = ", "))
# 
#   # actions
#   observe(session %>% updateselectInput("cell_anno", choices = anno_list()))
#   observe({
#     genes <- intersect(rownames(cell_raw()), heat_gene())
#     session %>% updateselectInput("cell_gene", choices = genes, selected = if (!is.null(genes) && length(genes)) genes[1])
#   })
#   observe(session %>% updateselectInput("cell_rduc", choices = cell_rduc()))
#   observeEvent(c(uniq_anno(), uniq_sets()), {
#     removeNotification(store$note_over)
#     store$note_over <- if (length(uniq_anno()) + length(uniq_sets()))
#       showNotification("Incomplete overlap between database and annotations", duration = NULL, type = "warning")
#   })
#   observeEvent(gene_find_no(), {
#     removeNotification(store$note_gene)
#     store$note_gene <- if (length(gene_find_no()) && input$gene_from == "gene")
#       showNotification(
#         "Some genes were not recognised", duration = NULL, type = "warning",
#         action = a(href = "https://www.genenames.org/tools/multi-symbol-checker/", "Standardise names here"))
#   })
# 
#   # compute data
#   hits <- reactive(stat_pre()$hits)
#   stat <- reactive(calculate_post(stat_pre()$stat, length(gene_list()), input$gene_univ))
#   stat_pre <- reactive(calculate_pre(gene_list(), anno_list(), anno_proc()$sets_anno, sets_proc()$sets_gene))
# 
#   stat_order <- function(stat, order) stat %>% arrange(if (order == "Annotation") str_to_lower("Annotation") else rev(order))
#   stat_bars <- reactive(stat() %>% stat_order(input$bars_orda))
# 
#   over_anno <- reactive(stat() %>% stat_order(input$over_orda) %>% pull("Annotation"))
#   over_gene <- reactive(
#     if (input$over_ordg == "input") gene_list()
#     else if (input$over_ordg == "names") gene_list() %>% str_sort(numeric = T)
#     else if (input$over_ordg == "value") gene_vals() %>% arrange(value) %>% pull(gene)
#   )
# 
#   # tertiary controls
#   init_range <- function(ui, limit, items) {
#     isolate(if (!ui %in% names(store)) store[[ui]] <- c(0, 0))
#     output[[ui]] <- renderUI(sliderTextInput(ui, NULL, req(items()), selected = items()[c(1, min(length(items()), limit))], force_edges = T))
#     outputOptions(output, ui, suspendWhenHidden = F)
#   }
#   init_range("bars_anno", bars_anno_max, function() stat_bars()$Annotation)
#   init_range("over_anno", over_anno_max, function() req(over_anno()))
#   init_range("over_gene", over_gene_max, function() req(over_gene()))
# 
#   auto_range <- function(ui, limit, items) {
#     lpos <- which(items == input[[ui]][1])
#     rpos <- which(items == input[[ui]][2])
# 
#     store[[ui]] <- if (!length(lpos) || !length(rpos)) c(1, min(length(items), limit))
#     else if (rpos - lpos <= limit) c(lpos, rpos)
#     else if (lpos != store[[ui]][1]) c(lpos, lpos + limit)
#     else if (rpos != store[[ui]][2]) c(rpos - limit, rpos)
#     updateSliderTextInput(session, ui, selected = items[store[[ui]]])
#   }
#   observeEvent(input$bars_anno, auto_range("bars_anno", bars_anno_max, stat_bars()$Annotation))
#   observeEvent(input$over_anno, auto_range("over_anno", over_anno_max, over_anno()))
#   observeEvent(input$over_gene, auto_range("over_gene", over_gene_max, over_gene()))
# 
#   pick_cols <- function(ui, data) updateVarselectInput(session, ui, data = data, selected = names(data))
#   observe(pick_cols("info_cols", sets_info()))
#   observe(pick_cols("stat_cols", stat()))
# 
#   # crop data
#   view_bars <- reactive(stat_bars()[store$bars_anno[1]:store$bars_anno[2], ])
#   view_over_anno <- reactive(over_anno()[store$over_anno[1]:store$over_anno[2]])
#   view_over_gene <- reactive(over_gene()[store$over_gene[1]:store$over_gene[2]])
# 
#   # secondary information
#   cell_gene <- reactive(input$cell_gene) %>% debounce(debounce_time)
#   view_table <- function(table, cols, split)
#     table %>% select(!!!cols) %>% mutate(across(!!split, str_replace_all, "_", ifelse(input$name_long, " ", "_")))
#   output$bars <- renderPlot(display_stats(view_bars(), input$bars_vals, input$bars_colr, input$bars_xfvl, input$bars_xfcl, input$bars_sort))
#   output$over <- renderPlot(display_overlap(view_over_gene(), view_over_anno(), hits(), input$over_xfcl))
#   output$cell <- renderPlot(DimPlot(cell_raw(), reduction = input$cell_rduc, label = T, repel = T))
#   output$heat <- renderPlot({
#     feats <- req(cell_gene())
#     nfeat <- feats %>% length %>% sqrt %>% ceiling
#     gene <- if (!input$cell_clst) feats else
#       FindAllMarkers(cell_raw(), features = feats, only.pos = T) %>% group_by(cluster) %>% pull(gene)
# 
#     if (input$cell_type == "dot") {
#       DotPlot(cell_raw(), features = gene)
#     } else if (input$cell_type == "feat") {
#       FeaturePlot(cell_raw(), features = gene, ncol = nfeat)
#     } else if (input$cell_type == "heat") {
#       DoHeatmap(subset(cell_raw(), downsample = input$cell_down), features = gene)
#     } else if (input$cell_type == "ridge") {
#       RidgePlot(cell_raw(), features = gene, ncol = nfeat)
#     } else if (input$cell_type == "violin") {
#       VlnPlot(cell_raw(), features = gene, ncol = nfeat)
#     }
#   })
#   output$stat <- renderDataTable(view_table(stat(), input$stat_cols, "Annotation"), datatable_opt)
#   output$info <- renderDataTable(view_table(sets_info(), input$info_cols, "Gene Set"), datatable_opt)
# 
#   # downloads
#   output$down_load <- downloadHandler("glacier_results.csv", . %>% write_csv(stat(), .))
}
