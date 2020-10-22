library(shiny)
library(shinyjs)
library(shinyWidgets)

library(glacier)

library(dplyr)
library(magrittr)
library(purrr)
library(readr)
library(stringr)
library(tools)

# if pushing to shinyapps.io
# devtools::install_github("AhmedMehdiLab/glacier")
# setRepositories()
# library(limma)
# library(Seurat)

source("upload.R")

# options
LOAD_EXAMPLES <- TRUE
BARS_ANNO_MAX <- 80
OVER_ANNO_MAX <- 80
OVER_GENE_MAX <- 100
DEBOUNCE_TIME <- 1000

RAND_SEED <- 444
LARGE_DT <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 235px)")
SMALL_DT <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 220.3px)")
DIM_RED <- c("Principal component analysis" = "pca", "Independent component analysis" = "ica", "t-distributed Stochastic Neighbor Embedding" = "tsne", "Uniform Manifold Approximation and Projection" = "umap")
options(shiny.maxRequestSize = 5 * 1024 ^ 3)

server <- function(input, output, session) {
  store <- reactiveValues(anno = list(), cell = list(), data = list(), proc = reactiveValues())
  showNotification(str_c("Welcome to glacier (v", packageVersion("glacier"), ")!", collapse = ""), type = "message")
  
  if (LOAD_EXAMPLES) isolate({
    # if pushing to shinyapps.io
    # store$anno$GEPATH <- readRDS("GEA.rds")
    # store$anno[["MSigDB C7 example"]] <- readRDS("MSA.rds")
    # store$data$GEPATH <- readRDS("GED.rds")
    
    # load private data
    store$anno$`E.PAGE` <- readRDS(system.file("extdata", "private", "epage_anno.rds", package = "glacier"))
    store$anno$`MSigDB C7` <- readRDS(system.file("extdata", "private", "msigdb_anno_orig.rds", package = "glacier"))
    store$anno$`MSigDB C7 <MODIFIED>` <- readRDS(system.file("extdata", "private", "msigdb_anno.rds", package = "glacier"))
    store$data$`E.PAGE` <- readRDS(system.file("extdata", "private", "epage_data.rds", package = "glacier"))
    store$data$`MSigDB 7.2` <- readRDS(system.file("extdata", "private", "msigdb_data.rds", package = "glacier"))
    
    # store$anno$Example <- import_annotations(system.file("extdata", "ex_anno.csv", package = "glacier"), ",", TRUE, c(2, 4), 5)
    # store$data$Example <- import_database(system.file("extdata", "ex_data.csv", package = "glacier"), ",", FALSE, c(2, 4), 0)
    store$cell$Example <- readRDS(system.file("extdata", "ex_seurat.rds", package = "glacier"))
  })
  
  # process upload
  observeEvent(input$file.up, showModal(uploadUI("file")))
  uploadServer("file", reactive(input$file.up), store$proc)
  observeEvent(input$confirm, {
    store[[store$proc$type]][[store$proc$name]] <- store$proc$proc()
    removeModal()
  })
  
  # primary controls
  observe(updateSelectInput(session, "anno.source", NULL, names(store$anno)))
  observe(updateSelectInput(session, "cell.source", NULL, names(store$cell)))
  observe(updateSelectInput(session, "data.source", NULL, names(store$data)))
  observe(toggleState("input.source", input$cell.source != ""))
  observe(toggleState("input.text", input$input.source == "text"))
  observe(if (input$input.source == "cell") updateTextAreaInput(session, "input.text", value = str_c(input_proc()$gene, input_proc()$value, sep = " \t", collapse = "\n")))
  
  # load data
  anno_raw <- reactive(store$anno[[req(input$anno.source)]])
  cell_raw <- reactive(if (requireNamespace("Seurat", quietly = T)) store$cell[[req(input$cell.source)]])
  data_raw <- reactive(store$data[[req(input$data.source)]])
  info_raw <- reactive(if (input$info.source == "anno") anno_raw() else data_raw()$gs_info)
  universe <- reactive(data_raw()$gs_genes %>% unlist(use.names = F) %>% c(input_proc()$gene) %>% unique %>% length)
  
  cgname <- reactive(if ("grp" %in% names(cell_raw()@meta.data)) "grp" else "group")
  clusts <- reactive(levels(cell_raw()))
  cgroup <- reactive(unique(cell_raw()@meta.data %>% pull(cgname())))
  info <- reactive(info_raw() %>% select("name", "info"))
  
  # secondary controls
  observe({
    opts <- tryCatch(
      if (is.null(cgroup())) c("Inter-cluster" = "_inter_")
      else c("Inter-cluster" = "_inter_", "Intra-cluster" = "_intra_", setNames(clusts(), str_c("Intra-cluster ", clusts()))),
      error = function(e) c("Inter-cluster" = "_inter_")
    )
    
    updateSelectInput(session, "cell.cluster", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- setNames(clusts(), str_c("Compare cluster ", clusts()))
      else opts <- setNames(cgroup(), str_c("Compare group ", cgroup())),
      error = function(e) NULL
    )
    
    updateSelectInput(session, "cell.select", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- c("Against all others" = "_all_", setNames(clusts(), str_c("Against cluster ", clusts())))
      else opts <- c("Against all others" = "_all_", setNames(cgroup(), str_c("Against group ", cgroup()))),
      error = function(e) NULL
    )
    
    updateSelectInput(session, "cell.compare", NULL, opts)
  })
  observe(updateSelectInput(session, "cell.view.cluster", NULL, clusts(), clusts()))
  observe(updateSelectInput(session, "data.categories", NULL, levels(data_raw()$gs_info$category), levels(data_raw()$gs_info$category)[1]))
  observe(updateSelectInput(session, "data.organisms", NULL, levels(data_raw()$gs_info$organism), levels(data_raw()$gs_info$organism)[1]))
  observe(updateNumericInput(session, "input.universe", NULL, universe(), universe()))
  
  # process data
  anno_proc <- reactive(glacier:::process_annotations(anno_raw(), info(), input$anno.types))
  cell_proc <- reactive({
    tryCatch(
      process_input_seurat(cell_raw(), input$cell.select, if (input$cell.compare != "_all_") input$cell.compare, if (input$cell.cluster != "_inter_") "grp", if (!input$cell.cluster %in% c("_inter_", "_intra_")) input$cell.cluster),
      error = function(e) tibble(gene = character(), value = numeric())
    )
  })
  data_proc <- reactive(glacier:::process_database(data_raw(), input$data.categories, input$data.organisms))
  text_proc <- reactive(input$input.text %>% process_input_text) %>% debounce(DEBOUNCE_TIME)
  cell_reductions <- reactive(DIM_RED[DIM_RED %in% names(cell_raw()@reductions)])
  input_proc <- reactive(if (input$input.source == "text") text_proc() else cell_proc())
  
  # derive data
  anno_list <- reactive(tryCatch(anno_proc()$annos %>% str_subset(input$anno.regex) %>% as.character, error = function(e) character()))
  anno_only <- reactive(setdiff(anno_sets(), data_sets()))
  cell_gene <- reactive(rownames(cell_raw()))
  anno_sets <- reactive(anno_proc()$gs_anno$name)
  data_list <- reactive(data_proc()$genes)
  data_info <- reactive(data_proc()$gs_info %>% select(`Gene Set` = name, Information = info, Description = any_of("desc"), Category = category, Organism = organism))
  data_only <- reactive(setdiff(data_sets(), anno_sets()))
  data_sets <- reactive(data_proc()$gs_info$name)
  
  cell_gene_list <- reactive(if (input$cell.gene.match) intersect(input_proc()$gene, cell_gene()) else cell_gene())
  cell_anno_list <- reactive(calc_post() %>% filter(`Adjusted P-value` <= 0.05) %>% pull(Annotation) %>% str_sort(numeric = TRUE) %>% setNames(., nm = .) %>% map(~glacier:::explore_annotation(., anno_proc()$gs_annos, data_proc()$gs_genes, cell_gene_list())$genes) %>% compact())
  cell_anno_gene <- reactive(if (input$cell.anno == "") character() else cell_anno_list()[[input$cell.anno]] %>% str_sort(numeric = TRUE))
  
  gene_ok <- reactive(intersect(input_proc()$gene, data_list()))
  gene_no <- reactive(setdiff(input_proc()$gene, data_list()))
  vals_ok <- reactive(input_proc() %>% filter(!is.na(value), .preserve = T) %>% pull(gene))
  vals_no <- reactive(input_proc() %>% filter(is.na(value), .preserve = T) %>% pull(gene))
  
  # primary information
  output$anno.count <- renderText(str_c(length(anno_sets()), " gene sets annotated\n", length(anno_list()), " unique annotations"))
  output$cell.count <- renderText(str_c(nrow(cell_raw()), " genes\n", ncol(cell_raw()), " cells\n", length(clusts()), " clusters"))
  output$data.count <- renderText(str_c(length(data_sets()), " gene sets loaded\n", length(data_list()), " unique genes"))
  output$input.count <- renderText(str_c(nrow(input_proc()), " unique genes\n", length(gene_ok()), " genes recognised\n", length(vals_ok()), " values entered"))
  
  output$cell.gene <- renderDataTable(tibble(`Genes in Seurat file` = rownames(cell_raw())), SMALL_DT)
  output$data.gene <- renderDataTable(tibble(`Genes in database` = data_list()), SMALL_DT)
  output$anno.only <- renderDataTable(tibble(`Gene sets only in annotations` = anno_only()), SMALL_DT)
  output$data.only <- renderDataTable(tibble(`Gene sets only in database` = data_only()), SMALL_DT)
  output$gene.ok <- renderDataTable(tibble(`Genes recognised` = gene_ok()), SMALL_DT)
  output$gene.no <- renderDataTable(tibble(`Genes not recognised` = gene_no()), SMALL_DT)
  output$vals.ok <- renderDataTable(tibble(`Genes with values` = vals_ok()), SMALL_DT)
  output$vals.no <- renderDataTable(tibble(`Genes without values` = vals_no()), SMALL_DT)
  
  # actions
  observe(updateSelectInput(session, "cell.anno", NULL, names(cell_anno_list())))
  observe(updateSelectInput(session, "cell.genes", NULL, cell_anno_gene(), head(cell_anno_gene(), 16)))
  observe(updateSelectInput(session, "cell.overview", NULL, cell_reductions()))
  observeEvent(c(anno_only(), data_only()), {
    removeNotification(store$note.overlap)
    store$note.overlap <- if (length(anno_only()) + length(data_only())) showNotification("Gene sets of selected annotations and database do not match", duration = NULL, type = "warning")
  })
  
  # compute data
  matches <- reactive(calc_pre()$matches)
  calc_pre <- reactive(glacier:::calculate_pre(input_proc(), anno_list(), anno_proc()$gs_annos, data_proc()$gs_genes))
  calc_post <- reactive({
    stats <- glacier:::calculate_post(calc_pre()$stats_pre, nrow(input_proc()), max(input$input.universe, universe()))
    
    if (input$stat.sigonly) {
      stats <- stats %>% filter(`Adjusted P-value` <= 0.05)
    }
    
    return(stats)
  })
  
  bars_stat <- reactive(calc_post() %>% arrange(
    if (input$bars.anno.order == "Annotation") str_to_lower(Annotation)
    else if (input$bars.anno.order %in% c("P-value", "Adjusted P-value")) .data[[input$bars.anno.order]]
    else desc(.data[[input$bars.anno.order]])
  ))
  over_stat <- reactive(calc_post() %>% arrange(
    if (input$over.anno.order == "Annotation") str_to_lower(Annotation)
    else if (input$over.anno.order %in% c("P-value", "Adjusted P-value")) .data[[input$over.anno.order]]
    else desc(.data[[input$over.anno.order]])
  ))
  over_input <- reactive(
    if (input$over.gene.order == "input") input_proc()
    else if (input$over.gene.order == "names") input_proc() %>% arrange(gene)
    else if (input$over.gene.order == "value") input_proc() %>% arrange(value)
  )
  
  # tertiary controls
  init_text_range <- function(ui, limit, items) {
    isolate(if (!ui %in% names(store)) store[[ui]] <- c(0, 0))
    output[[ui]] <- renderUI(sliderTextInput(ui, NULL, req(items()), selected = items()[c(1, min(length(items()), limit))], force_edges = T))
    outputOptions(output, ui, suspendWhenHidden = F)
  }
  init_text_range("bars.anno.select", BARS_ANNO_MAX, function() bars_stat()$Annotation)
  init_text_range("over.anno.select", OVER_ANNO_MAX, function() over_stat()$Annotation)
  init_text_range("over.gene.select", OVER_GENE_MAX, function() over_input()$gene)
  bars_anno_select <- reactive(input$bars.anno.select) %>% debounce(DEBOUNCE_TIME)
  over_anno_select <- reactive(input$over.anno.select) %>% debounce(DEBOUNCE_TIME)
  over_gene_select <- reactive(input$over.gene.select) %>% debounce(DEBOUNCE_TIME)
  
  auto_text_range <- function(ui, limit, items) {
    lpos <- which(items == input[[ui]][1])
    rpos <- which(items == input[[ui]][2])
    
    store[[ui]] <- if (!length(lpos) || !length(rpos)) c(1, min(length(items), limit + 1))
    else if (rpos - lpos <= limit) c(lpos, rpos)
    else if (lpos != store[[ui]][1]) c(lpos, lpos + limit)
    else if (rpos != store[[ui]][2]) c(rpos - limit, rpos)
    updateSliderTextInput(session, ui, NULL, items[store[[ui]]])
  }
  observeEvent(bars_anno_select(), auto_text_range("bars.anno.select", BARS_ANNO_MAX, bars_stat()$Annotation))
  observeEvent(over_anno_select(), auto_text_range("over.anno.select", OVER_ANNO_MAX, over_stat()$Annotation))
  observeEvent(over_gene_select(), auto_text_range("over.gene.select", OVER_GENE_MAX, over_input()$gene))
  observe(updateVarSelectInput(session, "info.columns", NULL, data_info(), names(data_info())))
  observe(updateVarSelectInput(session, "stat.columns", NULL, calc_post(), names(calc_post())))
  observe(toggleState("cell.gene.cluster", input$cell.gene.match))
  
  # crop data
  view_bars <- reactive(bars_stat()[store$bars.anno.select[1]:store$bars.anno.select[2], ])
  view_over <- reactive(over_stat()[store$over.anno.select[1]:store$over.anno.select[2], ])
  view_over_gene <- reactive(over_input()[store$over.gene.select[1]:store$over.gene.select[2], ])
  view_cell_gene <- reactive(input$cell.genes) %>% debounce(DEBOUNCE_TIME)
  view_cell_gene_clust <- reactive(cell_raw() %>% Seurat::FindAllMarkers(features = view_cell_gene()) %>% group_by(cluster) %>% pull(gene))
  view_cell_clust <- reactive(input$cell.view.cluster) %>% debounce(DEBOUNCE_TIME)
  
  # secondary information
  view_table <- function(table, cols, split) table %>% select(!!!cols) %>% mutate(across(!!split, str_replace_all, "_", ifelse(input$name.fix, " ", "_")))
  output$bars <- renderPlot(plot_stats(view_bars(), input$bars.value, input$bars.color, input$bars.value.trans, input$bars.color.trans, input$bars.anno.sort))
  output$over <- renderPlot(plot_overlap(calc_pre()$matches, input$over.color, view_over_gene(), view_over(), input$over.color.trans))
  output$cell <- renderPlot(if (requireNamespace("Seurat", quietly = T)) Seurat::DimPlot(cell_raw(), reduction = input$cell.overview, label = T, repel = T))
  output$heat <- renderPlot({
    if (!requireNamespace("Seurat", quietly = T) || !length(view_cell_gene())) return(NULL)
    
    feats <- if (input$cell.gene.match && input$cell.gene.cluster) view_cell_gene_clust() else view_cell_gene()
    if (input$cell.plot != "heat") feats <- unique(feats)
    
    width <- feats %>% length %>% sqrt %>% ceiling
    sample <- subset(cell_raw(), downsample = input$cell.downsample, idents = view_cell_clust(), seed = RAND_SEED)
    if (input$cell.plot == "dot") Seurat::DotPlot(sample, features = feats)
    else if (input$cell.plot == "feat") Seurat::FeaturePlot(sample, features = feats, ncol = width, label = TRUE, repel = TRUE)
    else if (input$cell.plot == "heat") Seurat::DoHeatmap(sample, features = feats)
    else if (input$cell.plot == "ridge") Seurat::RidgePlot(sample, features = feats, ncol = width)
    else if (input$cell.plot == "violin") Seurat::VlnPlot(sample, features = feats, ncol = width)
  })
  output$stat <- renderDataTable(view_table(calc_post(), input$stat.columns, "Annotation"), LARGE_DT)
  output$info <- renderDataTable(view_table(data_info(), input$info.columns, "Gene Set"), LARGE_DT)
  
  output$file.down <- downloadHandler("glacier_results.csv", . %>% write_csv(calc_post(), .))
}
