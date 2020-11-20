library(glacier)
requireNamespace("biomaRt")
source("score.R")
source("upload.R")

library(dplyr)
library(magrittr)
library(purrr)
library(readr)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(stringr)
library(tibble)
library(tools)

EXAMPLE <- TRUE
PRIVATE <- TRUE
SHINYIO <- FALSE

if (SHINYIO) {
  # devtools::install_github("AhmedMehdiLab/glacier")
  # setRepositories()
  
  library(limma)
  library(Seurat)
}

# options
BARS_ANNO_MAX <- 60
OVER_ANNO_MAX <- 60
OVER_GENE_MAX <- 60
DEBOUNCE_TIME <- 1500

SEED <- 444
DT_LARGE <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 235px)")
DT_SMALL <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 220.3px)")
DIMREDUC <- c("Principal component analysis" = "pca",
              "Independent component analysis" = "ica",
              "t-distributed Stochastic Neighbor Embedding" = "tsne",
              "Uniform Manifold Approximation and Projection" = "umap")
options(shiny.maxRequestSize = 5 * 1024 ^ 3) # 5 GiB max upload size

server <- function(input, output, session) {
  store <- reactiveValues(anno = list(), cell = list(), data = list(), mart = list(), note = list(), proc = reactiveValues(), time = list())
  showNotification(str_c("Welcome to glacier (", packageVersion("glacier"), ")!", collapse = ""), type = "message")
  
  if (PRIVATE) isolate({
    store$anno$`E.PAGE` <- readRDS(system.file("extdata", "private", "epage_anno.rds", package = "glacier"))
    store$anno$`MSigDB C7` <- readRDS(system.file("extdata", "private", "msigdb_anno_orig.rds", package = "glacier"))
    store$anno$`MSigDB C7 (MODIFIED)` <- readRDS(system.file("extdata", "private", "msigdb_anno.rds", package = "glacier"))
    
    store$data$`E.PAGE` <- readRDS(system.file("extdata", "private", "epage_data.rds", package = "glacier"))
    store$data$`MSigDB 7.2` <- readRDS(system.file("extdata", "private", "msigdb_data.rds", package = "glacier"))
    
    store$cell$GSE150728_nano <- readRDS(system.file("extdata", "private", "GSE150728_cell.rds", package = "glacier"))
    store$time$GSE150728 <- readRDS(system.file("extdata", "private", "GSE150728_time.rds", package = "glacier"))
  })
  
  if (SHINYIO) isolate({
    store$anno$`E.PAGE` <- readRDS("EPAGE_A.rds")
    store$anno$`MSigDB C7` <- readRDS("MSIG_A.rds")
    store$data$`E.PAGE` <- readRDS("EPAGE_D.rds")
  })
  
  if (EXAMPLE) isolate({
    store$anno$Example <- import_annotations(system.file("extdata", "ex_anno.csv", package = "glacier"), ",", TRUE, c(2, 4), 5)
    store$data$Example <- import_database(system.file("extdata", "ex_data.csv", package = "glacier"), ",", FALSE, c(2, 4), 0)
    store$cell$Example <- readRDS(system.file("extdata", "ex_seurat.rds", package = "glacier"))
  })
  
  # process upload
  observeEvent(input$file.up, showModal(uploadUI("file")))
  uploadServer("file", reactive(input$file.up), store$proc)
  observeEvent(input$confirm, {store[[store$proc$type]][[store$proc$name]] <- store$proc$proc(); removeModal()})
  
  # primary controls
  observe(updateSelectInput(session, "anno.source", NULL, names(store$anno)))
  observe(updateSelectInput(session, "cell.source", NULL, names(store$cell)))
  observe(updateSelectInput(session, "data.source", NULL, names(store$data)))
  observe(updateSelectInput(session, "time.source", NULL, names(store$time)))
  observe(toggleState("input.source", input$cell.source != ""))
  observe(toggleState("input.text", input$input.source == "text"))
  observe(if (input$input.source == "cell") updateTextAreaInput(session, "input.text", value = str_c(input_proc()$gene, input_proc()$value, sep = " \t", collapse = "\n")))
  
  # load data
  anno_raw <- reactive(store$anno[[req(input$anno.source)]])
  cell_raw <- reactive(if (requireNamespace("Seurat", quietly = T)) store$cell[[req(input$cell.source)]])
  data_raw <- reactive(store$data[[req(input$data.source)]])
  info_raw <- reactive(if (input$info.source == "anno") anno_raw() else data_raw()$gs_info)
  time_raw <- reactive(store$time[[req(input$time.source)]])
  universe <- reactive(data_raw()$gs_genes %>% unlist(use.names = F) %>% c(input_proc()$gene) %>% unique %>% length)
  
  cell_group_name <- reactive(if ("grp" %in% names(cell_raw()@meta.data)) "grp" else "group")
  cell_clust_list <- reactive(levels(cell_raw()))
  cell_group_list <- reactive(unique(cell_raw()@meta.data %>% pull(cell_group_name())))
  info <- reactive(info_raw() %>% select("name", "info"))
  
  # secondary controls
  observe({
    opts <- tryCatch(
      if (is.null(cell_group_list())) c("Inter-cluster" = "_inter_")
      else c("Inter-cluster" = "_inter_", "Intra-cluster" = "_intra_", setNames(cell_clust_list(), str_c("Intra-cluster ", cell_clust_list()))),
      error = function(e) c("Inter-cluster" = "_inter_")
    )
    
    updateSelectInput(session, "cell.cluster", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- setNames(cell_clust_list(), str_c("Compare cluster ", cell_clust_list()))
      else opts <- setNames(cell_group_list(), str_c("Compare group ", cell_group_list())),
      error = function(e) NULL
    )
    
    updateSelectInput(session, "cell.select", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- c("Against all others" = "_all_", setNames(cell_clust_list(), str_c("Against cluster ", cell_clust_list())))
      else opts <- c("Against all others" = "_all_", setNames(cell_group_list(), str_c("Against group ", cell_group_list()))),
      error = function(e) NULL
    )
    
    updateSelectInput(session, "cell.compare", NULL, opts)
  })
  observe(updateSelectInput(session, "cell.view.cluster", NULL, cell_clust_list(), cell_clust_list()))
  observe(updateSelectInput(session, "data.categories", NULL, levels(data_raw()$gs_info$category), levels(data_raw()$gs_info$category)[1]))
  observe(updateSelectInput(session, "data.organisms", NULL, levels(data_raw()$gs_info$organism), levels(data_raw()$gs_info$organism)[1]))
  observe(updateNumericInput(session, "input.universe", NULL, universe(), universe()))
  
  # process data
  anno_proc <- reactive(glacier:::process_annotations(anno_raw(), info(), input$anno.types))
  cell_proc <- reactive(tryCatch(
    process_input_seurat(cell_raw(), input$cell.select, if (input$cell.compare != "_all_") input$cell.compare, if (input$cell.cluster != "_inter_") "grp", if (!input$cell.cluster %in% c("_inter_", "_intra_")) input$cell.cluster),
    error = function(e) tibble(gene = character(), value = numeric())
  ))
  data_proc <- reactive(glacier:::process_database(data_raw(), input$data.categories, input$data.organisms))
  text_proc <- reactive(input$input.text %>% process_input_text) %>% debounce(DEBOUNCE_TIME)
  cell_reductions <- reactive(DIMREDUC[DIMREDUC %in% names(cell_raw()@reductions)])
  input_proc <- reactive(if (input$input.source == "text") text_proc() else cell_proc())
  
  # derive data
  anno_list <- reactive(tryCatch(anno_proc()$annos %>% str_subset(input$anno.regex) %>% as.character, error = function(e) character()))
  cell_gene <- reactive(rownames(cell_raw()))
  anno_sets <- reactive(anno_proc()$gs_anno$name)
  data_list <- reactive(data_proc()$genes)
  data_info <- reactive(data_proc()$gs_info %>% select(`Gene Set` = name, Information = info, Description = any_of("desc"), Category = category, Organism = organism))
  data_sets <- reactive(data_proc()$gs_info$name)
  
  cell_gene_list <- reactive(if (input$cell.gene.match) intersect(input_proc()$gene, cell_gene()) else cell_gene())
  cell_anno_list <- reactive(calc_filt() %>% filter(`Adjusted P-value` <= 0.05) %>% pull(Annotation) %>% str_sort(numeric = TRUE) %>% setNames(., nm = .) %>% map(~glacier:::explore_annotation(., anno_proc()$gs_annos, data_proc()$gs_genes, cell_gene_list())$genes) %>% compact())
  cell_anno_gene <- reactive(if (input$cell.anno == "") character() else cell_anno_list()[[input$cell.anno]] %>% str_sort(numeric = TRUE))
  
  score_gene_list <- reactive(if (input$score.gene.match) intersect(input_proc()$gene, cell_gene()) else cell_gene())
  score_anno_list <- reactive(calc_filt() %>% filter(`Adjusted P-value` <= 0.05) %>% pull(Annotation) %>% str_sort(numeric = TRUE) %>% setNames(., nm = .) %>% map(~glacier:::explore_annotation(., anno_proc()$gs_annos, data_proc()$gs_genes, score_gene_list())$genes) %>% compact())
  score_anno_gene <- reactive(if (input$score.anno == "") character() else score_anno_list()[[input$score.anno]])
  
  cont_sets <- reactive({
    anno_part <- anno_sets() %>% tibble(name = .) %>% add_column(anno = TRUE)
    data_part <- data_sets() %>% tibble(name = .) %>% add_column(data = TRUE)
    results <- anno_part %>% full_join(data_part) %>% transmute("Set" = str_replace_all(name, "_", ifelse(input$name.fix, " ", "_")), "Status" = ifelse(!is.na(anno) & !is.na(data), "OK", ifelse(is.na(anno), "Database only", "Annotations only")))
    
    if (all(results$Status == "OK")) removeNotification(store$note$overlap)
    else if (is.null(store$note$overlap)) store$note$overlap <- showNotification("Gene set mismatch between annotations and database, check 'Quality' tab", duration = NULL, type = "warning")
    return(results)
  })
  cont_input <- reactive({
    results <- tryCatch(
      input_proc() %>% rename(Input = "gene", Value = "value") %>% mutate(Recognised = Input %in% data_list()),
      error = function(e) tibble(Input = character(), Value = numeric(), Recognised = logical())
    )
    
    if (all(results$Recognised)) removeNotification(store$note$recognised)
    else if (is.null(store$note$recognised)) store$note$recognised <- showNotification("Some genes were not recognised, check 'Quality' tab", duration = NULL, type = "warning")
    
    if (!any(is.na(results$Value)) || all(is.na(results$Value))) removeNotification(store$note$values)
    else if (is.null(store$note$values)) store$note$values <- showNotification("Some genes do not have values, check 'Quality' tab", duration = NULL, type = "warning")
    return(results)
  })
  
  # primary information
  output$anno.count <- renderText(str_c(length(anno_sets()), " gene sets annotated\n", length(anno_list()), " unique annotations"))
  output$cell.count <- renderText(str_c(nrow(cell_raw()), " genes\n", ncol(cell_raw()), " cells\n", length(cell_clust_list()), " clusters"))
  output$data.count <- renderText(str_c(length(data_sets()), " gene sets loaded\n", length(data_list()), " unique genes"))
  output$input.count <- renderText(str_c(nrow(input_proc()), " unique genes\n", sum(cont_input()$Recognised), " genes recognised\n", sum(!is.na(cont_input()$Value)), " values entered"))
  
  output$cont.anno <- renderDataTable(calc()$stats %>% select("Annotation"), DT_SMALL)
  output$cont.sets <- renderDataTable(cont_sets(), DT_SMALL)
  output$cont.cell <- renderDataTable(tibble(Seurat = rownames(cell_raw())), DT_SMALL)
  output$cont.data <- renderDataTable(tibble(Database = data_list()), DT_SMALL)
  output$cont.input <- renderDataTable(cont_input(), DT_SMALL)

  # actions
  observe(updateSelectInput(session, "cell.anno", NULL, names(cell_anno_list())))
  observe(updateSelectInput(session, "cell.genes", NULL, cell_anno_gene(), head(cell_anno_gene(), 16)))
  observe(updateSelectInput(session, "cell.overview", NULL, cell_reductions()))
  observe(updateSelectInput(session, "score.anno", NULL, names(cell_anno_list())))
  
  # compute data
  calc <- reactive({
    withProgress(message = "Calculating statistics", {
      setProgress(value = 0, detail = "Finding matches")
      pre <- glacier:::calculate_pre(input_proc(), anno_list(), anno_proc()$gs_annos, data_proc()$gs_genes)
      
      setProgress(value = 0.8, detail = "Computing significance")
      post <- glacier:::calculate_post(pre$stats_pre, nrow(input_proc()), max(input$input.universe, universe()))
      
      list(stats = post, matches = pre$matches)
    })
  })
  calc_filt <- reactive({
    stats <- calc()$stats
    
    if (input$stat.sigonly) {
      stats <- stats %>% filter(`Adjusted P-value` <= 0.05)
    }
    
    return(stats)
  })
  scores_stat <- reactive({
    if (!requireNamespace("mixOmics", quietly = T) ||
        !requireNamespace("pROC", quietly = T)) {
      return(NULL)
    }
    
    if (input$score.type == "time" && input$time.source != "") func <- score_exp
    else func <- score_seurat
    
    withProgress(message = "Calculating scores", {
      if (input$score.type == "time") {
        if (input$time.source != "") {
          time_raw() %>% mutate(bin = cut_interval(timeM, 20)) %>% score_exp(intersect(score_anno_gene(), colnames(time_raw())))
        }
      } else {
        cell_raw() %>% score_seurat(score_anno_gene())
      }
    })
  })
  
  bars_stat <- reactive(calc_filt() %>% arrange(
    if (input$bars.anno.order == "Annotation") str_to_lower(Annotation)
    else if (input$bars.anno.order %in% c("P-value", "Adjusted P-value")) .data[[input$bars.anno.order]]
    else desc(.data[[input$bars.anno.order]])
  ))
  over_stat <- reactive(calc_filt() %>% arrange(
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
  observe(updateVarSelectInput(session, "stat.columns", NULL, calc_filt(), names(calc_filt())))
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
  output$over <- renderPlot(plot_overlap(calc()$matches, input$over.color, view_over_gene(), view_over(), input$over.color.trans))
  output$cell <- renderPlot(if (requireNamespace("Seurat", quietly = T)) Seurat::DimPlot(cell_raw(), reduction = input$cell.overview, label = T, repel = T))
  output$heat <- renderPlot({
    if (!requireNamespace("Seurat", quietly = T) || !length(view_cell_gene())) return(NULL)
    
    feats <- if (input$cell.gene.match && input$cell.gene.cluster) view_cell_gene_clust() else view_cell_gene()
    if (input$cell.plot != "heat") feats <- unique(feats)
    
    width <- feats %>% length %>% sqrt %>% ceiling
    sample <- subset(cell_raw(), downsample = input$cell.downsample, idents = view_cell_clust(), seed = SEED)
    if (input$cell.plot == "dot") Seurat::DotPlot(sample, features = feats)
    else if (input$cell.plot == "feat") Seurat::FeaturePlot(sample, features = feats, ncol = width, label = TRUE, repel = TRUE)
    else if (input$cell.plot == "heat") Seurat::DoHeatmap(sample, features = feats)
    else if (input$cell.plot == "ridge") Seurat::RidgePlot(sample, features = feats, ncol = width)
    else if (input$cell.plot == "violin") Seurat::VlnPlot(sample, features = feats, ncol = width)
  })
  output$score <- renderPlot({
    if (input$score.type == "time") {
      if (input$time.source != "") {
        scores_stat() %>%
          pluck("scores") %>%
          score_summary(input$score.style, c("grp", "bin")) %>%
          show_summary(input$score.anno, "bin")
      }
    } else {
      scores_stat() %>%
        pluck("scores") %>%
        score_summary(input$score.style, c("grp", "seurat_clusters")) %>%
        show_summary(input$score.anno, "seurat_clusters")
    }
  })
  output$rocs <- renderPlot({
    if (input$score.type == "time") {
      if (input$time.source != "") {
        scores_stat() %>%
          pluck("rocs") %>%
          add_column(cluster = 0) %>%
          show_rocs(input$score.style)
      }
    } else {
      scores_stat() %>%
        pluck("rocs") %>%
        show_rocs(input$score.style)
    }
  })
  output$trans.out <- renderDataTable({
    withProgress(message = "Connecting to Ensembl", {
      setProgress(value = 0.1, detail = "Retrieving Homo sapiens data")
      while (is.null(store$mart$hs)) store$mart$hs <- tryCatch(biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl"), error = function(e) {print("Retrying"); NULL})
      
      setProgress(value = 0.5, detail = "Retrieving Mus musculus data")
      while (is.null(store$mart$mm)) store$mart$mm <- tryCatch(biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl"), error = function(e) {print("Retrying"); NULL})
    })
    
    withProgress(message = "Converting genes", {
      genes <- ""
      vals <- input_proc()$gene
      if (length(vals) == 0) return(NULL)
      if (input$trans == "mh") genes <- biomaRt::getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = vals, mart = store$mart$mm, attributesL = "hgnc_symbol", martL = store$mart$hs)
      if (input$trans == "hm") genes <- biomaRt::getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = vals, mart = store$mart$hs, attributesL = "mgi_symbol", martL = store$mart$mm)
    })
    
    return(genes)
  }, DT_LARGE)
  output$stat <- renderDataTable(view_table(calc_filt(), input$stat.columns, "Annotation"), DT_LARGE)
  output$info <- renderDataTable(view_table(data_info(), input$info.columns, "Gene Set"), DT_LARGE)
  
  output$file.down <- downloadHandler("glacier_results.csv", . %>% write_csv(calc_filt(), .))
}
