library(glacier)
requireNamespace("biomaRt")
source("upload.R")

library(dplyr)
library(ggplot2)
library(magrittr)
library(purrr)
library(readr)
library(shiny)
library(shinyjs)
library(shinyWidgets)
library(stringr)
library(tibble)
library(tidyr)
library(tools)

EXAMPLE <- TRUE
SHINYIO <- FALSE

POSITIVE <- "Yes"
NEGATIVE <- ""

if (SHINYIO) {
  # devtools::install_github("AhmedMehdiLab/glacier")
  # setRepositories()

  library(limma)
  library(Seurat)
}

# options
BARS_ANNO_MAX <- 80
OVER_ANNO_MAX <- 80
OVER_GENE_MAX <- 80
DEBOUNCE_TIME <- 1500

DT_LARGE <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 235px)")
DT_SMALL <- list(dom = "tr", paging = F, scrollCollapse = T, scrollY = "calc(100vh - 220.3px)")
REDUCTIONS <- c("Principal component analysis" = "pca",
              "Independent component analysis" = "ica",
              "t-distributed Stochastic Neighbor Embedding" = "tsne",
              "Uniform Manifold Approximation and Projection" = "umap")
options(shiny.maxRequestSize = 5 * 1024 ^ 3) # 5 GiB max upload size

warn_modal <- modalDialog(
  "This operation may take some time. Continue?",
  title = "Confirm action",
  footer = tagList(actionButton("warn.cont", "Continue"), modalButton("Cancel"))
)

sym_unify <- function(..., .true = T, .false = F) {
  data <- list(...)
  syms <- tibble::tibble(Symbol = unique(unlist(data, use.names = FALSE)))
  
  for (name in names(data)) syms[[name]] <- ifelse(syms$Symbol %in% data[[name]], .true, .false)
  return(syms)
}

server <- function(input, output, session) {
  showNotification(str_c("Welcome to glacier (", packageVersion("glacier"), ")!", collapse = ""), type = "message")
  store <- reactiveValues()
  
  # load data
  isolate({
    store$anno <- list() # Annotations
    store$cell <- list() # Seurat data
    store$data <- list() # Gene set contents
    store$expr <- list() # Expression data
    store$mart <- list() # biomaRt marts
    store$note <- list() # Notification handlers
    store$file <- reactiveValues() # Upload slot
    store$pause <- TRUE
    
    if (requireNamespace("E.PATH", quietly = TRUE)) {
      store$anno$`E.PATH` <- E.PATH::annotations
      store$data$`E.PATH` <- E.PATH::database
    }

    if (EXAMPLE) {
      store$anno$Example <- import_annotations(system.file("extdata", "ex_anno.csv", package = "glacier"), ",", TRUE, c(2, 4), 5)
      store$data$Example <- import_database(system.file("extdata", "ex_data.csv", package = "glacier"), ",", FALSE, c(2, 4), 0)
      store$cell$Example <- system.file("extdata", "ex_seurat.rds", package = "glacier") %>% readRDS
      store$expr$Example <- system.file("extdata", "ex_expr.rds", package = "glacier") %>% readRDS
    }
  })

  # process upload
  observeEvent(input$file.up, showModal(uploadUI("file")))
  uploadServer("file", reactive(input$file.up), store$file)
  observeEvent(input$confirm, {store[[store$file$type()]][[store$file$name()]] <- store$file$proc(); removeModal()})

  # controls [1]: source selection
  observe(updateSelectInput(session, "anno.source", NULL, names(store$anno)))
  observe(updateSelectInput(session, "cell.source", NULL, names(store$cell)))
  observe(updateSelectInput(session, "cell.group", NULL, names(cell_raw()@meta.data), c("group", "grp")))
  observe(updateSelectInput(session, "data.source", NULL, names(store$data)))
  observe(updateSelectInput(session, "expr.source", NULL, names(store$expr)))
  observe(toggleState("input.source", input$cell.source != ""))
  observe(toggleState("input.text", input$input.source == "text"))
  observe(toggleState("cell.group", input$cell.cluster != "_inter_"))
  observe(if (input$input.source == "cell") updateTextAreaInput(session, "input.text", value = str_c(input_proc()$gene, input_proc()$value, sep = " \t", collapse = "\n")))

  # process [1]: load from source
  anno_raw <- reactive(store$anno[[req(input$anno.source)]])
  cell_raw <- reactive({
    if (requireNamespace("Seurat", quietly = T)) {
      cell <- store$cell[[req(input$cell.source)]]
      if ("seurat_clusters" %in% names(cell[[]])) Seurat::Idents(cell) <- "seurat_clusters"
      return(cell)
  }})
  data_raw <- reactive(store$data[[req(input$data.source)]])
  expr_raw <- reactive(store$expr[[req(input$expr.source)]])
  info <- reactive((if (input$info.source == "anno") anno_raw() else data_raw()$gs_info) %>% select("name", "info"))
  
  # process [2]: extract information from source
  cell_clusts <- reactive({
    columns <- names(cell_raw()[[]])
    if ("anno" %in% columns) cell_raw()$anno %>% levels
    else if ("cell_type" %in% columns) cell_raw()$cell_type %>% levels
    else cell_raw()$seurat_clusters %>% levels
  })
  cell_groups <- reactive(cell_raw()@meta.data %>% pull(input$cell.group) %>% unique)
  universe <- reactive(data_raw()$gs_genes %>% unlist(use.names = F) %>% c(input_proc()$gene) %>% unique %>% length)

  # controls [2]: source configuration
  observe({
    opts <- tryCatch(
      if (is.null(cell_groups())) c("Inter-cluster" = "_inter_")
      else c("Inter-cluster" = "_inter_", "Intra-cluster" = "_intra_", setNames(cell_clusts(), str_c("Intra-cluster ", cell_clusts()))),
      error = function(e) c("Inter-cluster" = "_inter_")
    )
    
    updateSelectInput(session, "cell.cluster", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- setNames(cell_clusts(), str_c("Compare cluster ", cell_clusts()))
      else opts <- setNames(cell_groups(), str_c("Compare ", input$cell.group, " ", cell_groups())),
      error = function(e) NULL
    )

    updateSelectInput(session, "cell.select", NULL, opts)
  })
  observe({
    opts <- tryCatch(
      if (input$cell.cluster == "_inter_") opts <- c("Against all others" = "_all_", setNames(cell_clusts(), str_c("Against cluster ", cell_clusts())))
      else opts <- c("Against all others" = "_all_", setNames(cell_groups(), str_c("Against ", input$cell.group, " ", cell_groups()))),
      error = function(e) NULL
    )

    updateSelectInput(session, "cell.compare", NULL, opts)
  })
  observe({
    score_opts <- c()
    if (input$cell.source != "") score_opts <- c(score_opts, c("Seurat" = "cell"))
    if (input$expr.source != "") score_opts <- c(score_opts, c("Expression" = "expr"))
    
    updateSelectInput(session, "score.type", NULL, score_opts)
  })
  observe(updateSelectInput(session, "cell.view.cluster", NULL, cell_clusts(), cell_clusts()))
  observe({
    selected <- if (input$data.source == "MSigDB 7.4") c("C7 IMMUNESIGDB", "C7 VAX") else levels(data_raw()$gs_info$category)[1]
    updateSelectInput(session, "data.categories", NULL, levels(data_raw()$gs_info$category), selected)
  })
  observe({
    selected <- if (input$data.source == "MSigDB 7.4") "Homo sapiens" else levels(data_raw()$gs_info$organism)[1]
    updateSelectInput(session, "data.organisms", NULL, levels(data_raw()$gs_info$organism), selected)
  })
  observe(updateNumericInput(session, "input.universe", NULL, universe(), universe()))

  # process [3]: process source using inputs
  anno_proc <- reactive(glacier:::process_annotations(anno_raw(), info(), input$anno.types))
  cell_proc <- reactive(tryCatch(
    process_input_seurat(cell_raw(), input$cell.select, if (input$cell.compare != "_all_") input$cell.compare, if (input$cell.cluster != "_inter_") "grp", if (!input$cell.cluster %in% c("_inter_", "_intra_")) input$cell.cluster),
    error = function(e) tibble(gene = character(), value = numeric())
  ))
  data_proc <- reactive(glacier:::process_database(data_raw(), input$data.categories, input$data.organisms))
  text_proc <- reactive(input$input.text %>% process_input_text) %>% debounce(DEBOUNCE_TIME)
  cell_reductions <- reactive(REDUCTIONS[REDUCTIONS %in% names(cell_raw()@reductions)])
  input_proc <- reactive({
    input <- if (input$input.source == "text") text_proc() else cell_proc()
    
    if (!any(is.na(input$value)) || all(is.na(input$value))) removeNotification(store$note$values)
    else if (is.null(store$note$values)) store$note$values <- showNotification("Some genes do not have values, check 'Quality' tab", duration = NULL, type = "warning")
    return(input)
  })
  
  # process [4]: derive values from processed data
  anno_list <- reactive(tryCatch(anno_proc()$annos %>% str_subset(input$anno.regex) %>% as.character, error = function(e) character()))
  anno_sets <- reactive(anno_proc()$gs_anno$name)
  cell_gene <- reactive(rownames(cell_raw()))
  data_list <- reactive(data_proc()$genes)
  data_info <- reactive(data_proc()$gs_info %>% select(`Gene Set` = name, Information = info, Description = any_of("desc"), Category = category, Organism = organism))
  data_sets <- reactive(data_proc()$gs_info$name)
  expr_gene <- reactive(colnames(expr_raw()))

  cell_gene_list <- reactive(if (input$cell.gene.match) intersect(input_proc()$gene, cell_gene()) else cell_gene())
  cell_anno_list <- reactive(stats() %>% pull(Annotation) %>% str_sort(numeric = TRUE) %>% setNames(., nm = .) %>% map(~glacier:::explore_annotation(., anno_proc()$gs_annos, data_proc()$gs_genes, cell_gene_list())$genes) %>% compact())
  cell_anno_gene <- reactive(if (input$cell.anno == "") character() else cell_anno_list()[[input$cell.anno]] %>% str_sort(numeric = TRUE))

  score_anno_list <- reactive(stats() %>% pull(Annotation) %>% str_sort(numeric = TRUE) %>% setNames(., nm = .) %>% map(~glacier:::explore_annotation(., anno_proc()$gs_annos, data_proc()$gs_genes)) %>% compact())
  score_anno_gene <- reactive({if (input$score.anno == "") character() else score_anno_list()[[input$score.anno]]$genes})

  cont_sets <- reactive({
    sets <- sym_unify(Annotations = anno_sets(), Database = data_sets(), .true = POSITIVE, .false = NEGATIVE) %>%
      mutate(`Gene Set` = str_replace_all(Symbol, "_", ifelse(input$name.fix, " ", "_")), .keep = "unused") %>%
      select(`Gene Set`, Annotations, Database)
    
    if (all(sets[-1] == POSITIVE)) removeNotification(store$note$overlap)
    else if (is.null(store$note$overlap)) store$note$overlap <- showNotification("Gene set mismatch between annotations and database, check 'Quality' tab", duration = NULL, type = "warning")
    return(sets)
  })
  cont_gene <- reactive({
    input_proc() %>%
      mutate(Input = as.character(value)) %>%
      replace_na(list(Input = "No value")) %>%
      select(Symbol = gene, Input) %>%
      full_join(sym_unify(Database = data_list(), Seurat = cell_gene(), .true = POSITIVE, .false = NEGATIVE)) %>%
      rename(Gene = Symbol)
  })
  input_recognised <- reactive({
    number <- cont_gene() %>% filter(!is.na(Input) & (Database == POSITIVE)) %>% nrow
    
    if (number == nrow(input_proc())) removeNotification(store$note$recognised)
    else if (is.null(store$note$recognised)) store$note$recognised <- showNotification("Some genes were not recognised, check 'Quality' tab", duration = NULL, type = "warning")
    return(number)
  })

  # display [1]: input counts
  output$count <- renderText(str_c(
    "Input:      ", input_recognised(), " / ", nrow(input_proc()), " : ", sum(!is.na(input_proc()$value)), "\n",
    "Annotation: ", length(anno_sets()), " : ", length(anno_list()), "\n",
    "Database:   ", length(data_sets()), " : ", length(data_list()), "\n",
    "Seurat:     ", nrow(cell_raw()), " × ", ncol(cell_raw()), " : ", length(cell_clusts()), "\n",
    "Expression: ", nrow(expr_raw()), " × ", ncol(expr_raw()), "\n"
  ))

  output$cont.anno <- renderDataTable({if (store$pause) return(NULL); stats_raw()$stats %>% select("Annotation")}, DT_SMALL)
  output$cont.sets <- renderDataTable({if (store$pause) return(NULL); cont_sets()}, DT_SMALL)
  output$cont.gene <- renderDataTable({if (store$pause) return(NULL); cont_gene()}, DT_SMALL)

  # controls [3]: display configuration
  observe(updateSelectInput(session, "cell.anno", NULL, names(cell_anno_list())))
  observe(updateSelectInput(session, "cell.genes", NULL, cell_anno_gene(), head(cell_anno_gene(), 16)))
  observe(updateSelectInput(session, "cell.overview", NULL, cell_reductions()))
  observe(updateSelectInput(session, "score.anno", NULL, names(cell_anno_list())))

  # process [5]: calculate statistics
  stats_raw <- reactive({
    withProgress(message = "Calculating statistics", {
      setProgress(detail = "Finding matches")
      pre <- glacier:::calculate_pre(input_proc(), anno_list(), anno_proc()$gs_annos, data_proc()$gs_genes)

      setProgress(value = 0.8, detail = "Computing significance")
      post <- glacier:::calculate_post(pre$stats_pre, nrow(input_proc()), max(input$input.universe, universe()))

      list(stats = post, matches = pre$matches)
    })
  })
  stats <- reactive({
    if (input$stat.sigonly) stats_raw()$stats %>% filter(`Adjusted P-value` <= 0.05)
    else stats_raw()$stats
  })
  matches <- reactive({stats_raw()$matches})

  scores_stat <- reactive({
    if (!requireNamespace("mixOmics", quietly = T) || !requireNamespace("pROC", quietly = T)) return(NULL)

    validate(
      need("grp" %in% names(cell_raw()[[]]), "No groups found in Seurat object")
    )
    
    withProgress(message = "Calculating scores", {
      if (input$score.type == "cell") score_seurat(cell_raw(), input$cell.group, score_anno_gene())
      else if (input$score.type == "expr") score_expr(expr_raw(), "grp", score_anno_gene())
      else stop("scores_stat: invalid option chosen")
    })
  })

  bars_stat <- reactive(stats() %>% arrange(
    if (input$bars.anno.order == "Annotation") str_to_lower(Annotation)
    else if (input$bars.anno.order %in% c("P-value", "Adjusted P-value")) .data[[input$bars.anno.order]]
    else desc(.data[[input$bars.anno.order]])
  ))
  over_stat <- reactive(stats() %>% arrange(
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
  observe(updateVarSelectInput(session, "stat.columns", NULL, stats(), names(stats())))
  observe(toggleState("cell.gene.cluster", input$cell.gene.match))

  # warn user about time-consuming operations
  observe({
    if (input$main_view %in% c("score_pane", "conv_pane", "cont_pane")) {
      store$pause <- TRUE
      showModal(warn_modal)
    }
  })
  observeEvent(input$warn.cont, {if (input$warn.cont) store$pause <- FALSE; removeModal()})
  
  # crop data
  view_bars <- reactive(bars_stat()[store$bars.anno.select[1]:store$bars.anno.select[2], ])
  view_over <- reactive(over_stat()[store$over.anno.select[1]:store$over.anno.select[2], ])
  view_over_gene <- reactive(over_input()[store$over.gene.select[1]:store$over.gene.select[2], ])
  view_cell_gene <- reactive(input$cell.genes) %>% debounce(DEBOUNCE_TIME)
  view_cell_gene_clust <- reactive(cell_raw() %>% Seurat::FindAllMarkers(features = view_cell_gene()) %>% group_by(cluster) %>% pull(gene))
  view_cell_clust <- reactive(input$cell.view.cluster) %>% debounce(DEBOUNCE_TIME)

  # secondary information
  view_table <- function(table, cols, split) table %>% select(!!!cols) %>% mutate(across(!!split, str_replace_all, "_", ifelse(input$name.fix, " ", "_")))
  output$bars <- renderPlot({
    plot_stats(view_bars(), input$bars.value, input$bars.color, input$bars.anno.sort) +
      geom_col(colour = "black") + ylab(NULL) + scale_fill_gradient(low = "black", high = "white") +
      theme_classic()
  })
  output$over <- renderPlot({
    plot_overlap(matches(), input$over.color, view_over_gene(), view_over(), "top") +
      xlab(NULL) + ylab(NULL) +
      scale_fill_gradientn(na.value = "transparent", colours = grDevices::hcl.colors(3, palette = "Blue-Red 2")) +
      theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5))
  })
  output$cell <- renderPlot(if (requireNamespace("Seurat", quietly = T)) Seurat::DimPlot(cell_raw(), reduction = input$cell.overview, label = T, repel = T))
  output$heat <- renderPlot({
    if (!requireNamespace("Seurat", quietly = T) || !length(view_cell_gene())) return(NULL)

    feats <- if (input$cell.gene.match && input$cell.gene.cluster) view_cell_gene_clust() else view_cell_gene()
    if (input$cell.plot != "heat") feats <- unique(feats)

    width <- feats %>% length %>% sqrt %>% ceiling
    sample <- subset(cell_raw(), idents = view_cell_clust())
    if (input$cell.plot == "dot") Seurat::DotPlot(sample, features = feats)
    else if (input$cell.plot == "feat") Seurat::FeaturePlot(sample, features = feats, ncol = width, label = TRUE, repel = TRUE)
    else if (input$cell.plot == "heat") Seurat::DoHeatmap(sample, features = feats)
    else if (input$cell.plot == "ridge") Seurat::RidgePlot(sample, features = feats, ncol = width)
    else if (input$cell.plot == "violin") Seurat::VlnPlot(sample, features = feats, ncol = width)
  })
  output$score <- renderPlot({
    if (store$pause) return(NULL)
    
    x <- if (input$score.type == "expr") "bin" else "seurat_clusters"
    plot_scores(scores_stat()$scores, x, input$score.method, "grp", input$score.plot) +
      xlab(NULL) + ylab(str_c(input$score.anno, " scores")) + theme(legend.position = "top")
  })
  output$rocs <- renderPlot({if (store$pause) return(NULL); plot_auc(scores_stat()$aucs, input$score.method)})
  output$conv <- renderDataTable({
    if (store$pause) return(NULL)
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
      if (input$conv == "mh") genes <- biomaRt::getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = vals, mart = store$mart$mm, attributesL = "hgnc_symbol", martL = store$mart$hs)
      if (input$conv == "hm") genes <- biomaRt::getLDS(attributes = "hgnc_symbol", filters = "hgnc_symbol", values = vals, mart = store$mart$hs, attributesL = "mgi_symbol", martL = store$mart$mm)
    })

    return(genes)
  }, DT_LARGE)
  output$stat <- renderDataTable(view_table(stats(), input$stat.columns, "Annotation"), DT_LARGE)
  output$info <- renderDataTable(view_table(data_info(), input$info.columns, "Gene Set"), DT_LARGE)

  output$file.down <- downloadHandler("glacier_results.csv", . %>% write_csv(stats(), .))
}
