# ============================================================================
# Complete WGCNA Analysis Pipeline
# ============================================================================
# Full workflow from dissimilarity matrix to enrichment analysis
# ============================================================================

library(shiny)
library(data.table)
library(DT)
library(ggplot2)
library(plotly)
library(gplots)
library(dplyr)

# Check Bioconductor packages
check_bioc_packages <- function() {
  required <- c("clusterProfiler", "org.Hs.eg.db")
  missing <- c()
  for (pkg in required) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  if (length(missing) > 0) {
    return(list(available = FALSE, missing = missing))
  }
  return(list(available = TRUE, missing = c()))
}

# ==============================================================================
# UI
# ==============================================================================
ui <- fluidPage(
  titlePanel("Complete WGCNA Analysis Pipeline"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      h4("Step 1: Load & Cluster"),
      actionButton("loadBtn", "Load Dissimilarity Matrix", 
                   class = "btn-primary btn-lg",
                   style = "width: 100%; margin-bottom: 10px;"),
      actionButton("clusterBtn", "Run Clustering",
                   class = "btn-success btn-lg",
                   style = "width: 100%; margin-bottom: 20px;"),
      
      conditionalPanel(
        condition = "output.clusteringDone",
        hr(),
        
        h4("Step 2: Detect Modules"),
        sliderInput("cutHeight", "Cut Height:", 
                    min = 0, max = 1, value = 0.25, step = 0.01),
        numericInput("minModuleSize", "Min Module Size:",
                     value = 30, min = 10, max = 200, step = 10),
        actionButton("cutBtn", "Detect Modules",
                     class = "btn-info",
                     style = "width: 100%; margin-bottom: 20px;"),
        
        hr(),
        
        h4("Step 3: Deep Analysis"),
        actionButton("runAnalysis", "Run Module Analysis",
                     class = "btn-warning btn-lg",
                     style = "width: 100%; margin-bottom: 10px;"),
        actionButton("runEnrichment", "Run GO/KEGG Enrichment",
                     class = "btn-danger btn-lg",
                     style = "width: 100%; margin-bottom: 20px;"),
        
        hr(),
        
        h4("Export"),
        downloadButton("downloadModules", "Download Modules",
                       style = "width: 100%; margin-bottom: 10px;"),
        downloadButton("downloadHub", "Download Hub Genes",
                       style = "width: 100%; margin-bottom: 10px;"),
        downloadButton("downloadStats", "Download Statistics",
                       style = "width: 100%;")
      )
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "mainTabs",
        
        # Tab 1: Overview & Log
        tabPanel("Overview",
                 icon = icon("info-circle"),
                 br(),
                 h3("Pipeline Status"),
                 verbatimTextOutput("statusLog"),
                 
                 hr(),
                 
                 conditionalPanel(
                   condition = "output.modulesDetected",
                   fluidRow(
                     column(6,
                            h3("Data Summary"),
                            tableOutput("dataSummary")
                     ),
                     column(6,
                            h3("Module Size Distribution"),
                            plotlyOutput("moduleSizePlot", height = "300px")
                     )
                   )
                 )
        ),
        
        # Tab 2: Dendrogram
        tabPanel("Dendrogram",
                 icon = icon("tree"),
                 br(),
                 
                 conditionalPanel(
                   condition = "!output.clusteringDone",
                   div(style = "padding: 50px; text-align: center;",
                       h3("Clustering Not Performed"),
                       p("Click 'Run Clustering' in the sidebar."))
                 ),
                 
                 conditionalPanel(
                   condition = "output.clusteringDone",
                   h3("Hierarchical Clustering Dendrogram"),
                   plotOutput("dendroPlot", height = "700px")
                 )
        ),
        
        # Tab 3: Module Table
        tabPanel("Module Genes",
                 icon = icon("list"),
                 br(),
                 
                 conditionalPanel(
                   condition = "!output.modulesDetected",
                   div(style = "padding: 50px; text-align: center;",
                       h3("Modules Not Detected"),
                       p("Complete clustering and module detection first."))
                 ),
                 
                 conditionalPanel(
                   condition = "output.modulesDetected",
                   h3("Gene-Module Assignments"),
                   DTOutput("moduleTable"),
                   
                   hr(),
                   
                   h4("Module Statistics"),
                   tableOutput("moduleStats")
                 )
        ),
        
        # Tab 4: Hub Genes
        tabPanel("Hub Genes",
                 icon = icon("star"),
                 br(),
                 
                 conditionalPanel(
                   condition = "!output.analysisComplete",
                   div(style = "padding: 50px; text-align: center;",
                       h3("Analysis Not Complete"),
                       p("Click 'Run Module Analysis' in the sidebar."))
                 ),
                 
                 conditionalPanel(
                   condition = "output.analysisComplete",
                   h3("Hub Genes by Module"),
                   
                   fluidRow(
                     column(3,
                            selectInput("hubModule", "Select Module:",
                                        choices = NULL)
                     ),
                     column(3,
                            numericInput("topN", "Top N Hubs:",
                                         value = 10, min = 5, max = 100)
                     )
                   ),
                   
                   hr(),
                   
                   DTOutput("hubTable")
                 )
        ),
        
        # Tab 5: Module Quality
        tabPanel("Module Quality",
                 icon = icon("chart-bar"),
                 br(),
                 
                 conditionalPanel(
                   condition = "!output.analysisComplete",
                   div(style = "padding: 50px; text-align: center;",
                       h3("Analysis Not Complete"),
                       p("Click 'Run Module Analysis' in the sidebar."))
                 ),
                 
                 conditionalPanel(
                   condition = "output.analysisComplete",
                   h3("Module Quality Metrics"),
                   DTOutput("statsTable"),
                   
                   hr(),
                   
                   fluidRow(
                     column(6,
                            h4("Separation Score"),
                            plotlyOutput("separationPlot", height = "400px")
                     ),
                     column(6,
                            h4("Intra vs Inter Similarity"),
                            plotlyOutput("similarityPlot", height = "400px")
                     )
                   )
                 )
        ),
        
        # Tab 6: Enrichment
        tabPanel("Enrichment",
                 icon = icon("dna"),
                 br(),
                 
                 conditionalPanel(
                   condition = "!output.enrichmentComplete",
                   div(style = "padding: 50px; text-align: center;",
                       h3("Enrichment Not Run"),
                       p("Click 'Run GO/KEGG Enrichment' in the sidebar."))
                 ),
                 
                 conditionalPanel(
                   condition = "output.enrichmentComplete",
                   h3("GO and KEGG Enrichment Results"),
                   
                   fluidRow(
                     column(3,
                            selectInput("enrichModule", "Select Module:",
                                        choices = NULL)
                     ),
                     column(3,
                            radioButtons("enrichType", "Type:",
                                         choices = c("GO" = "go", "KEGG" = "kegg"),
                                         selected = "go")
                     )
                   ),
                   
                   hr(),
                   
                   h4("Top Enriched Terms"),
                   DTOutput("enrichTable")
                 )
        )
      )
    )
  )
)

# ==============================================================================
# Server
# ==============================================================================
server <- function(input, output, session) {
  
  # Reactive values
  vals <- reactiveValues(
    dissim_matrix = NULL,
    sim_matrix = NULL,
    gene_names = NULL,
    hclustObj = NULL,
    modules = NULL,
    
    hub_genes = NULL,
    module_stats = NULL,
    kME_all = NULL,
    
    go_results = NULL,
    kegg_results = NULL,
    
    log = "Welcome! Click 'Load Dissimilarity Matrix' to begin.\n",
    dataLoaded = FALSE,
    clusteringDone = FALSE,
    modulesDetected = FALSE,
    analysisComplete = FALSE,
    enrichmentComplete = FALSE
  )
  
  # Output flags
  output$dataLoaded <- reactive({ vals$dataLoaded })
  output$clusteringDone <- reactive({ vals$clusteringDone })
  output$modulesDetected <- reactive({ vals$modulesDetected })
  output$analysisComplete <- reactive({ vals$analysisComplete })
  output$enrichmentComplete <- reactive({ vals$enrichmentComplete })
  
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  outputOptions(output, "clusteringDone", suspendWhenHidden = FALSE)
  outputOptions(output, "modulesDetected", suspendWhenHidden = FALSE)
  outputOptions(output, "analysisComplete", suspendWhenHidden = FALSE)
  outputOptions(output, "enrichmentComplete", suspendWhenHidden = FALSE)
  
  # --------------------------------------------------------------------------
  # Step 1: Load Data
  # --------------------------------------------------------------------------
  observeEvent(input$loadBtn, {
    vals$log <- paste0(vals$log, "\n=== Loading Data ===\n")
    
    withProgress(message = 'Loading dissimilarity matrix...', value = 0, {
      tryCatch({
        if (!file.exists("dissimilarity_matrix.csv")) {
          stop("dissimilarity_matrix.csv not found!")
        }
        
        incProgress(0.3, detail = "Reading CSV...")
        vals$log <- paste0(vals$log, "Loading dissimilarity_matrix.csv...\n")
        
        data <- fread("dissimilarity_matrix.csv", data.table = FALSE)
        
        incProgress(0.6, detail = "Processing...")
        vals$gene_names <- data[[1]]
        data <- data[, -1]
        rownames(data) <- vals$gene_names
        colnames(data) <- vals$gene_names
        vals$dissim_matrix <- as.matrix(data)
        vals$sim_matrix <- 1 - vals$dissim_matrix
        
        incProgress(1.0, detail = "Done!")
        
        vals$log <- paste0(vals$log, sprintf("Loaded %d x %d matrix\n", 
                                             nrow(vals$dissim_matrix), 
                                             ncol(vals$dissim_matrix)))
        vals$log <- paste0(vals$log, "Next: Click 'Run Clustering'\n")
        
        vals$dataLoaded <- TRUE
        
      }, error = function(e) {
        vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      })
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 2: Clustering
  # --------------------------------------------------------------------------
  observeEvent(input$clusterBtn, {
    req(vals$dataLoaded)
    
    vals$log <- paste0(vals$log, "\n=== Running Clustering ===\n")
    vals$log <- paste0(vals$log, "This may take 5-15 minutes...\n")
    
    withProgress(message = 'Clustering...', value = 0, {
      tryCatch({
        incProgress(0.4, detail = "Running hclust...")
        
        start_time <- Sys.time()
        dist_obj <- as.dist(vals$dissim_matrix)
        vals$hclustObj <- hclust(dist_obj, method = "average")
        end_time <- Sys.time()
        
        time_taken <- as.numeric(difftime(end_time, start_time, units = "mins"))
        tree_height <- max(vals$hclustObj$height)
        suggested_cut <- tree_height * 0.75
        
        incProgress(1.0, detail = "Done!")
        
        vals$log <- paste0(vals$log, sprintf("Completed in %.1f minutes\n", time_taken))
        vals$log <- paste0(vals$log, sprintf("Tree height: %.4f\n", tree_height))
        vals$log <- paste0(vals$log, sprintf("Suggested cut: %.4f\n", suggested_cut))
        vals$log <- paste0(vals$log, "Next: Adjust cut height and detect modules\n")
        
        updateSliderInput(session, "cutHeight",
                          min = 0, max = tree_height,
                          value = suggested_cut,
                          step = tree_height / 100)
        
        vals$clusteringDone <- TRUE
        
      }, error = function(e) {
        vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      })
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 3: Detect Modules
  # --------------------------------------------------------------------------
  observeEvent(input$cutBtn, {
    req(vals$clusteringDone)
    
    vals$log <- paste0(vals$log, "\n=== Detecting Modules ===\n")
    vals$log <- paste0(vals$log, sprintf("Cut height: %.4f\n", input$cutHeight))
    vals$log <- paste0(vals$log, sprintf("Min size: %d\n", input$minModuleSize))
    
    withProgress(message = 'Detecting modules...', value = 0, {
      tryCatch({
        incProgress(0.5, detail = "Cutting tree...")
        
        module_labels <- cutree(vals$hclustObj, h = input$cutHeight)
        initial_modules <- length(unique(module_labels))
        
        vals$log <- paste0(vals$log, sprintf("Initial: %d clusters\n", initial_modules))
        
        incProgress(0.8, detail = "Filtering...")
        
        module_sizes <- table(module_labels)
        small_modules <- as.numeric(names(module_sizes)[module_sizes < input$minModuleSize])
        
        for (mod in small_modules) {
          module_labels[module_labels == mod] <- 0
        }
        
        unique_modules <- sort(unique(module_labels))
        unique_modules <- unique_modules[unique_modules != 0]
        
        new_labels <- module_labels
        for (i in seq_along(unique_modules)) {
          new_labels[module_labels == unique_modules[i]] <- i
        }
        
        vals$modules <- data.frame(
          GeneID = vals$gene_names,
          ModuleID = new_labels,
          stringsAsFactors = FALSE
        )
        
        incProgress(1.0, detail = "Done!")
        
        n_modules <- length(unique(new_labels[new_labels != 0]))
        n_unclustered <- sum(new_labels == 0)
        
        vals$log <- paste0(vals$log, sprintf("Detected %d modules\n", n_modules))
        vals$log <- paste0(vals$log, sprintf("Unclustered: %d genes\n", n_unclustered))
        
        if (n_modules == 0) {
          vals$log <- paste0(vals$log, "WARNING: No modules! Try different parameters.\n")
        } else {
          vals$log <- paste0(vals$log, "Next: Run Module Analysis\n")
        }
        
        vals$modulesDetected <- TRUE
        
      }, error = function(e) {
        vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      })
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 4: Module Analysis
  # --------------------------------------------------------------------------
  observeEvent(input$runAnalysis, {
    req(vals$modulesDetected)
    
    vals$log <- paste0(vals$log, "\n=== Module Analysis ===\n")
    
    withProgress(message = 'Analyzing modules...', value = 0, {
      tryCatch({
        modules <- vals$modules[vals$modules$ModuleID != 0, ]
        module_ids <- sort(unique(modules$ModuleID))
        
        vals$log <- paste0(vals$log, sprintf("Analyzing %d modules\n", length(module_ids)))
        
        # Calculate kME
        incProgress(0.2, detail = "Calculating connectivity...")
        kME_results <- list()
        
        for (mod_id in module_ids) {
          mod_genes <- modules$GeneID[modules$ModuleID == mod_id]
          mod_indices <- which(vals$gene_names %in% mod_genes)
          
          kME <- numeric(length(vals$gene_names))
          for (i in 1:length(vals$gene_names)) {
            kME[i] <- mean(vals$sim_matrix[i, mod_indices])
          }
          
          kME_results[[as.character(mod_id)]] <- data.frame(
            GeneID = vals$gene_names,
            ModuleID = mod_id,
            kME = kME,
            stringsAsFactors = FALSE
          )
        }
        
        vals$kME_all <- do.call(rbind, kME_results)
        
        # Hub genes
        incProgress(0.5, detail = "Identifying hub genes...")
        hub_genes <- list()
        
        for (mod_id in module_ids) {
          mod_data <- vals$kME_all[vals$kME_all$ModuleID == mod_id, ]
          mod_genes_only <- modules$GeneID[modules$ModuleID == mod_id]
          mod_data <- mod_data[mod_data$GeneID %in% mod_genes_only, ]
          mod_data <- mod_data[order(mod_data$kME, decreasing = TRUE), ]
          
          n_hubs <- max(5, ceiling(nrow(mod_data) * 0.1))
          hubs <- head(mod_data, n_hubs)
          hubs$HubRank <- 1:nrow(hubs)
          
          hub_genes[[as.character(mod_id)]] <- hubs
        }
        
        vals$hub_genes <- do.call(rbind, hub_genes)
        
        # Module statistics
        incProgress(0.8, detail = "Calculating quality...")
        module_stats <- data.frame()
        
        for (mod_id in module_ids) {
          mod_genes <- modules$GeneID[modules$ModuleID == mod_id]
          mod_indices <- which(vals$gene_names %in% mod_genes)
          n_genes <- length(mod_indices)
          
          if (n_genes > 1) {
            mod_sim_matrix <- vals$sim_matrix[mod_indices, mod_indices]
            diag(mod_sim_matrix) <- NA
            
            intra_sim <- mean(mod_sim_matrix, na.rm = TRUE)
            intra_sim_sd <- sd(as.vector(mod_sim_matrix), na.rm = TRUE)
            
            other_indices <- setdiff(1:length(vals$gene_names), mod_indices)
            if (length(other_indices) > 0) {
              inter_sim <- mean(vals$sim_matrix[mod_indices, other_indices])
            } else {
              inter_sim <- NA
            }
            
            separation <- intra_sim - inter_sim
          } else {
            intra_sim <- NA
            intra_sim_sd <- NA
            inter_sim <- NA
            separation <- NA
          }
          
          mod_hubs <- vals$hub_genes[vals$hub_genes$ModuleID == mod_id, ]
          top_hub <- if(nrow(mod_hubs) > 0) mod_hubs$GeneID[1] else NA
          
          module_stats <- rbind(module_stats, data.frame(
            ModuleID = mod_id,
            GeneCount = n_genes,
            IntraModuleSimilarity = intra_sim,
            IntraModuleSimilarity_SD = intra_sim_sd,
            InterModuleSimilarity = inter_sim,
            Separation = separation,
            TopHubGene = top_hub,
            stringsAsFactors = FALSE
          ))
        }
        
        vals$module_stats <- module_stats[order(module_stats$Separation, decreasing = TRUE), ]
        
        incProgress(1.0, detail = "Done!")
        
        vals$log <- paste0(vals$log, sprintf("Found %d hub genes\n", nrow(vals$hub_genes)))
        vals$log <- paste0(vals$log, "Analysis complete!\n")
        vals$log <- paste0(vals$log, "Next: (Optional) Run GO/KEGG Enrichment\n")
        
        vals$analysisComplete <- TRUE
        
        # Update selectors
        updateSelectInput(session, "hubModule", 
                          choices = module_ids, 
                          selected = module_ids[1])
        
      }, error = function(e) {
        vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      })
    })
  })
  
  # --------------------------------------------------------------------------
  # Step 5: Enrichment
  # --------------------------------------------------------------------------
  observeEvent(input$runEnrichment, {
    req(vals$analysisComplete)
    
    vals$log <- paste0(vals$log, "\n=== Enrichment Analysis ===\n")
    
    pkg_check <- check_bioc_packages()
    if (!pkg_check$available) {
      vals$log <- paste0(vals$log, "ERROR: Missing packages: ", 
                         paste(pkg_check$missing, collapse = ", "), "\n")
      vals$log <- paste0(vals$log, "\nInstall with:\n")
      vals$log <- paste0(vals$log, "  install.packages('BiocManager')\n")
      vals$log <- paste0(vals$log, "  BiocManager::install(c('", 
                         paste(pkg_check$missing, collapse = "','"), "'))\n")
      
      showModal(modalDialog(
        title = "Missing Packages",
        paste("Please install:", paste(pkg_check$missing, collapse = ", ")),
        easyClose = TRUE
      ))
      return()
    }
    
    vals$log <- paste0(vals$log, "Loading libraries...\n")
    
    tryCatch({
      library(clusterProfiler)
      library(org.Hs.eg.db)
    }, error = function(e) {
      vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      return()
    })
    
    withProgress(message = 'Running enrichment...', value = 0, {
      tryCatch({
        modules <- vals$modules[vals$modules$ModuleID != 0, ]
        modules$ENSEMBL <- gsub("\\..*", "", modules$GeneID)
        
        incProgress(0.2, detail = "Converting IDs...")
        gene_mapping <- bitr(modules$ENSEMBL, 
                             fromType = "ENSEMBL",
                             toType = c("ENTREZID", "SYMBOL"),
                             OrgDb = org.Hs.eg.db)
        
        modules <- merge(modules, gene_mapping, by = "ENSEMBL", all.x = TRUE)
        modules <- modules[!is.na(modules$ENTREZID), ]
        
        vals$log <- paste0(vals$log, sprintf("Mapped %d genes\n", nrow(modules)))
        
        module_ids <- sort(unique(modules$ModuleID))
        go_results <- list()
        kegg_results <- list()
        
        for (i in seq_along(module_ids)) {
          mod_id <- module_ids[i]
          incProgress(0.3 + (i / length(module_ids)) * 0.5, 
                      detail = sprintf("Module %d...", mod_id))
          
          mod_genes <- modules$ENTREZID[modules$ModuleID == mod_id]
          if (length(mod_genes) < 5) next
          
          # GO
          go_bp <- tryCatch({
            enrichGO(gene = mod_genes, OrgDb = org.Hs.eg.db,
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                     readable = TRUE)
          }, error = function(e) NULL)
          
          if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
            go_results[[paste0("Module_", mod_id)]] <- as.data.frame(go_bp)
          }
          
          # KEGG
          kegg <- tryCatch({
            enrichKEGG(gene = mod_genes, organism = 'hsa',
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2)
          }, error = function(e) NULL)
          
          if (!is.null(kegg) && nrow(as.data.frame(kegg)) > 0) {
            kegg_results[[paste0("Module_", mod_id)]] <- as.data.frame(kegg)
          }
        }
        
        vals$go_results <- go_results
        vals$kegg_results <- kegg_results
        
        incProgress(1.0, detail = "Done!")
        
        vals$log <- paste0(vals$log, sprintf("GO: %d modules\n", length(go_results)))
        vals$log <- paste0(vals$log, sprintf("KEGG: %d modules\n", length(kegg_results)))
        vals$log <- paste0(vals$log, "Enrichment complete!\n")
        
        vals$enrichmentComplete <- TRUE
        
        # Update selector
        enrich_modules <- unique(c(
          gsub("Module_", "", names(go_results)),
          gsub("Module_", "", names(kegg_results))
        ))
        updateSelectInput(session, "enrichModule", 
                          choices = enrich_modules, 
                          selected = enrich_modules[1])
        
      }, error = function(e) {
        vals$log <- paste0(vals$log, "ERROR: ", e$message, "\n")
      })
    })
  })
  
  # --------------------------------------------------------------------------
  # Outputs
  # --------------------------------------------------------------------------
  output$statusLog <- renderText({ vals$log })
  
  output$dataSummary <- renderTable({
    req(vals$modulesDetected)
    
    module_counts <- table(vals$modules$ModuleID)
    module_counts_no_zero <- module_counts[names(module_counts) != "0"]
    
    data.frame(
      Metric = c("Total Genes", "Modules", "Largest", "Smallest", "Unclustered"),
      Value = c(
        length(vals$gene_names),
        length(module_counts_no_zero),
        if(length(module_counts_no_zero) > 0) max(module_counts_no_zero) else 0,
        if(length(module_counts_no_zero) > 0) min(module_counts_no_zero) else 0,
        sum(vals$modules$ModuleID == 0)
      )
    )
  }, striped = TRUE, hover = TRUE)
  
  output$moduleSizePlot <- renderPlotly({
    req(vals$modulesDetected)
    
    module_counts <- as.data.frame(table(vals$modules$ModuleID))
    module_counts <- module_counts[module_counts$Var1 != 0, ]
    
    if (nrow(module_counts) == 0) return(plotly_empty())
    
    colnames(module_counts) <- c("Module", "GeneCount")
    module_counts <- module_counts[order(module_counts$GeneCount, decreasing = TRUE), ]
    module_counts$Module <- factor(module_counts$Module, levels = module_counts$Module)
    
    p <- ggplot(module_counts, aes(x = Module, y = GeneCount, fill = GeneCount)) +
      geom_bar(stat = "identity") +
      scale_fill_gradient(low = "lightblue", high = "darkblue") +
      theme_minimal() +
      labs(x = "Module ID", y = "Genes") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p)
  })
  
  output$dendroPlot <- renderPlot({
    req(vals$clusteringDone)
    
    max_display <- min(10000, length(vals$hclustObj$order))
    
    if (length(vals$hclustObj$order) > max_display) {
      subset_idx <- 1:max_display
      subset_dist <- as.dist(vals$dissim_matrix[subset_idx, subset_idx])
      hc_display <- hclust(subset_dist, method = "average")
      title_suffix <- sprintf(" (first %d genes)", max_display)
    } else {
      hc_display <- vals$hclustObj
      title_suffix <- ""
    }
    
    plot(hc_display,
         main = paste0("Dendrogram", title_suffix),
         xlab = "Genes", ylab = "Height",
         sub = "", labels = FALSE,
         ylim = c(0, max(hc_display$height) * 1.05))
    
    abline(h = input$cutHeight, col = "red", lty = 2, lwd = 2)
    text(x = par("usr")[1] + (par("usr")[2] - par("usr")[1]) * 0.02,
         y = input$cutHeight,
         labels = sprintf("Cut = %.2f", input$cutHeight),
         pos = 3, col = "red", cex = 1.2, font = 2)
  })
  
  output$moduleTable <- renderDT({
    req(vals$modulesDetected)
    datatable(vals$modules,
              options = list(pageLength = 25, scrollX = TRUE),
              filter = 'top', rownames = FALSE)
  })
  
  output$moduleStats <- renderTable({
    req(vals$modulesDetected)
    module_summary <- as.data.frame(table(vals$modules$ModuleID))
    colnames(module_summary) <- c("ModuleID", "GeneCount")
    module_summary[order(module_summary$GeneCount, decreasing = TRUE), ]
  }, striped = TRUE, hover = TRUE)
  
  output$hubTable <- renderDT({
    req(vals$hub_genes, input$hubModule)
    hub_data <- vals$hub_genes[vals$hub_genes$ModuleID == input$hubModule, ]
    hub_data <- head(hub_data, input$topN)
    datatable(hub_data[, c("GeneID", "kME", "HubRank")],
              options = list(pageLength = 25), rownames = FALSE)
  })
  
  output$statsTable <- renderDT({
    req(vals$module_stats)
    datatable(vals$module_stats,
              options = list(pageLength = 25, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = c("IntraModuleSimilarity", "InterModuleSimilarity", 
                              "Separation"), digits = 4)
  })
  
  output$separationPlot <- renderPlotly({
    req(vals$module_stats)
    p <- ggplot(vals$module_stats, 
                aes(x = reorder(factor(ModuleID), -Separation), y = Separation)) +
      geom_point(size = 3, color = "steelblue") +
      geom_hline(yintercept = 0.15, linetype = "dashed", color = "red") +
      theme_minimal() +
      labs(x = "Module", y = "Separation")
    ggplotly(p)
  })
  
  output$similarityPlot <- renderPlotly({
    req(vals$module_stats)
    p <- ggplot(vals$module_stats, 
                aes(x = IntraModuleSimilarity, y = InterModuleSimilarity, 
                    color = Separation, size = GeneCount)) +
      geom_point(alpha = 0.7) +
      scale_color_gradient(low = "red", high = "blue") +
      theme_minimal() +
      labs(x = "Intra-Module Similarity", y = "Inter-Module Similarity")
    ggplotly(p)
  })
  
  output$enrichTable <- renderDT({
    req(vals$enrichmentComplete, input$enrichModule, input$enrichType)
    
    if (input$enrichType == "go") {
      result <- vals$go_results[[paste0("Module_", input$enrichModule)]]
    } else {
      result <- vals$kegg_results[[paste0("Module_", input$enrichModule)]]
    }
    
    if (is.null(result)) {
      return(data.frame(Message = "No enrichment found"))
    }
    
    datatable(result[, c("Description", "GeneRatio", "pvalue", "qvalue")],
              options = list(pageLength = 25), rownames = FALSE) %>%
      formatSignif(columns = c("pvalue", "qvalue"), digits = 3)
  })
  
  # --------------------------------------------------------------------------
  # Downloads
  # --------------------------------------------------------------------------
  output$downloadModules <- downloadHandler(
    filename = function() { paste0("modules_", Sys.Date(), ".csv") },
    content = function(file) {
      req(vals$modules)
      write.csv(vals$modules, file, row.names = FALSE)
    }
  )
  
  output$downloadHub <- downloadHandler(
    filename = function() { paste0("hub_genes_", Sys.Date(), ".csv") },
    content = function(file) {
      req(vals$hub_genes)
      write.csv(vals$hub_genes, file, row.names = FALSE)
    }
  )
  
  output$downloadStats <- downloadHandler(
    filename = function() { paste0("module_stats_", Sys.Date(), ".csv") },
    content = function(file) {
      req(vals$module_stats)
      write.csv(vals$module_stats, file, row.names = FALSE)
    }
  )
}

# Run app
shinyApp(ui = ui, server = server)