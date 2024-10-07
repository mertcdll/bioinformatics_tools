suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(EnhancedVolcano))
suppressPackageStartupMessages(library(ggplot2))




DiffexpWedgeR <- function(sample_info_file, output_dir, data_source, interaction, cutoff_lfc, cutoff_p) {
  
  cutoff_p <- as.double(cutoff_p)
  cutoff_lfc <- as.double(cutoff_lfc)
  interaction <- interaction
  sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  sample_names <- sample_info$names # sample names 
  count_dirs <- sample_info$count_dirs # count directories
  factors <- sample_info[, !names(sample_info) %in% c("names", "subject", "count_dirs")]
  rownames(factors) <- sample_names
  subjects <- sample_info$subject #subjects
  factor_names <- names(factors) #factor names
  factor_count <- length(factor_names) #number of factors
  
  exp_dict <- {}
  
  # Create raw counts dictionary # do not change
  for (i in 1:length(sample_names)) {
    count_table <- read.table(count_dirs[i], header = TRUE, row.names = 1, sep = "\t")
    exp_dict[[sample_names[i]]] <- count_table
  }
  
  # Convert it to table # do not change
  exp_table <- data.frame(exp_dict)
  
  
  # Convert floats to integers # do not change
  for (col in names(exp_table)) {   
    if (is.numeric(exp_table[[col]])) {
      exp_table[[col]] <- as.integer(exp_table[[col]])
    }
  } 
  
  
  raw_counts_table <- file.path(output_dir, "raw_counts_table.txt") # do not change
  
  # write raw unfiltered counts # do not change
  write.table(as_tibble(exp_table, rownames = "GeneID"), file = raw_counts_table, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Get groups # change this one 
  factors <- data.frame(lapply(factors, as.factor)) #converts all columns to factors
  grp <- as.factor(do.call(paste, c(factors, sep="_"))) #takes factor combinations. (if exist)
  
  
  for (col_name in names(factors)) {
    factor_object <- factors[[col_name]] 
    assign(col_name, factor_object)  
  }
  
  
  
  # Get entrez ids # do not change
  if (data_source == "ncbi") { 
    keytype <- "SYMBOL"
  } else if (data_source == "ensembl") { 
    keytype <- "ENSEMBL"
    rownames(exp_table) <- sub("\\.\\d+$", "", rownames(exp_table))
  } else {
    stop("Invalid data_source. Use 'ncbi' or 'ensembl'")
  }
  
  # Check if the specified keytype is valid # do not change
  valid_keys <- keys(org.Hs.eg.db, keytype = keytype)
  if (is.null(valid_keys)) {
    stop(paste("Invalid keytype:", keytype))
  }
  
  tryCatch({
    entrez_ids <- mapIds(org.Hs.eg.db,
                         keys = row.names(exp_table),
                         keytype = keytype,
                         column = "ENTREZID")
  }, error = function(e) {
    stop("Error while mapping IDs: ", conditionMessage(e))
  })
  
  # Remove duplicates # do not change
  symbol_id  <- data.frame(SYMBOL = names(entrez_ids), ENTREZID = as.character(entrez_ids), stringsAsFactors = FALSE)
  
  symbol_id <- symbol_id[!duplicated(symbol_id$SYMBOL),]  
  
  # Create dgelist object 
  y <- DGEList(exp_table, genes = symbol_id, group = grp)
  
  # Filter by expression (Remove genes having low counts)
  keep <- filterByExpr(y)
  yfiltered <- y[keep, , keep.lib.sizes=FALSE]
  
  # Write raw filtered counts
  raw_filtered_counts_table <- file.path(output_dir, "raw_filtered_counts_table.txt")
  
  write.table(as_tibble(yfiltered$counts, rownames = "GeneID"), file = raw_filtered_counts_table, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Calculate normalization factors
  
  yfiltered <- calcNormFactors(yfiltered, method = "TMM")
  
  # Write samples table
  samples_table <- file.path(output_dir, "samples_table.txt")
  
  write.table(as_tibble(yfiltered$samples, rownames = "Samples"), file = samples_table, sep = '\t', quote = FALSE, row.names = FALSE)
  
  #Create MD plots of each sample 
  
  md_plot <- file.path(output_dir, "MDPlots.jpg")
  
  jpeg(md_plot, width = 750, height = 750, quality = 100)
  par(mfrow = c(ceiling(length(sample_names)/2), 2), mar = c(3, 3, 3, 3))
  sapply(1:length(sample_names), function(i) {
    plotMD(cpm(yfiltered, log=TRUE), column = i, 
           xlab = "Log CPM", ylab = "Log-fold-change") 
    abline(h=0, col="red", lty=2, lwd=2)
  })
  dev.off()
  
  #Explore the data with principal components
  
  pca_plot <- file.path(output_dir, "PCA_plot.jpg")
  
  jpeg(pca_plot, width = 750, height = 500, quality = 100)
  points <- seq(1,length(levels(grp)))
  colors <- seq(1,length(levels(grp)))
  plotMDS(yfiltered, col=colors[grp], pch=points[grp])
  legend("top", legend=levels(grp), pch=points, col=colors, ncol=2)
  dev.off()
  
  # Design for single factorial model/multi factorial model
  
  if (factor_count == 1) {
    single_design <- model.matrix(~ 0 + grp)
    colnames(single_design) <- levels(grp)
    yfiltered <- estimateDisp(yfiltered, single_design, robust = TRUE)
    
    dispersion_plot <- file.path(output_dir, "Dispersion_plot.jpg")
    jpeg(dispersion_plot, width = 750, height = 500, quality = 100)
    plotBCV(yfiltered)
    dev.off()
    
    single_fit <- glmQLFit(yfiltered, single_design, robust = TRUE)
    mean_ql_disp <- file.path(output_dir, "Fitted_mean_ql_dispersion_plot.jpg")
    jpeg(mean_ql_disp, width = 750, height = 500, quality = 100)
    plotQLDisp(single_fit)
    dev.off()
    
    contrast_combinations <- combn(levels(grp), 2, simplify = FALSE)
    
    # Create tables and plot MD graphs for single factor (all comparisons between factor levels)
    
    for (combo in contrast_combinations) {
      treatment1 <- combo[1]
      treatment2 <- combo[2]
      
      contrast_name <- paste(treatment1, "-", treatment2, sep = "")
      
      cmd <- paste("con <- makeContrasts(", contrast_name, ", levels=single_design)", sep ='"')
      eval(parse(text = cmd))
      
      qlf <- glmQLFTest(single_fit, contrast = con)
      
      #top_tags_qlf <- topTags(qlf, n = nrow(qlf$table))
      
      tr <- glmTreat(single_fit, contrast = con, lfc = log2(cutoff_lfc))
      
      top_tags_tr <- topTags(tr, n = nrow(tr$table))
      
      #qlf_tp <- file.path(output_dir, paste(contrast_name, "result_qlftest.txt", sep = "_"))
      tr_tp <- file.path(output_dir, paste(contrast_name, "_result_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      
      #write.table(as_tibble(top_tags_qlf)[1], file = qlf_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(as_tibble(top_tags_tr)[1], file = tr_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      
      #summary_qlf <- summary(decideTests(qlf))
      summary_tr <- summary(decideTests(tr, p.value = cutoff_p))
      
      #sum_qlf <- file.path(output_dir, paste(contrast_name, "summary_qlf_test.txt", sep = "_"))
      sum_tr <- file.path(output_dir, paste(contrast_name, "_summary_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      
      #write.table(summary_qlf, file = sum_qlf, quote = FALSE, sep = "\t")
      write.table(summary_tr, file = sum_tr, quote = FALSE, sep = "\t")
      
      
      #qlf_plt <- file.path(output_dir, paste(contrast_name, "MDPlot.jpg", sep = "_"))
      tr_plt <- file.path(output_dir, paste(contrast_name, "_MDPlot_w", cutoff_lfc, "cutoff.jpg", sep = ""))
      
      #jpeg(qlf_plt, width = 750, height=500, quality = 100)
      #plotMD(qlf)
      #abline(h=c(-1,1), col = "blue")
      #dev.off()
      
      jpeg(tr_plt, width = 750, height=500, quality = 100)
      plotMD(tr)
      abline(h=c(-cutoff_lfc, cutoff_lfc), col = "blue")
      dev.off()
      
      most_de_genes = rownames(topTags(qlf, n= 20))
      an_df <- data.frame(treatments = grp[grp%in%combo])
      row.names(an_df) <- row.names(yfiltered$samples[which(grp %in% grp[grp %in% combo]),])
      heatmap_plt <- file.path(output_dir, paste(contrast_name, "heatmap_log2_transformed.jpg", sep="_"))
      
      jpeg(heatmap_plt, width = 750, height=500, quality = 100)
      
      pheatmap(cpm(yfiltered, log = TRUE)[most_de_genes,which(grp %in% grp[grp %in% combo])], cluster_rows=FALSE, show_rownames=TRUE,
              cluster_cols=FALSE, annotation_col=an_df)
      dev.off()
      
      
      volcano_data <- data.frame(log2FoldChange = top_tags_tr$table$logFC, pvalue_adj = top_tags_tr$table$FDR)
      volcano_plt <- file.path(output_dir, paste(contrast_name, "volcano_plot.jpg", sep = "_"))

      volcano_plot <- EnhancedVolcano(volcano_data, lab=row.names(top_tags_tr$table), x="log2FoldChange", y="pvalue_adj", 
                      pCutoff = cutoff_p, FCcutoff = cutoff_lfc, pointSize = 3)
      ggsave(filename = volcano_plt, plot = volcano_plot, width = 7.5, height = 10, dpi = 200)

      #Gene ontologies for each combination of comparisons.
      
      go <- goana(tr, species ="Hs", geneid = tr$genes$ENTREZID, FDR = cutoff_p)
      godata <- topGO(go, n=30, truncate=30)
      
      godata_t <- file.path(output_dir, paste(contrast_name, "GO_results.txt", sep = "_"))
      
      write.table(as_tibble(godata, rownames = "GO_ids"), file = godata_t, quote = FALSE, sep = "\t", row.names = FALSE)
      
    }
    
    
  } else {
    if (interaction == "True") {
      formula <- as.formula(paste("~", paste(names(factors), collapse = "*")))
    } else {
      formula <- as.formula(paste("~", paste(names(factors), collapse = "+")))
    }
    
    multi_design <- model.matrix(formula)
    single_design <- model.matrix(~ 0 + grp)
    colnames(single_design) <- levels(grp)
    rownames(multi_design) <- colnames(yfiltered)
    yfiltered <- estimateDisp(yfiltered, single_design, robust = TRUE)
    yfiltered_multi <- estimateDisp(yfiltered, multi_design)
    
    dispersion_plot <- file.path(output_dir, "Dispersion_plot.jpg")
    jpeg(dispersion_plot, width = 750, height = 500, quality = 100)
    plotBCV(yfiltered_multi)
    dev.off()
    
    single_fit <- glmQLFit(yfiltered, single_design, robust = TRUE)
    mean_ql_disp <- file.path(output_dir, "Fitted_mean_ql_dispersion_plot_single.jpg")
    jpeg(mean_ql_disp, width = 750, height = 500, quality = 100)
    plotQLDisp(single_fit)
    dev.off()
    
    multi_fit <- glmQLFit(yfiltered_multi, multi_design, robust = TRUE)
    mean_ql_disp <- file.path(output_dir, "Fitted_mean_ql_dispersion_plot_multi.jpg")
    jpeg(mean_ql_disp, width = 750, height = 500, quality = 100)
    plotQLDisp(multi_fit)
    dev.off()
    
    contrast_combinations <- combn(levels(grp), 2, simplify = FALSE)
    
    # Create tables and plot MD graphs for of all multifactor combinations (comparison between each combination)
    
    for (combo in contrast_combinations) {
      treatment1 <- combo[1]
      treatment2 <- combo[2]
      
      contrast_name <- paste(treatment1, "-", treatment2, sep = "")
      
      cmd <- paste("con <- makeContrasts(", contrast_name, ", levels=single_design)", sep ='"')
      eval(parse(text = cmd))
      
      qlf <- glmQLFTest(single_fit, contrast = con)
      
      #top_tags_qlf <- topTags(qlf, n = nrow(qlf$table))
      
      tr <- glmTreat(single_fit, contrast = con, lfc = log2(cutoff_lfc))
      
      top_tags_tr <- topTags(tr, n = nrow(tr$table))
      
      #qlf_tp <- file.path(output_dir, paste(contrast_name, "result_qlftest.txt", sep = "_"))
      tr_tp <- file.path(output_dir, paste(contrast_name, "_result_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      
      #write.table(as_tibble(top_tags_qlf)[1], file = qlf_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      write.table(as_tibble(top_tags_tr)[1], file = tr_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      
      #summary_qlf <- summary(decideTests(qlf))
      summary_tr <- summary(decideTests(tr, p.value = cutoff_p))
      
      #sum_qlf <- file.path(output_dir, paste(contrast_name, "summary_qlf_test.txt", sep = "_"))
      sum_tr <- file.path(output_dir, paste(contrast_name, "_summary_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      
      #write.table(summary_qlf, file = sum_qlf, quote = FALSE, sep = "\t")
      write.table(summary_tr, file = sum_tr, quote = FALSE, sep = "\t")
      
      
      #qlf_plt <- file.path(output_dir, paste(contrast_name, "MDPlot.jpg", sep = "_"))
      tr_plt <- file.path(output_dir, paste(contrast_name, "_MDPlot_w", cutoff_lfc, "cutoff.jpg", sep = ""))
      
      #jpeg(qlf_plt, width = 750, height=500, quality = 100)
      #plotMD(qlf)
      #abline(h=c(-1, 1), col = "blue")
      #dev.off()
      
      jpeg(tr_plt, width = 750, height=500, quality = 100)
      plotMD(tr)
      abline(h=c(-cutoff_lfc, cutoff_lfc), col = "blue")
      dev.off()

      most_de_genes = rownames(topTags(qlf, n= 20))
      an_df <- data.frame(treatments = grp[grp%in%combo])
      row.names(an_df) <- row.names(yfiltered$samples[which(grp %in% grp[grp %in% combo]),])
      heatmap_plt <- file.path(output_dir, paste(contrast_name, "heatmap_log2_transformed.jpg", sep="_"))
      
      jpeg(heatmap_plt, width = 750, height=500, quality = 100)
      
      pheatmap(cpm(yfiltered, log = TRUE)[most_de_genes,which(grp %in% grp[grp %in% combo])], cluster_rows=FALSE, show_rownames=TRUE,
              cluster_cols=FALSE, annotation_col=an_df)
      dev.off()
      
      volcano_data <- data.frame(log2FoldChange = top_tags_tr$table$logFC, pvalue_adj = top_tags_tr$table$FDR)
      volcano_plt <- file.path(output_dir, paste(contrast_name, "volcano_plot.jpg", sep = "_"))
      
      volcano_plot <- EnhancedVolcano(volcano_data, lab=row.names(top_tags_tr$table), x="log2FoldChange", y="pvalue_adj", 
                                      pCutoff = cutoff_p, FCcutoff = cutoff_lfc, pointSize = 3)
      ggsave(filename = volcano_plt, plot = volcano_plot, width = 7.5, height = 10, dpi = 200)
      #Gene ontologies for each combination of comparisons.
      
      go <- goana(tr, species ="Hs", geneid = qlf$genes$ENTREZID, FDR = cutoff_p)
      godata <- topGO(go, n=30, truncate=30)
      
      godata_t <- file.path(output_dir, paste(contrast_name, "GO_results.txt", sep = "_"))
      
      write.table(as_tibble(godata, rownames = "GO_ids"), file = godata_t, quote = FALSE, sep = "\t", row.names = FALSE)
      
    }
    
    for (i in 2:length(colnames(multi_design))){
      qlf <- glmQLFTest(multi_fit, coef = i)
      tr <- glmTreat(multi_fit, coef = i, lfc = log(cutoff_lfc))
      
      #top_tags_qlf <- topTags(qlf, n = nrow(qlf$table))
      #qlf_tp <- file.path(output_dir, paste(colnames(multi_design)[i], "result_qlftest.txt", sep = "_"))
      #write.table(as_tibble(top_tags_qlf)[1], file = qlf_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      
      top_tags_tr <- topTags(tr, n = nrow(tr$table))
      tr_tp <- file.path(output_dir, paste(colnames(multi_design)[i], "_result_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      write.table(as_tibble(top_tags_tr)[1], file = tr_tp, sep = "\t", quote = FALSE, row.names = FALSE)
      
      #summary_qlf <- summary(decideTests(qlf))
      summary_tr <-summary(decideTests(tr, p.value=cutoff_p))
      
      sum_tr <- file.path(output_dir, paste(colnames(multi_design)[i], "_summary_qlf_test_", "w", cutoff_lfc ,"cutoff.txt", sep = ""))
      #sum_qlf <- file.path(output_dir, paste(colnames(multi_design)[i], "summary_qlf_test.txt", sep = "_"))
      
      #write.table(summary_qlf, file = sum_qlf, quote = FALSE, sep = "\t")
      write.table(summary_tr, file = sum_tr, quote = FALSE, sep = "\t")
      
      #qlf_plt <- file.path(output_dir, paste(colnames(multi_design)[i], "MDPlot.jpg", sep = "_"))
      
      #jpeg(qlf_plt, width = 750, height=500, quality = 100)
      #plotMD(qlf)
      #abline(h=c(-1, 1), col = "blue")
      #dev.off()
      
      
      tr_plt <- file.path(output_dir, paste(colnames(multi_design)[i], "_MDPlot_w", cutoff_lfc, "cutoff.jpg", sep = ""))
      
      jpeg(tr_plt, width = 750, height=500, quality = 100)
      plotMD(tr)
      abline(h=c(-cutoff_lfc, cutoff_lfc), col = "blue")
      dev.off()

      volcano_data <- data.frame(log2FoldChange = top_tags_tr$table$logFC, pvalue_adj = top_tags_tr$table$FDR)
      volcano_plt <- file.path(output_dir, paste(colnames(multi_design)[i], "volcano_plot.jpg", sep = "_"))
      
      volcano_plot <- EnhancedVolcano(volcano_data, lab=row.names(top_tags_tr$table), x="log2FoldChange", y="pvalue_adj", 
                                      pCutoff = cutoff_p, FCcutoff = cutoff_lfc, pointSize = 3)
      ggsave(filename = volcano_plt, plot = volcano_plot, width = 7.5, height = 10, dpi = 200)
      
      
      go <- goana(tr, species ="Hs", geneid = qlf$genes$ENTREZID, FDR = cutoff_p)
      godata <- topGO(go, n=30, truncate=30)
      
      godata_t <- file.path(output_dir, paste(colnames(multi_design)[i], "GO_results.txt", sep = "_"))
      
      write.table(as_tibble(godata, rownames = "GO_ids"), file = godata_t, quote = FALSE, sep = "\t", row.names = FALSE)

      most_de_genes = rownames(topTags(tr, n= 20))
      
      if (interaction == TRUE) {
        if (i < length(colnames(multi_design))) {
          tryCatch({
            an_df <- data.frame(treatments = factors[, i - 1])
          }, error = function(e) {
            cat("Error occurred in creating annotation df", i, ":", conditionMessage(e), "\n")
          })    
          
          row.names(an_df) <- row.names(yfiltered$samples)
          heatmap_plt <- file.path(output_dir, paste(colnames(multi_design)[i], "heatmap_log2_transformed.jpg", sep="_"))
          
          tryCatch({
            jpeg(heatmap_plt, width = 750, height = 500, quality = 100)
            
            pheatmap(cpm(yfiltered, log = TRUE)[most_de_genes, ], cluster_rows = FALSE, show_rownames = TRUE,
                    cluster_cols = FALSE, annotation_col = an_df)
            
            dev.off()
          }, error = function(e) {
            cat("Error occurred in plotting for iteration", i, ":", conditionMessage(e), "\n")
          })
          
        }
      } else {
        tryCatch({
          an_df <- data.frame(treatments = factors[, i - 1])
        }, error = function(e) {
          cat("Error occurred in creating annotation df", i, ":", conditionMessage(e), "\n")
        })    
        
        row.names(an_df) <- row.names(yfiltered$samples)
        heatmap_plt <- file.path(output_dir, paste(colnames(multi_design)[i], "heatmap_log2_transformed.jpg", sep="_"))
        
        tryCatch({
          jpeg(heatmap_plt, width = 750, height = 500, quality = 100)
          
          pheatmap(cpm(yfiltered, log = TRUE)[most_de_genes, ], cluster_rows = FALSE, show_rownames = TRUE,
                  cluster_cols = FALSE, annotation_col = an_df)
          
          dev.off()
        }, error = function(e) {
          cat("Error occurred in plotting for iteration", i, ":", conditionMessage(e), "\n")
        })
        
      }

    }
    
  }
  
}

parser <- ArgumentParser(description = "From Gene Level Counts to Identifying Differentially Expressed Genes")

parser$add_argument("--sample_info_file", help = "File containing sample information (names, treatments and directories)")
parser$add_argument("--output_dir", help = "Output Directory")
parser$add_argument("--data_source", help = "Data source (Can be either ncbi or ensembl)")
parser$add_argument("--interaction", help = "Find out interaction between the factors if design is multifactorial (TRUE or FALSE)")
parser$add_argument("--lfc_cutoff", help = "Program makes calculations and creates charts with respect to this log fold change value")
parser$add_argument("--p_cutoff", help = "Cutoff p value - i.e 0.05, or 0.01")

args <- parser$parse_args()


DiffexpWedgeR(args$sample_info_file, args$output_dir, args$data_source, args$interaction, args$lfc_cutoff, args$p_cutoff)
  