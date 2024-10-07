suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(statmod))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(GO.db)) 
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(argparse))
suppressPackageStartupMessages(library(dplyr))


DiffexpWedgeR <- function(sample_info_file, output_dir, data_source) {
  
  sample_info <- read.table(sample_info_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  sample_names <- sample_info$names
  count_dirs <- sample_info$count_dirs
  
  
  exp_dict <- {}
  
  # Create raw counts dictionary 
  for (i in 1:length(sample_names)) {
    count_table <- read.table(count_dirs[i], header = TRUE, row.names = 1, sep = "\t")
    exp_dict[[sample_names[i]]] <- count_table
  }
  
  # Convert it to table
  exp_table <- data.frame(exp_dict)
  

  # Convert floats to integers
  for (col in names(exp_table)) {   
    if (is.numeric(exp_table[[col]])) {
      exp_table[[col]] <- as.integer(exp_table[[col]])
    }
  } 
  
  
  raw_counts_table <- file.path(output_dir, "raw_counts_table.txt")
  
  # write raw unfiltered counts
  write.table(as_tibble(exp_table, rownames = "GeneID"), file = raw_counts_table, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Get groups
  group <- as.factor(sample_info$treatments)
  
  # Get entrez ids
  if (data_source == "ncbi") { 
    keytype <- "SYMBOL"
  } else if (data_source == "ensembl") { 
    keytype <- "ENSEMBL"
    rownames(exp_table) <- sub("\\.\\d+$", "", rownames(exp_table))
  } else {
    stop("Invalid data_source. Use 'ncbi' or 'ensembl'")
  }

  # Check if the specified keytype is valid
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
  
  # Remove duplicates
  symbol_id  <- data.frame(SYMBOL = names(entrez_ids), ENTREZID = as.character(entrez_ids), stringsAsFactors = FALSE)
  
  symbol_id <- symbol_id[!duplicated(symbol_id$SYMBOL),]  
  
  # Create dgelist object
  y <- DGEList(exp_table, group = group, genes = symbol_id)
  
  # Filter by expression (Remove genes having low counts)
  keep <- filterByExpr(y)
  yfiltered <- y[keep, , keep.lib.sizes=FALSE]
  
  # Write raw filtered counts
  raw_filtered_counts_table <- file.path(output_dir, "raw_filtered_counts_table.txt")
  
  write.table(as_tibble(yfiltered$counts, rownames = "GeneID"), file = raw_filtered_counts_table, sep = '\t', quote = FALSE, row.names = FALSE)
  
  # Calculate normalization factors
  
  yfiltered <- calcNormFactors(yfiltered)
  
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
  points <- c(0,1)
  colors <- rep(c("blue", "red"), times = ceiling(length(sample_names)/2))
  plotMDS(yfiltered, col=colors[group], pch=points[group])
  legend("topleft", legend=levels(group), pch=points, col=colors, ncol=2)
  dev.off()
  
  # Design for single factorial model
  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  # Dispersion estimation
  
  yfiltered <- estimateDisp(yfiltered, design, robust = TRUE)
  
  dispersion_plot <- file.path(output_dir, "Dispersion_plot.jpg")
  
  jpeg(dispersion_plot, width = 750, height = 500, quality = 100)
  plotBCV(yfiltered)
  dev.off()
  
  # Differential expression
  
  fit <- glmQLFit(yfiltered, design, robust=TRUE)
  
  mean_ql_disp <- file.path(output_dir, "Fitted_mean_ql_dispersion_plot.jpg")
  
  jpeg(mean_ql_disp, width = 750, height = 500, quality = 100)
  plotQLDisp(fit)
  dev.off()
  
  # Creating qlf and qlf-tr table/plots for each contrast combinations: there can me more than 2 groups.
  
  contrast_combinations <- combn(levels(group), 2, simplify = FALSE)
  
  # Create tables and plot MD graphs
  
  for (combo in contrast_combinations) {
    treatment1 <- combo[1]
    treatment2 <- combo[2]
    
    contrast_name <- paste(treatment1, "-", treatment2, sep = "")
    
    cmd <- paste("con <- makeContrasts(", contrast_name, ", levels =design)", sep ='"')
    eval(parse(text = cmd))
    
    qlf <- glmQLFTest(fit, contrast = con)
    
    top_tags_qlf <- topTags(qlf, n = nrow(qlf$table))
    
    tr <- glmTreat(fit, contrast = con, lfc = log2(1.2))
    
    top_tags_tr <- topTags(tr, n = nrow(tr$table))
    
    qlf_tp <- file.path(output_dir, paste(contrast_name, "result_qlftest.txt", sep = "_"))
    tr_tp <- file.path(output_dir, paste(contrast_name, "result_qlf_test_w1.2cutoff.txt", sep = "_"))
    
    write.table(as_tibble(top_tags_qlf)[1], file = qlf_tp, sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(as_tibble(top_tags_tr)[1], file = tr_tp, sep = "\t", quote = FALSE, row.names = FALSE)
    
    summary_qlf <- summary(decideTests(qlf))
    summary_tr <- summary(decideTests(tr))
    
    sum_qlf <- file.path(output_dir, paste(contrast_name, "summary_qlf_test.txt", sep = "_"))
    sum_tr <- file.path(output_dir, paste(contrast_name, "summary_qlf_test_w1.2cutoff.txt", sep = "_"))
    
    write.table(summary_qlf, file = sum_qlf, quote = FALSE, sep = "\t")
    write.table(summary_tr, file = sum_tr, quote = FALSE, sep = "\t")
    
    
    qlf_plt <- file.path(output_dir, paste(contrast_name, "MDPlot.jpg", sep = "_"))
    tr_plt <- file.path(output_dir, paste(contrast_name, "MDPlot_wcutoff.jpg", sep = "_"))
    
    jpeg(qlf_plt, width = 750, height=500, quality = 100)
    plotMD(qlf)
    dev.off()
    
    jpeg(tr_plt, width = 750, height=500, quality = 100)
    plotMD(tr)
    dev.off()
    
    #Gene ontologies for each combination of comparisons.
    
    go <- goana(qlf, species ="Hs", geneid = qlf$genes$ENTREZID, FDR = 0.05)
    godata <- topGO(go, n=30, truncate=30)
    
    godata_t <- file.path(output_dir, paste(contrast_name, "GO_results.txt", "_"))
    
    write.table(as_tibble(godata, rownames = "GO_ids"), file = godata_t, quote = FALSE, sep = "\t", row.names = FALSE)
    
  }
  
  
}

parser <- ArgumentParser(description = "From Gene Level Counts to Identifying Differentially Expressed Genes")

parser$add_argument("--sample_info_file", help = "File containing sample information (names, treatments and directories)")
parser$add_argument("--output_dir", help = "Output Directory")
parser$add_argument("--data_source", help = "Data source (Can be either ncbi or ensembl)")


args <- parser$parse_args()


DiffexpWedgeR(args$sample_info_file, args$output_dir, args$data_source)

## Add heatmaps...


