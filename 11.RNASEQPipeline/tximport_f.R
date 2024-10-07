library(tximport, quietly = TRUE)
library(ensembldb, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(argparse, quietly =TRUE)

processRNASeq <- function(gtf_file, abundance_file, source, output_dir) {

  transcript_data <- read.table(abundance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ens_db <- makeTxDbFromGFF(gtf_file, format = "gtf")
  
  if (source == "ensembl") {
    transcript_data$target_id <- gsub("\\.\\d+$", "", transcript_data$target_id)
  }
  
  transcript_ids <- transcript_data$target_id
  
  gene_ids <- mapIds(ens_db, keys = transcript_ids, column = "GENEID", keytype = "TXNAME")
  
  tx2gene_df <- data.frame(gene_ids)
  tx2gene_df$transcript_id <- rownames(tx2gene_df)
  rownames(tx2gene_df) <- NULL
  
  tx2gene_df <- tx2gene_df %>%
    select(transcript_id, gene_ids)
  
  tx2gene_df_nona <- tx2gene_df %>%
    filter(!is.na(gene_ids))
  

  if (source == "ensembl") {
    ignoreTxV <- TRUE
  } else if (source == "ncbi") {
    ignoreTxV <- FALSE
  } else {
    stop("Unknown source:", source)
  }
  
  txi_list <- tximport(abundance_file,
                       type = "kallisto",
                       tx2gene = tx2gene_df_nona, ignoreTxVersion = ignoreTxV)
  

  output_abundance <- file.path(output_dir, "gene_level_abundance.txt")
  output_counts <- file.path(output_dir, "gene_level_counts.txt")
  output_lengths <- file.path(output_dir, "gene_lengths.txt")
  
  write.table(txi_list$abundance, file = output_abundance, sep = "\t", quote = FALSE)
  write.table(txi_list$counts, file = output_counts, sep = "\t", quote = FALSE)
  write.table(txi_list$length, file = output_lengths, sep = "\t", quote = FALSE)
}



parser <- ArgumentParser(description = "Transcript to Gene Level Counts")

parser$add_argument("--gtf_file", help = "Path to GTF file")
parser$add_argument("--abundance_file", help = "Path to abundance file")
parser$add_argument("--source", help = "Data source (ensembl or ncbi)")
parser$add_argument("--output_dir", help = "Output directory")


args <- parser$parse_args()


processRNASeq(args$gtf_file, args$abundance_file, args$source, args$output_dir)


