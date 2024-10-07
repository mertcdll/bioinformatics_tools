suppressPackageStartupMessages(library(tximport))
suppressPackageStartupMessages(library(ensembldb))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(argparse))

processRNASeq <- function(id_file, abundance_file , output_dir, sample_name) {
  
  transcript_data <- read.table(abundance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  tx2gene_df <- read.table(id_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  names(abundance_file) <- sample_name

  txi_list <- tximport(abundance_file,
                       type = "kallisto",
                       tx2gene = tx2gene_df, ignoreTxVersion = FALSE)
  
  
  output_abundance <- file.path(output_dir, "gene_level_abundance.txt")
  output_counts <- file.path(output_dir, "gene_level_counts.txt")
  output_lengths <- file.path(output_dir, "gene_lengths")
  
  write.table(txi_list$abundance, file = output_abundance, sep = "\t", quote = FALSE)
  write.table(txi_list$counts, file = output_counts, sep = "\t", quote = FALSE)
  write.table(txi_list$length, file = output_lengths, sep = "\t", quote = FALSE)
}


parser <- ArgumentParser(description = "Transcript to Gene Level Counts")

parser$add_argument("--id_file", help = "Path to ID file")
parser$add_argument("--abundance_file", help = "Path to abundance file")
parser$add_argument("--output_dir", help = "Output directory")
parser$add_argument("--sample_name", help = "Name of the Sample")

args <- parser$parse_args()


processRNASeq(args$id_file, args$abundance_file , args$output_dir, args$sample_name)





