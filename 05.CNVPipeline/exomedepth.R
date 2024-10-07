library(ExomeDepth)
library(doParallel)
library(foreach)

setwd("/home/genwork2/Mert/CNVpip/newtrial")

bamfiles <- list.files("/home/genwork2/Mert/CNVpip/newtrial", pattern = ".bam$")

segments <- read.table("/home/genwork2/Mert/CNVpip/hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.bed", sep = "\t", as.is = TRUE, fill = TRUE)


colnames(segments) = c("chromosome","start","end", "name")

fasta <- "/home/genwork2/Mert/hg38.fa"

num_cores <- 40 

cstart_time <- Sys.time()

registerDoParallel(cores = num_cores)

getCountsFromBAM <- function(bamfile, segments) {
  getBamCounts(bed.frame = segments, bam.files = bamfile, include.chr = FALSE, referenceFasta = fasta)
}

results <- foreach(bamfile = bamfiles, .combine = c) %dopar% {
  getCountsFromBAM(bamfile, segments)
}

stopImplicitCluster()

cend_time <- Sys.time()


my.counts <- data.frame(results)

my.counts <- my.counts[!duplicated(as.list(my.counts))]

count_matrix <- as.matrix(my.counts[,!names(my.counts) %in% c("chromosome", "start", "end", "exon", "GC")])


samples <- colnames(count_matrix)
samples <- sub("\\.bam$", "", samples)

nsamples <- ncol(count_matrix)

pstart_time <- Sys.time()

registerDoParallel(cores = num_cores)

iterate <- function(i) {
 
  my.choice <- select.reference.set(test.counts = count_matrix[,i], reference.counts = count_matrix[,-i], bin.length = (my.counts$end - my.counts$start)/1000,n.bins.reduced = 10000)

  my.reference.selected <- apply(X = count_matrix[, my.choice$reference.choice, drop = FALSE], MAR = 1, FUN = sum)

  all.exons <- new('ExomeDepth', test = count_matrix[,i], reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')

  all.exons <- CallCNVs(x = all.exons, transition.probability = 10^-4, chromosome = my.counts$chromosome, start = my.counts$start, end = my.counts$end, name = my.counts$exon, expected.CNV.length = 10000000)

  output_file <- paste0(samples[i], "_cnvcalls.txt")

  write.table(file = output_file, x = all.exons@CNV.calls, row.names = FALSE, sep = "\t")
}

indices = 1:nsamples

foreach(i = indices, .packages = c("ExomeDepth"), .combine = c) %dopar% {
  iterate(i)
}

stopImplicitCluster()

pend_time <- Sys.time()

print(paste("Count time:", (cend_time - cstart_time)))
print(paste("CNV calling time:", (pend_time-pstart_time)))