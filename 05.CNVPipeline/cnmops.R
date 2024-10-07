
### CN MOPS


setwd("/home/genwork2/Mert/CNVpip/newtrial")

library(cn.mops)
library(dplyr)
bamfiles <- list.files("/home/genwork2/Mert/CNVpip/newtrial", pattern = ".bam$")

segments <- read.table("/home/genwork2/Mert/CNVpip/hg38_exome_v2.0.2_targets_sorted_validated.annotated_0.bed", sep = "\t", as.is = TRUE, fill =TRUE)

segments <- segments[-nrow(segments),] # exclude M

gr <- GRanges(segments[,1],IRanges(segments[,2],segments[,3]))

X <- getSegmentReadCountsFromBAM(bamfiles,GR=gr, parallel = 40)

chrY_indices <- which(seqnames(X) == "chrY") # exclude chrY

X_filtered <- X[-chrY_indices, ]

resCNMOPS <- exomecn.mops(X_filtered, segAlgorithm = "DNAcopy", parallel = 40, lowerThreshold = -0.4) 
?exomecn.mops
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPS)

CNVS <- as.data.frame(cnvs(resCNMOPS))

CNVR <- as.data.frame(cnvr(resCNMOPS))

SEGMENT <- segmentation(resCNMOPS)


write.table(CNVS, file = "cnvs.txt", sep = "\t", row.names = FALSE, quote = FALSE)


write.table(CNVR, file = "cnv_regions.txt", sep = "\t", row.names = FALSE, quote = FALSE)


write.table(SEGMENT, file = "cnv_segments.txt", sep = "\t", row.names = FALSE, quote = FALSE)