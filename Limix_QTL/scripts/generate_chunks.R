log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(matrixStats)
library(readr)
library(dplyr)

## gene annotation ####
#Example files for testing
#infile <- "Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot_chr1_small.txt"
#nGenes <- 5
#startPos <- 0
#endOffset <- 1000000000
#outfile <- "OutGeuvadis/chunks.txt"

infile <- snakemake@input[['annotation']]
nGenes <- snakemake@params[['nGenes']]
startPos <- snakemake@params[['startPos']]
endOffset <- snakemake@params[['endOffset']]
outfile <- snakemake@output[['chunks']]


message("Reading annotation file: ", infile)
gene_anno <- read_delim(infile)

message("Process annotations, group features in chunks of ", nGenes, ",
        save as range")

gr <- c()
for(chr in sort(unique(gene_anno$chromosome))){
  
  annotationRel <- gene_anno %>%
    filter(chromosome == chr) %>%
    arrange(start, end)
  
  ## First go through the list to fix genes so they all touch.
  annotationRel$start[1] <- startPos
  for(i in 2:nrow(annotationRel)){
    if(i == nrow(annotationRel)){
      annotationRel$end[i] <- annotationRel$end[i] + endOffset
    }
    #If "overlapping" than we don't need to do anything.
    if((annotationRel$start[i] > annotationRel$end[i-1])){
      distance <- (annotationRel$start[i]-annotationRel$end[i-1])/2
      annotationRel$start[i] <- ceiling(annotationRel$start[i] - distance)
      annotationRel$end[i-1] <- floor(annotationRel$end[i-1] + distance)
    }
  }
  chunks <- seq(1, nrow(annotationRel), nGenes)
  
  # Need to add the last one as a separate entry.
  if(chunks[length(chunks)] < nrow(annotationRel)){
    chunks <- c(chunks, (nrow(annotationRel)+1))
  }
  for(i in 1:(length(chunks)-1)){
    gr <- c(gr, paste(chr,":",
                annotationRel$start[chunks[i]], "-",
                annotationRel$end[(chunks[i+1]-1)],
                sep=""))
  }
}
message("Features process and chunking finished, write to chunk file: ",
        outfile)

write_delim(as.data.frame(gr), outfile, col_names = FALSE)
message("Finished")