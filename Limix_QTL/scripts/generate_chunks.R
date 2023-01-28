library(matrixStats)
library(readr)
library(dplyr)

## gene annotation ####
infile <- "Limix_QTL/test_data/Expression/Geuvadis_CEU_Annot.txt"
nGenes <- 50
startPos <- 0
endOffset <- 1000000000
outfile <- "Limix_QTL/test_data/Expression/chunks_Geuvadis_CEU_Annot.txt"
  
#infile <- snakemake@input['annotation']
#nGenes <- snakemake@params['ngenes']
#startPos <- snakemake@params['startPos']
#endOffset <- snakemake@params['endOffset']
#outfile <- snakemake@output['chunks']

gene_anno <- read.delim(infile, as.is=TRUE)

if(file.exists(outfile)) file.remove(outfile)

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
    write_delim(data.frame(paste(chr,":",
                annotationRel$start[chunks[i]], "-",
                annotationRel$end[(chunks[i+1]-1)],
                sep="")),
                outfile,
                append=TRUE)
  }
}
