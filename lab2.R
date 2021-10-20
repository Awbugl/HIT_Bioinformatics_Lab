data <- read.table(
    "C:/Users/hit/Desktop/refGene.hg19.sorted.bed",
    header=F,    sep = "\t",
    col.names = c("chrom","chromStart","chromEnd","name","score","strand","thickStart",
        "thickEnd","geneName","blockCount","blockSizes","blockStarts") # http://genome.ucsc.edu/FAQ/FAQformat.html
  )

result <- head(data[which(data$chrom == "chr22"),],20) # data is sorted, get top20 


promoterFunc <- function(strand,chromStart,chromEnd){ # 5' -> 3' 
    ifelse(strand == "+",    chromStart - 2000,    chromEnd + 2000)
}

result$promoterStart <- mapply(promoterFunc,
    result$strand,result$chromStart,result$chromEnd)


exonFunc <- function(strand,blockStarts,chromStart,blockSizes){
    paste(
      mapply("sum",
        chromStart,mapply(as.numeric,strsplit(blockStarts, ",")), # chromStart + blockStarts
        
        ifelse(strand == "-",0,
               mapply(as.numeric,strsplit(blockSizes, ","))) # if strand == "+" then plus blockSizes
        
      ),collapse =",") # ToString()
}

result$exon_3_pos <- mapply(exonFunc,
    result$strand,result$blockStarts,result$chromStart,result$blockSizes)


final_result <- result[,c('name','strand','promoterStart','exon_3_pos')]

View(final_result)
