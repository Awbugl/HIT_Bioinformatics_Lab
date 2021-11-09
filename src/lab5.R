sam_data <- read.table("C:/Users/littl/Desktop/lab5/ERR1594329.chrI.sorted.slim.sam",sep = "\t",
                       col.names = c("QNAME","FLAG","RNAME","POS","MAPQ","CIGAR"))

bed_data <- read.table(
   "C:/Users/littl/Desktop/lab5/sacCer3.chrI.bed",sep = "\t",
   col.names = c("chrom","chromStart","chromEnd","name","score","strand","thickStart","thickEnd","itemRgb",
                 "blockCount","blockSizes","blockStarts")) # http://genome.ucsc.edu/FAQ/FAQformat.html)

#' @title calc_length
#' @description calc_length
#' @param sam_data
calc_length <- function(sam_data) {
   split_str <- unlist(strsplit(sam_data["CIGAR"], "[MDIS]"))
   
   index <- grep("[MDIS]",unlist(strsplit(sam_data["CIGAR"],"")),value = T)
   
   sum(mapply(function(index,split_str) ifelse(index == "M" || index == "I" ,as.numeric(split_str), 0),index, split_str))
}

sam_data$length <- apply(sam_data, 1, calc_length)

bed_data$reads_fall <- 
   mapply(
      function(chromStart,chromEnd) 
         length(which(chromStart <= sam_data[,"POS"] & chromEnd >= sam_data[,"POS"] + sam_data[,"length"])),
      bed_data$chromStart, bed_data$chromEnd
      )

bed_data$rpkm <- ((10^9)*bed_data$reads_fall)/(length(sam_data[,1])*(bed_data$chromEnd - bed_data$chromStart))

result <- bed_data[,c("name","rpkm")]

View(result)
