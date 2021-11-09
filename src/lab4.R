expr <- read.table("C:/Users/hit/Desktop/lab4/GSE66597.expr.txt",sep = "\t",header = T)
sample_info <- read.table("C:/Users/hit/Desktop/lab4/GSE66597.sample_info.txt",sep = "\t")

scale_data <- scale(log(expr[-(1:3)]))

#' @title findDiffGene 
#' @description find differentiated expression genes (p-value < 0.01 & abs(log2(FC)) >= 1)
#' @param scale_data the scaled data
findDiffGene <- function(scale_data,sample_info)
{
   con_data <- sample_info[which(sample_info$V2 == "con"),]$V1
   hm_data <- sample_info[which(sample_info$V2 == "hm"),]$V1
      
   #' @title getPValues 
   #' @description Calc P-values
   #' @param scale_data the scaled data
   getPValues <- function(scale_data)
   {
      apply(scale_data, 1, function(scale_data_item) t.test(scale_data_item[con_data], scale_data_item[hm_data])$p.value)
   }
   
   #' @title getlog2FC 
   #' @description Calc log2(FC)
   #' @param scale_data the scaled data
   #' @param sample_info the sample_info file
   getlog2FC <- function(scale_data,sample_info)
   {
      log2(abs(apply(scale_data[,hm_data], 1, mean)/apply(scale_data[,con_data], 1, mean)))
   }
   
   p_values <- getPValues(scale_data)
   
   log2FC <- getlog2FC(scale_data,sample_info)
   
   diff_gene <- cbind(expr[1:3], p_values, log2FC)
   
   diff_gene$expr = "Normal"
   diff_gene[which(diff_gene$log2FC <= -1),"expr"] = "Down"
   diff_gene[which(diff_gene$log2FC >= 1) ,"expr"] = "Up"
   diff_gene[which(diff_gene$p_values >= 0.01),"expr"] = "Not significant"
    
   diff_gene
}


#' @title plotDiffGenes 
#' @description provided by weiweiwei010119
#' @param diff_gene the differentiated expression genes data
plotDiffGenes <- function(diff_gene)
{
   library(ggplot2)
   p <- qplot(x = diff_gene$log2FC,y =-log10(diff_gene$p_values),
              xlab = "log2FC",ylab ="-log10(p.values)",color = factor(diff_gene$expr))
   #Normally, -log10(FDR) is used as the ordinate, 
   #but in this question, use the p-value as the ordinate (provided by weiweiwei010119)
   p <- p + geom_vline(xintercept = c(-1, 1),lty = 2,color ="grey11") 
   p <- p + geom_hline(yintercept = 2,lty = 2,color ="grey11") 
   p <- p + theme_bw() + 
      theme(panel.background = element_rect(color = "black", size = 1, fill ="white")
            ,panel.grid = element_blank()) 
   p
}

diff_gene <- findDiffGene(scale_data,sample_info)
plotDiffGenes(diff_gene)
  
significant_diff_gene <- diff_gene[-which(diff_gene$expr == "Not significant"),]
View(significant_diff_gene)
  
plot(hclust(dist(t(scale_data[which(expr$ID %in% significant_diff_gene$ID), ]))))
