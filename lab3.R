#' @title Download_JASPAR_PFM
#' @description Download pfm file from JASPAR
#' @param JASPAR_ID the Matrix ID from JASPAR. #string
#' @return  data from JASPAR. #data.frame
Download_JASPAR_PFM <- function(JASPAR_ID){
  filename <- paste(JASPAR_ID , ".pfm",sep = "",collapse = "")
  url <- paste("https://jaspar.genereg.net/api/v1/matrix/",filename,sep = "",collapse = "")
  download.file(url, filename)
  data <- read.table(filename,skip = 1, sep = "",row.names = c("A","C","G","T"))
  file.remove(filename)
  data
}

#' @title PFM_2_PWM
#' @description Convert PFM matrix to PWM matrix.
#' @param PFMmatrix the PFM Matrix. #matrix
#' @param BgFreq the background frequency of nucleotide.
#' @return the PWM matrix.
PFM_2_PWM <- function(PFMmatrix,BgFreq = c(0.3,0.2,0.2,0.3)){
  
  PFM2PWM_foreach <- function(PFMmatrix){
    mapply(w_calcFunc,PFMmatrix,sum(PFMmatrix),BgFreq)
  }
  
  w_calcFunc <- function(data,N,BgFreq){
    sqrt_N <- sqrt(N)
    log2((data + sqrt_N/4)/(N+sqrt_N)/BgFreq)
  }
  
  PWMmatrix <- mapply(PFM2PWM_foreach,PFMmatrix)
  rownames(PWMmatrix) <- c("A","C","G","T")
  
  PWMmatrix
}

#' @title Calc_Relative_Score
#' @description Get the subseq and its relative score.
#' @param seqStr the sequence. #string
#' @param PFMmatrix the PFM Matrix. #matrix
#' @return #data.frame
Calc_Relative_Score <- function(seqStr,PWMmatrix){
  
  relativeScore <- function(PWMmatrix,absolute_score){
    min_score <- sum(apply(PWMmatrix,2,min))
    max_score <- sum(apply(PWMmatrix,2,max)) 
    (absolute_score-min_score)/(max_score-min_score)
  }
  
  seq <- strsplit(seqStr,split <- "")[[1]]
  result <- data.frame(row.names <- c("seq","score"))
  
  for(i in 1:(length(seq) - length(data))) # 28 - 19
  {
    score <- 0
    seqTmp <- seq[i:(i + length(data) - 1)] # i : i + 18
    
    for(j in 1:length(data))
      score <- score + PWMmatrix[seqTmp[j],j]
    
    result[i] <- c(paste(seqTmp,collapse <- ""),relativeScore(PWMmatrix,score))
  }
  
  t(result)
}

data <- Download_JASPAR_PFM("MA0139.1")
PWMmatrix <- PFM_2_PWM(data)
result <- Calc_Relative_Score("CCCGGGGTCCAGTAGGGGGCGCACTCAC",PWMmatrix)
View(result)
