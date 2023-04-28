#' Normalize count data
#'
#' @param count raw data
#' @param method name which normalization method to use
#'
#' @return normalized data
#' @export

normcount <- function(count,method){
  if (method=='cpm'){
    count1<-t(t(count)/colSums(count)*1e6)
  } else if (method=='logcpm') {
    count1<-log2(t(t(count+0.5)/(colSums(count)+1)*1e6))
  } else if (method=='countNorm') {
    count1 <- t(t(count)/colSums(count))*median(colSums(count))
  } else if (method=='quantile') {
    count1 <- NOISeq::uqua(count)
  } else if (method=='tmm'){
    count1=NOISeq::tmm(count)
  } else if (method=='tpm'){
    DGEobj.utils::convertCounts(as.matrix(count),unit='TPM',
                                geneLength= rowMeans(count),log = FALSE,normalize = "none")
  } else { # no normalize
    count1=as.matrix(count)
  }
}


################################################################################################
#' identify cell types in each domain
#'
#' @param ctData estimated coefficients for each cell type within each domain
#' @param ilay ith domain
#' @param cutoff the minimum percent of nonzero entries, default value is 0
#' @param offset the value add on the threshold to binarize coefficient from lasso regression, if not provide by user, then use the mean of all coefficients
#'
#' @return the cell type names that pass the filtering
#' @export
#'

ctcut1=function(ctData,ilay,cutoff=0,offset=NULL){
  rep1_layer1=t(ctData[[ilay]])
  aa=abs(rep1_layer1)
  if (is.null(offset)){
    offset=mean(colMeans(aa))
  }
  b=EBImage::thresh(aa, w=ceiling(mean(aa)*0.6), h=1, offset=offset)
  b1=apply(b,2,function(x) sum(x!=0)/length(x))  # count nonzero
  names(which(b1>cutoff))
}

