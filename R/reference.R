
#' construct cell type reference from scRNA-seq data
#'
#' @param input_data an expressionSet class contains single cell gene expression and its metadata which contains annoted cell type
#' @param ct_var variable name for cell type column
#' @param sample_var variable name for sample column
#'
#' @return a matrix contains gene expression for each cell type
#' @export
#'

scRefer <- function(input_data, ct_var, sample_var) {
  library(Biobase)
  data_subset <- input_data[rowSums(exprs(input_data)) > 0,]
  counts <- exprs(data_subset)
  metad=input_data@phenoData@data
  sample_id <- unique(as.character(data_subset@phenoData@data[, sample_var]))
  ct_id <- unique((as.factor(data_subset@phenoData@data[, ct_var])))
  avg_exp=list()
  for (i in sample_id){
    avg_exp[[i]] <- sapply(ct_id, function(j) {
      a=metad[c(metad[,sample_var]==i & metad[,ct_var]==j),]
      b = as.matrix(counts[,rownames(a)])
      exp=apply(b, 1, sum, na.rm = TRUE) / sum(b)
      sf=sum(b) / ncol(b)
      return(list(meancount=exp,scalefac=sf))
    })
  }
  scalFac <- matrix(0, nrow=length(ct_id),ncol=length(sample_id))
  meanCount=c()
  for (i in 1:length(sample_id)){
    for (j in 1:length(ct_id)){
      scalFac[j,i]=avg_exp[[i]][[2*j]]
      mc=avg_exp[[i]][[2*j-1]]
      meanCount=cbind(meanCount,mc)
    }
  }
  rownames(scalFac) <- ct_id
  colnames(scalFac) <- sample_id
  avg_scalFac <- rowMeans(scalFac, na.rm = T)
  colnames(meanCount)=sapply(sample_id, paste, ct_id,sep='_')
  avg_scalFac2=rep(avg_scalFac,length(sample_id))

  scref <- t(t(meanCount) *avg_scalFac2)
  scref2<-matrix(0,nrow=nrow(scref),ncol=length(ct_id))
  colnames(scref2)=ct_id;rownames(scref2)=rownames(scref)
  for(i in ct_id){
    y <- scref[,grep(i, colnames(meanCount))]
    scref2[,i]=apply(y, 1, mean, na.rm = TRUE)
  }
  return(scref2)
}
