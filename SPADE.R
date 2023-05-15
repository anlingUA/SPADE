#' Estimate cell type coefficients at each location by LASSO regression
#'
#' @param loc data.frame contains location information
#' @param stcount spatial gene expression
#' @param scref cell type reference from single cell RNA-seq data
#' @param sign_list marker genes lists for each cell type
#' @param lasso whether to use lasso or elastic net, the default is lasso
#'
#' @return a matrix contains coefficients for each cell type at each location
#' @export
#'

CTperDom=function(loc,stcount,scref,sign_list,lasso=T){
  library(caret)
  library(glmnet)
  lapply(c(0:nlay),function(ilay){ #for each layer
    #loc contains:x/y array, x/y pixel, domain location
    locs=loc[loc$domain==ilay,][,'location']
    stcount2=stcount[,locs]
    genelist=Reduce(intersect,list(rownames(scref),rownames(stcount2),unique(sign_list)))
    spa=stcount2[genelist,] #select ct markers
    scref2=scref[genelist,] #subset ref
    print(identical(rownames(spa),rownames(scref2)))
    ct=colnames(scref)  # celltype name
    lo=colnames(spa) #location
    samples=c()
    pseudo.GCT=NULL  #gene by celltype matrix
    for (i in 1:length(lo)){
      samples=c(samples,lo[i],ct)
      pseudo=cbind(spa[,i],scref2) #combine location with ct,so loc1 ct loc2 ct....
      pseudo.GCT=cbind(pseudo.GCT,pseudo)
    }
    colnames(pseudo.GCT)=samples
    nct=ncol(scref)+1
    set.seed(2022)
    if (lasso==T){   #use lasso to select
      ct_select=lapply(c(0:(length(lo)-1)), function(id){
        df1=as.data.frame(pseudo.GCT[,c((id*nct+1):(id*nct+nct))])  # loc.i ct1 ct2....ct12
        X <- model.matrix(df1[,1]~.,data=df1)
        if (ncol(X)!=(nct-1)){
          X=X[,-c(1,2)]
        }
        Y <- df1[,1]
        # Y is same for all cause error, so add a very small number to Y
        if (quantile(Y,0.95)==0){
          Y=Y+runif(length(Y),1e-10,1e-6)
        }
        #find the amount of penalty,λ by cross-validation
        cv.lambda.lasso <- cv.glmnet(x=X, y=Y, alpha = 1)
        l.lasso.min <- cv.lambda.lasso$lambda.min
        lasso.model <- glmnet(x=X, y=Y,alpha  = 1,lambda = l.lasso.min)
        as.matrix(lasso.model$beta) #not contain intercept
      })
    } else {
      ct_select=lapply(c(0:(length(lo)-1)), function(id) { #for each location with layer0
        df1=as.data.frame(pseudo.GCT[,c((id*nct+1):(id*nct+nct))])
        cv_5 = trainControl(method = "cv", number = 5)
        y=df1[,1];x=df1[,-1]
        if (quantile(y,0.95)==0){
          y=y+runif(length(y),1e-10,1e-6)
        }
        elnet = train(
          x=x,y=y, data = df1,
          method = "glmnet",
          trControl = cv_5
        )
        bestm=elnet$finalModel
        coeff=coef(bestm,s=elnet$bestTune$lambda)
        as.matrix(coeff)[-1,] #contain intercept
      })
    }
    if(all(sapply(lapply(ct_select,FUN=rownames), FUN = identical, rownames(ct_select[[1]])))) {
      cts=rlist::list.cbind(ct_select)
      colnames(cts)=lo
      cts
    } else {ct_select}
  })
}

################################################################################

#' SPADE:constrained nonlinear optimization to estimate cell type proportion
#'
#' @param bulk  spatial gene expression
#' @param reference cell type reference
#' @param w weights for genes, default is 1
#' @param yNorm normalization method for spatial data, default is cpm
#' @param bNorm normalization method for reference, default is cpm
#'
#' @return
#' @export
#'
#' @examples
spadeDecon <- function(bulk,reference,w=1,yNorm='cpm',bNorm='cpm'){
  bulk1=normcount(bulk,yNorm)
  reference1=normcount(reference,bNorm)
  nct <- ncol(reference1)
  func <- function(P){
    return(sum(w*abs(reference1%*%P  - y), na.rm = TRUE))
  }
  Pinit <- rep(0,nct)
  estp=matrix(0,nrow=ncol(bulk1),ncol=ncol(reference1))
  rownames(estp)=colnames(bulk1);colnames(estp)=colnames(reference1)
  for (i in 1:ncol(bulk1)){
    y=bulk1[,i]
    hin <- function(P){
      P_mat <- matrix(P, nrow = nct)
      c(P)}
    heq <- function(P){ # equality constraint
      P_mat <- matrix(P, nrow = nct)
      c(1-colSums(P_mat))}
    aug_res <- alabama::auglag(Pinit, func, hin = hin,heq=heq,control.outer = list(NMinit=F,trace = F,itmax = 100))
    estp[i,]=aug_res$par
  }
  estp2=abs(estp)
  estp3=t(apply(estp2,1,function(x){x/sum(x)}))
  list(estp3)
}

########################################################################
#' SPADE deconvolution for each spatial domain
#'
#' @param stcount spatial gene expression
#' @param scref cell type reference from single cell RNA-seq data
#' @param sign_list marker genes lists for each cell type
#' @param loc data.frame contains location information
#' @param ctData estimated coefficients for each cell type within each domain
#' @param SVG differentially expressed genes between domains
#' @param cutoff the minimum percent of nonzero entries, default value is 0
#' @param offset the value add on the threshold to binarize coefficient from lasso regression, if not provide by user, then use the mean of all coefficients
#' @param yNorm the normalization method for spatial gene expression, default is cpm
#' @param bNorm the normalization method for reference, default is cpm.
#'
#' @return a matrix for cell type proportion
#' @export
#'

SPADE=function(stcount,scref,sign_list,loc,ctData,SVG,cutoff=0,offset,yNorm='cpm',bNorm='cpm'){
  CTest=lapply(c(0:nlay),function(ilay){
    locs=loc[loc[,'domain']==ilay,][,'location']
    marker1=unlist(sign_list[names(sign_list) %in% ctcut1(ctData,c(ilay+1),cutoff=0,offset=offset)]) #colnames(truepro)
    genelist=Reduce(intersect,list(rownames(scref),rownames(stcount),unique(marker1)))
    ref=scref[genelist,ctcut1(ctData,c(ilay+1),cutoff=0,offset=offset)]#ctcut2(i,c(ilay+1))]
    stcount3=stcount[genelist,locs] #filter location
    # run deconv
    est1=spadeDecon(stcount3,ref,w=1,yNorm,bNorm)
    est1[[1]]
  })
  return(CTest)
}

