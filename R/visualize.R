
#' Dotplot to compare estimated proportion with true proportion
#'
#' @param trueprop true proportion
#' @param estprop estimated proportion
#' @param nct number of cell types
#' @param colors defind colors for each cell type
#'
#' @return dot plot
#' @export
#'
dot_plot=function(trueprop,estprop,nct,colors){
  trueprop=as.data.frame(trueprop)
  trueprop$sample=rownames(trueprop)
  truep2=tidyr::gather(trueprop,celltype,proportion,c(1:nct))
  dat=as.data.frame(estprop)
  dat$sample=rownames(dat)
  dat=gather(dat,celltype,proportion,c(1:nct))
  dat1=cbind(truep2,est=dat$proportion)
  ggplot(dat1,aes(x=proportion,y=est,color=celltype)) +
    geom_point()+geom_abline(size=0.8)+xlab("True Proportion")+
    ylab("Estimated Proportion\n")+theme_classic()+
    theme(axis.text.y=element_blank(),
          axis.text=element_text(face='bold',size=10),
          axis.title=element_text(face='bold',size=12),
          strip.text=element_text(face='bold',size=12),
          legend.title = element_blank(),
          axis.ticks.y = element_blank(),
          legend.text = element_text(size=12,face='bold'),
          panel.spacing.x = unit(1, "lines"))+
    scale_color_manual(values=colors)
}

#' Scatter pie plot for cell type composition at each location
#'
#' @param fuldat data.frame contains cell type proportion and spatial location
#' @param dat data.frame for cell type proportion
#' @param fillcolors identify colors for each cell type
#'
#' @return scatter pie plot
#' @export
#'

scatter_pie=function(fuldat,dat,fillcolors){
  ggplot() + geom_scatterpie(aes(x=x_pixel, y=y_pixel), data=fuldat, cols=colnames(dat),pie_scale=0.5,color=NA) + coord_equal()+
    scale_fill_manual(values=fillcolors)+
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.title =element_blank(),
          strip.text=element_text(face='bold',size=12),
          legend.position='right',
          legend.title = element_blank(),
          legend.text = element_text(size=12,face='bold'))
    }


#' scatter plot for cell types at each location
#'
#' @param propdat cell type proportion
#' @param locdat spatial location
#' @param celltype what cell types to plot
#'
#' @return scatterplot shows the cell type proportion at each location
#' @export
#'

featureplot=function(propdat,locdat,celltype){
  ctprop = as.data.frame(propdat[,celltype])
  rownames(ctprop)=rownames(propdat)
  colnames(ctprop)=celltype
  ctprop2 = as.data.frame(apply(ctprop,2,function(x){
    (x - min(x)) / (max(x) - min(x))
  } ))
  ctprop2$x = as.numeric(locdat$x_pixel)
  ctprop2$y = as.numeric(locdat$y_pixel)
  ctprop3 = tidyr::gather(ctprop2,cellType,values,1:length(celltype))
  ggplot(ctprop3, aes(x, y)) +
    geom_point(aes(colour = values),size = 2) +
    scale_color_viridis_c(option = "A")+
    facet_wrap(~cellType)+
    theme(panel.background = element_blank(),
          plot.background = element_blank(),
          axis.text =element_blank(),
          axis.ticks =element_blank(),
          axis.title =element_blank(),
          strip.text = element_text(size = 12,face="bold"),
          legend.position = "bottom",
          legend.title =  element_blank())
}
