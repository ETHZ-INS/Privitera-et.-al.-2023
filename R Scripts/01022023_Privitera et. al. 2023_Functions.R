

dosvacor <- function(SE, form=NULL, form0=~1, ...){
  CD <- as.data.frame(colData(SE))
  mm <- model.matrix(form, data=CD)
  mm0 <- model.matrix(form0, data=CD)
  dds <- DESeqDataSetFromMatrix(round(assay(SE)), as.data.frame(colData(SE)), form)
  dds <- estimateSizeFactors(dds)
  en <- as.matrix(assay(vst(dds, blind=FALSE)))
  sv <- sva(en, mm, mm0, n.sv=NULL, ...)
  n.sv <- sv$n.sv
  sv <- sv$sv
  
  colnames(sv) <- paste0("SV",1:ncol(sv))
  X <- cbind(mm, sv)
  mm2 <- cbind(mm[,1,drop=F],sv,mm[,-1,drop=F])
  H <- solve(t(X)%*%X)%*%t(X)
  b <- (H%*%t(en))
  cn <- setdiff(colnames(X),setdiff(colnames(mm), colnames(mm0)))  
  cn <- setdiff(cn, "(Intercept)")
  encor <- en - t(as.matrix(X[,cn]) %*% b[cn,])
  SE <- SE[row.names(encor),]
  colData(SE) <- cbind(colData(SE), sv)
  assays(SE)$corrected <- encor
  return(SE)
}


#Base function for Boxplots----------------------------------------------------

BoxplotMP <- function(Plotdata, x, y, z) {
  
  ggplot(Plotdata,aes({{x}}, {{y}},color={{z}}))+
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge())+ 
    scale_color_manual(values=c("grey18", "steelblue")) +
    theme_light() + 
    theme(plot.title = element_text(hjust = 0), 
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid = element_line(linetype="dashed"), 
          strip.background = element_rect(fill = "grey18"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") 
     
}


#Base function for Volcanoplot-------------------------------------------------


VolcanoPlotsEdgeRMP <- function(res, FDRcutoff = 0.05, FCcutoff = 0, main = "Volcanoplots", limits = NULL, lab = NULL,PvalLimit = 0){
  library("scales")
  library("ggplot2")
  library("ggrepel")
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  
  PlotData <- NULL
  for(i in names(res)){
    PlotData <- rbind(PlotData, data.frame(res[[i]]$table[,c("logFC","PValue")], FDR = p.adjust(res[[i]]$table$PValue, method = "BH"), Test = i, Gene = rownames(res[[i]])))
  }
  
  change <- NULL
  for(i in 1:dim(PlotData)[1]){
    if(PlotData[i,"logFC"] > FCcutoff & PlotData[i,"FDR"] <= FDRcutoff){
      change<-append(change,"up-regulated")
    }
    else if(PlotData[i,"logFC"] < -FCcutoff & PlotData[i,"FDR"] <= FDRcutoff){
      change<-append(change,"down-regulated")
    }
    else{
      change<-append(change,"not significant")
    }
  }
  
  PlotData$lab <- rep(NA,nrow(PlotData))
  PlotData$change <- factor(change, levels = c("not significant","up-regulated","down-regulated"))
  PlotData$Test <- factor(PlotData$Test, levels = (names(res)))
  if(!is.null(lab)){
    PlotData[PlotData$Gene %in% lab,"lab"] <- PlotData[PlotData$Gene %in% lab,"Gene"]
  }
  
  p1 <- ggplot(PlotData,aes(logFC,PValue + PvalLimit,colour=change, label = lab))+
    geom_point()+
    guides(alpha = "none", size = "none") +
    #Change axis range
    scale_x_continuous(limits = limits) +
    scale_y_continuous(trans=reverselog_trans(10)) + #, limits = c(1e+00, 1e-10) can be added when necessary
    
    #Change axis titles
    labs(x = "Fold Change", y = "P-Value")+
    
    #Change colours/size/transparency
    scale_colour_manual(values = c("not significant" = "grey18","up-regulated" = "firebrick","down-regulated" = "royalblue")) + # Modify point colour
    scale_size_manual(values = c(1, 1, 1)) + # Modify point size
    scale_alpha_manual(values = c(1, 0.5, 1)) + # Modify point transparency
    
    #Change theme and theme/legend elements
    theme_light()+
    theme(panel.grid = element_line(linetype="dashed"), strip.background = element_rect(fill = "grey18"),legend.title = element_text(colour = "grey18", size = 10, face = "bold")) +
    ggtitle(main) +
    geom_text_repel(color = "black", size =3)+
    #Grid
    facet_grid(.~Test)
  
  if(!is.null(limits)){
    p1 <- p1 + scale_x_continuous(limits = limits)
  }
  
  print(p1)
}


```