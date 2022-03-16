MAEdgeR <- function(res, FDRcutoff = 0.05, main = "MAPlot"){
  res <- res$table
  res$FDR <- p.adjust(res$PValue, method = "BH")
  res$significant <- res$FDR < FDRcutoff
  res$significant <- factor(as.character(res$significant), levels = c("FALSE","TRUE"))
  res <- res[order(-res$FDR),]
  p1 <- ggplot(res,aes(logCPM,logFC, color = significant)) + 
    geom_point(size = 0.5) + 
    scale_color_manual(values = c("black","red")) + 
    theme_bw() +
    ggtitle(main)
  print(p1)
}


VolcanoPlotsEdgeR <- function(res, FDRcutoff = 0.05, FCcutoff = 0, main = "Volcanoplots", limits = NULL, lab = NULL,PvalLimit = 0){
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
    scale_colour_manual(values = c("black","red","blue")) + 
    scale_y_continuous(trans=reverselog_trans(10)) + 
    geom_text_repel(color = "black", size =3) +
    theme_bw() +
    facet_grid(.~Test) + 
    ggtitle(main) +
    ylab("p.value")
  
  if(!is.null(limits)){
    p1 <- p1 + scale_x_continuous(limits = limits)
  }
  
  print(p1)
}

VolcanoPlotsEdgeR_nopoints <- function(res, FDRcutoff = 0.05, FCcutoff = 0, main = "Volcanoplots", dpi = 300){
  library("scales")
  reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
  }
  
  PlotData <- NULL
  for(i in names(res)){
    PlotData <- rbind(PlotData, data.frame(res[[i]]$table[,c("logFC","PValue")], FDR = p.adjust(res[[i]]$table$PValue, method = "BH"), Test = i))
  }
  
  FCcutoff <- 0
  fdrCutoff <- 0.05
  change <- NULL
  for(i in 1:dim(PlotData)[1]){
    if(PlotData[i,"logFC"] > FCcutoff & PlotData[i,"FDR"] <= fdrCutoff){
      change<-append(change,"up-regulated")
    }
    else if(PlotData[i,"logFC"] < -FCcutoff & PlotData[i,"FDR"] <= fdrCutoff){
      change<-append(change,"down-regulated")
    }
    else{
      change<-append(change,"not significant")
    }
  }
  PlotData$change <- factor(change, levels = c("not significant","up-regulated","down-regulated"))
  
  p1 <- ggplot(PlotData,aes(logFC,PValue,colour=change))+
    geom_point(shape = NA)+
    scale_colour_manual(values = c("black","red","blue")) + 
    scale_y_continuous(trans=reverselog_trans(10)) + 
    theme_bw() +
    facet_grid(.~Test) + 
    ggtitle(main)
  print(p1)
}

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


CreateGoPlotData <- function(Results, max_n = 1000, nodeSize = 10, adjp = 0.05, fisher.p = 0.05, categories = c("CC","BP","MF")){
  sig <- topTags(Results,n = max_n, p.value = adjp)$table
  sg <- rownames(sig)
  
  gene_universe <- rownames(Results$table) %in% sg
  gene_universe <- as.numeric(gene_universe)
  gene_universe <- factor(gene_universe)
  names(gene_universe) <- rownames(Results$table)
  
  df <- NULL
  for(i in categories){
    go_data <- new("topGOdata",
                   ontology = i,
                   allGenes = gene_universe,
                   nodeSize = nodeSize,
                   annotationFun = annFUN.org,
                   mapping = "org.Mm.eg",
                   ID = "symbol")
    
    go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
    allGO = usedGO(object = go_data)
    allGOGenes = genesInTerm(go_data)
    go_table <- GenTable(go_data, weightFisher = go_test,
                              orderBy = "weightFisher", ranksOf = "weightFisher",
                              topNodes = length(allGO))
    go_table$weightFisher <- as.numeric(go_table$weightFisher)
    go_table  <- go_table[go_table$weightFisher <= fisher.p,]
    go_genes <- rep(" ",nrow(go_table))
    
    for(j in 1:nrow(go_table)){
      ggs <- allGOGenes[[go_table$GO.ID[j]]][allGOGenes[[go_table$GO.ID[j]]] %in% sg]
      go_genes[j] <- paste(ggs, collapse = ", ")
    }
    go_table$Genes <- go_genes
    go_table$category <- i
    go_table$ID <- go_table$GO.ID
    go_table$adj_pval <- go_table$weightFisher
    sig$ID <- rownames(sig)
    
    circ <- circle_dat(go_table,sig)
    df <- rbind(df,circ)
  }
  #df$adj_pval <- as.numeric(df$adj_pval)
  return(df)
}

CreateGoTableData <- function(genes, universe, nodeSize = 10, fisher.p = 0.05, categories = c("CC","BP","MF")){
  
  
  gene_universe <- universe %in% genes
  gene_universe <- as.numeric(gene_universe)
  gene_universe <- factor(gene_universe)
  names(gene_universe) <- universe
  
  df <- NULL
  for(i in categories){
    go_data <- new("topGOdata",
                   ontology = i,
                   allGenes = gene_universe,
                   nodeSize = nodeSize,
                   annotationFun = annFUN.org,
                   mapping = "org.Mm.eg",
                   ID = "symbol")
    
    go_test <- runTest(go_data, algorithm = "weight01", statistic = "fisher")
    allGO = usedGO(object = go_data)
    allGOGenes = genesInTerm(go_data)
    go_table <- GenTable(go_data, weightFisher = go_test,
                         orderBy = "weightFisher", ranksOf = "weightFisher",
                         topNodes = length(allGO))
    go_table$weightFisher <- as.numeric(go_table$weightFisher)
    go_table  <- go_table[go_table$weightFisher <= fisher.p,]
    go_genes <- rep(" ",nrow(go_table))
    
    for(j in 1:nrow(go_table)){
      ggs <- allGOGenes[[go_table$GO.ID[j]]][allGOGenes[[go_table$GO.ID[j]]] %in% sg]
      go_genes[j] <- paste(ggs, collapse = ", ")
    }
    go_table$category <- i
    go_table$Genes <- go_genes

    df <- rbind(df,go_table)
  }
  #df$adj_pval <- as.numeric(df$adj_pval)
  return(df)
}
