loadCounts <- function(filename) {
  cts_all <- read.table(gzfile(filename), sep='\t', header=T, stringsAsFactors = FALSE)
  rownames(cts_all) <- cts_all$Geneid
  cts_all$Geneid <- NULL
  rowdata <- cts_all[,1:5]
  cts_all <- as.matrix(cts_all[,-c(1:5)])
  colnames(cts_all) <- sapply(strsplit(colnames(cts_all), '\\.'), function(x) tail(x,2)[1])
  cts_all <- cts_all[,order(colnames(cts_all))]
  return(cts_all)
}

getRowdata <- function(filename) {
  cts_all <- read.table(gzfile(filename), sep='\t', header=T, stringsAsFactors = FALSE)
  rownames(cts_all) <- cts_all$Geneid
  cts_all$Geneid <- NULL
  rowdata <- cts_all[,1:5]
  colnames(rowdata)[colnames(rowdata) == 'Length'] <- 'basepairs'
  return(rowdata)
}

# plotClusterMeans <- function(mycpheatctsout, x_time, lwd=.1, col=1, ...) {
#   N <- max(mycpheatctsout$clusterID)
#   par(mfrow=c(N,2))
#   for(k in 1:N) {
#     x <- t(mycpheatctsout$countmeans[mycpheatctsout$clusterID==k,seq_along(x_time)])
#     y <- t(mycpheatctsout$countmeans[mycpheatctsout$clusterID==k,seq_along(x_time)+length(seq_along(x_time))])
#     matplot(x_time, x,lty=1, col=col, lwd=lwd, type='l', 
#             ylim=range(c(x,y)), ...)
#     points(x_time, rowMeans(x), pch='o', cex=1.2, col='red')
#     abline(h=0, lwd=2, lty=2)
#     matplot(x_time, y,lty=1, col=col, lwd=lwd, type='l', 
#             ylim=range(c(x,y)), ...)
#     points(x_time, rowMeans(y), pch='o', cex=1.2, col='red')
#     abline(h=0, lwd=2, lty=2)
#   }  
# }

plotClusterMeans <- function(mycpheatctsout, x_time, genotypes=2, lwd=.1, 
                             plot_means = TRUE,
                             lwdRenorm = FALSE, col=1, ...) {
  N <- max(mycpheatctsout$clusterID)
  if(lwdRenorm) {
    origlwd <- lwd
    lwd <- (1/log10(sapply(mycpheatctsout$clusters_gene_list, length)))^3
    lwd <- lwd/max(lwd)
    lwd <- lwd * origlwd
  } else lwd <- rep(lwd, N)
  par(mfrow=c(N,genotypes))
  for(k in 1:N) {
    X <- mycpheatctsout$countmeans[mycpheatctsout$clusterID==k,]
    for(j in 1:genotypes) {
      y <- t(X[,seq_along(x_time)+(j-1)*length(seq_along(x_time))])
      matplot(x_time, y,lty=1, col=col, lwd=lwd[k], type='l', 
              ylim=range(X), ...)
      if(plot_means) points(x_time, rowMeans(y), pch='o', cex=1.2, col='red')
      abline(h=0, lwd=2, lty=2)
    }
  }  
}

modClProf <- function(x, max.nchar=200) {
  # add the log10(p.adjust) column
  x@compareClusterResult$log10.p.adjust <- log10(x@compareClusterResult$p.adjust)
  # reduce the description length
  descr <- factor(gsub('\n', ' ', x@compareClusterResult$Description))
  descr.levels <- levels(descr)
  subsrt.stop <- nchar(descr.levels)
  substr.start <- sapply(subsrt.stop-max.nchar, function(x) max(c(x,1)))
  new.descr.levels <- substr(descr.levels, substr.start, subsrt.stop)
  if(any(duplicated(new.descr.levels))) {
    idx <- duplicated(new.descr.levels)
    new.descr.levels[idx] <- paste(new.descr.levels[idx], 1:length(which(idx)), sep='_')
  }
  levels(descr) <- new.descr.levels
  new.descr <- as.character(descr)
  x@compareClusterResult$Description <- new.descr
  # return the object
  return(x)
}

### reorder pheatmap clusters based on display order
# clst <- cutree(p$tree_row, k=22)
# clst <- factor(clst, unique(clst[p$tree_row$order]))
# clst <- as.numeric(clst)

# ItemsList <- venn(sc_ec_marker_genes, zcolor = 'style', box = F)
# sc_ec_marker_genes_intersections <- attributes(ItemsList)$intersections

stabilizedFpkms <- function(dds) {
  stopifnot(require(DESeq2))
  stopifnot('basepairs' %in% colnames(rowData(dds)))
  vsd <- vst(dds, blind=FALSE)
  feature_length_norm <- rowData(dds)$basepairs/1000
  libsize_norm <- mean(colSums(counts(dds)))/10^6
  stabilized_fpkms <- 2^assay(vsd)/feature_length_norm/libsize_norm
  return(stabilized_fpkms)
}

getClusterGenes <- function(pheat_out, n) {
  ec_cutree_row <- cutree(pheat_out$tree_row, n)
  clusters_gene <- split(names(ec_cutree_row), ec_cutree_row)
  clusters_gene <- clusters_gene[unique(ec_cutree_row[pheat_out$tree_row$order])]
  names(clusters_gene) <- 1:length(clusters_gene)
  return(clusters_gene)
}

non_redundant <- function(filename) {
  file_lines <- readLines(filename)
  gene_list <- unlist(strsplit(toupper(file_lines), ',|;'))
  gene_table <- table(gene_list)
  out <- lapply(length(file_lines):1, function(i) {
    x <- names(gene_table)[gene_table==i]
    if(length(x)==0) return('-') else
      return(paste(x, collapse=';'))
    })
  out <- unlist(out)
  cat(out, sep='\n')
}

log2log2Plot <- function(de1, de2, diffres = NULL,
                         l2fc_thr=c(0.6, Inf), 
                         padj_thr=1e-1,
                         replacePadjNAby = 1,
                         replaceLog2fcNAby = 0,
                         expr_filter=NULL,
                         gene_labels = '',
                         max.overlaps = 10,
                         seed = 1,
                         global_trend = FALSE,
                         yellow_trend = FALSE,
                         blue_trend = FALSE
                         #color_scale = c('grey70','blue','gold','forestgreen','red')
) {
  
  color_scale <- c('grey70','blue','gold','forestgreen','red')
  
  # cg <- intersect(rownames(de1), rownames(de2))
  # de1 <- data.frame(de1[cg,])
  # de2 <- data.frame(de2[cg,])
  
  l2fcmat <- mergeTwoRes(de1, de2, cols = 'log2FoldChange', 
                         suffixes = c('.x','.y'), replaceNAby = replaceLog2fcNAby)
  if(!is.null(diffres)) {
    l2fcmat <- mergeTwoRes(l2fcmat, diffres, 'log2FoldChange', replaceLog2fcNAby,
                           suffixes = '.diff')
  }
  padjmat <- mergeTwoRes(de1, de2, cols = 'padj', 
                         suffixes = c('.x','.y'), replaceNAby = replacePadjNAby)
  if(!is.null(diffres)) {
    padjmat <- mergeTwoRes(padjmat, diffres, 'padj', replacePadjNAby,
                           suffixes = '.diff')
  }
  de1 <- cbind(l2fcmat[,1,drop=F],padjmat[,1,drop=F])
  colnames(de1) <- sub('\\..*','',colnames(de1))
  de2 <- cbind(l2fcmat[,2,drop=F],padjmat[,2,drop=F])
  colnames(de2) <- sub('\\..*','',colnames(de2))
  if(!is.null(diffres)) {
    de3 <- cbind(l2fcmat[,3,drop=F],padjmat[,3,drop=F])
    colnames(de3) <- sub('\\..*','',colnames(de3))
  }
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  
  diffexpressed <- gene_class(de1, de2, padj_thr_min, l2fc_thr_min)$class
  if(!is.null(expr_filter)) diffexpressed[!expr_filter] <- 0
  if(!is.null(diffres)) {
    nondiff <- de3$padj>padj_thr_min | abs(de3$log2FoldChange)<l2fc_thr_min
    diffexpressed[nondiff] <- 0
  }

  deSat <- function(de) {
    satFCup <- de$log2FoldChange > l2fc_thr_max
    de$log2FoldChange[satFCup] <- l2fc_thr_max
    satFCdn <- de$log2FoldChange < -l2fc_thr_max
    de$log2FoldChange[satFCdn] <- -l2fc_thr_max
    de$saturated <- satFCup | satFCdn 
    de$saturated[is.na(de$saturated)] <- FALSE
    return(de)
  }
  
  de1 <- deSat(de1); colnames(de1) <- paste(colnames(de1),'1',sep='_')
  de2 <- deSat(de2); colnames(de2) <- paste(colnames(de2),'2',sep='_')
  if(!is.null(diffres)) {
    de3 <- deSat(de3); colnames(de3) <- paste(colnames(de3),'3',sep='_')
    de <- cbind(de1, de2, de3, diffexpressed=diffexpressed)
  } else {
    de <- cbind(de1, de2, diffexpressed=diffexpressed)
  }
  de$saturated <- de$saturated_1 | de$saturated_2

  de$delabel <- NA
  if(is.null(gene_labels)) {
    de$delabel[de$diffexpressed != 0] <- rownames(de)[de$diffexpressed != 0]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  n_gene_class <- table(de$diffexpressed)
  
  # trend lines
  if( length(which(diffexpressed == 1)) < 2 ) blue_trend <- FALSE
  if( length(which(diffexpressed == 2)) < 2 ) yellow_trend <- FALSE
  if( global_trend ) {
    # calculate lm based on x-sum-of-squares 
    global_coef <- coef(with(de, lm(log2FoldChange_2~log2FoldChange_1)))
  }
  if( blue_trend ) {
    # calculate lm based on x-sum-of-squares 
    blue_coef <- coef(with(subset(de, diffexpressed == 1),
                           lm(log2FoldChange_2~log2FoldChange_1)))
  }
  if( yellow_trend ) {
    # calculate lm based on y-sum-of-squares (and then invert the slope)
    yellow_coef <- coef(with(subset(de, diffexpressed == 2),
                             lm(log2FoldChange_1~log2FoldChange_2)))
  }
  assign('log2log2de', de, envir = globalenv())
  # plot adding up all layers we have seen so far
  g <- ggplot(data=de, aes(x=log2FoldChange_1, y=log2FoldChange_2, 
                      col=diffexpressed, label=delabel, shape=saturated)) +
    geom_point(size=1.5) + 
    # geom_bin2d() + 
    geom_point(data = subset(de, diffexpressed != 0), 
               aes(x=log2FoldChange_1, y=log2FoldChange_2, 
                   col=diffexpressed, shape=saturated))

  if( global_trend ) g <- g +
    geom_abline(intercept = global_coef[1],
                slope = global_coef[2],
                color = 'black') #'dodgerblue4')
  if( blue_trend ) g <- g +
    geom_abline(intercept = blue_coef[1],
                slope = blue_coef[2],
                color = 'royalblue4') #'dodgerblue4')
  if( yellow_trend ) g <- g +
    geom_abline(intercept = yellow_coef[1],
                slope = 1/yellow_coef[2],
                color = 'darkgoldenrod')
  
  g <- g + theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1.5)) +
    geom_text_repel(max.overlaps=max.overlaps, size = 3, colour = "black", seed=seed) +
    scale_x_continuous(expand = expansion(mult = .015)) +
    scale_y_continuous(expand = expansion(mult = c(.01, .015))) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    geom_hline(yintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("no-degs (",n_gene_class[1],")"), 
                           paste0("x-degs (",n_gene_class[2],")"), 
                           paste0("y-degs (",n_gene_class[3],")"),
                           paste0("same (",n_gene_class[4],")"),
                           paste0("opposite (",n_gene_class[5],")")
                         )[n_gene_class>0], 
                         type=color_scale[n_gene_class>0])
  # if(is.finite(l2fc_thr_max)) {
  #   g <- g + xlim(-l2fc_thr_max, l2fc_thr_max) + ylim(-l2fc_thr_max, l2fc_thr_max)
  # }
  return(g)
  
}

# with this function you can set l2fc_thr specific for x and y
log2log2Plot2 <- function(de1, de2, diffres = NULL,
                         l2fc_thr=c(0.6, Inf), 
                         l2fc_thr_ymax=NULL, 
                         padj_thr=1e-1,
                         replacePadjNAby = 1,
                         replaceLog2fcNAby = 0,
                         expr_filter=NULL,
                         gene_labels = '',
                         max.overlaps = 10,
                         seed = 1,
                         global_cor = FALSE,
                         yellow_trend = FALSE,
                         blue_trend = FALSE
                         #color_scale = c('grey70','blue','gold','forestgreen','red')
) {
  
  color_scale <- c('grey70','blue','gold','forestgreen','red')

  # cg <- intersect(rownames(de1), rownames(de2))
  # de1 <- data.frame(de1[cg,])
  # de2 <- data.frame(de2[cg,])
  
  l2fcmat <- mergeTwoRes(de1, de2, cols = 'log2FoldChange', 
                         suffixes = c('.x','.y'), replaceNAby = replaceLog2fcNAby)
  if(!is.null(diffres)) {
    l2fcmat <- mergeTwoRes(l2fcmat, diffres, 'log2FoldChange', replaceLog2fcNAby,
                           suffixes = '.diff')
  }
  padjmat <- mergeTwoRes(de1, de2, cols = 'padj', 
                         suffixes = c('.x','.y'), replaceNAby = replacePadjNAby)
  if(!is.null(diffres)) {
    padjmat <- mergeTwoRes(padjmat, diffres, 'padj', replacePadjNAby,
                           suffixes = '.diff')
  }
  de1 <- cbind(l2fcmat[,1,drop=F],padjmat[,1,drop=F])
  colnames(de1) <- sub('\\..*','',colnames(de1))
  de2 <- cbind(l2fcmat[,2,drop=F],padjmat[,2,drop=F])
  colnames(de2) <- sub('\\..*','',colnames(de2))
  if(!is.null(diffres)) {
    de3 <- cbind(l2fcmat[,3,drop=F],padjmat[,3,drop=F])
    colnames(de3) <- sub('\\..*','',colnames(de3))
  }
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  if(is.null(l2fc_thr_ymax)) 
    l2fc_thr_ymax <- l2fc_thr_max
  padj_thr_min <- padj_thr[1]
  
  diffexpressed <- gene_class(de1, de2, padj_thr_min, l2fc_thr_min)$class
  if(!is.null(expr_filter)) diffexpressed[!expr_filter] <- 0
  if(!is.null(diffres)) {
    nondiff <- de3$padj>padj_thr_min | abs(de3$log2FoldChange)<l2fc_thr_min
    diffexpressed[nondiff] <- 0
  }
  
  deSat <- function(de, l2fc_thr_max) {
    satFCup <- de$log2FoldChange > l2fc_thr_max
    de$log2FoldChange[satFCup] <- l2fc_thr_max
    satFCdn <- de$log2FoldChange < -l2fc_thr_max
    de$log2FoldChange[satFCdn] <- -l2fc_thr_max
    de$saturated <- satFCup | satFCdn 
    de$saturated[is.na(de$saturated)] <- FALSE
    return(de)
  }
  
  de1 <- deSat(de1, l2fc_thr_max); colnames(de1) <- paste(colnames(de1),'1',sep='_')
  de2 <- deSat(de2, l2fc_thr_ymax); colnames(de2) <- paste(colnames(de2),'2',sep='_')
  if(!is.null(diffres)) {
    de3 <- deSat(de3); colnames(de3) <- paste(colnames(de3),'3',sep='_')
    de <- cbind(de1, de2, de3, diffexpressed=diffexpressed)
  } else {
    de <- cbind(de1, de2, diffexpressed=diffexpressed)
  }
  de$saturated <- de$saturated_1 | de$saturated_2
  
  de$delabel <- NA
  if(is.null(gene_labels)) {
    de$delabel[de$diffexpressed != 0] <- rownames(de)[de$diffexpressed != 0]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  n_gene_class <- table(de$diffexpressed)
  
  # trend lines
  if( length(which(diffexpressed == 1)) < 2 ) blue_trend <- FALSE
  if( length(which(diffexpressed == 2)) < 2 ) yellow_trend <- FALSE
  if( blue_trend ) {
    # calculate lm based on x-sum-of-squares 
    blue_coef <- coef(with(subset(de, diffexpressed == 1),
                           lm(log2FoldChange_2~log2FoldChange_1)))
  }
  if( yellow_trend ) {
    # calculate lm based on y-sum-of-squares (and then invert the slope)
    yellow_coef <- coef(with(subset(de, diffexpressed == 2),
                             lm(log2FoldChange_1~log2FoldChange_2)))
  }
  assign('log2log2de', de, envir = globalenv())
  # plot adding up all layers we have seen so far
  g <- ggplot(data=de, aes(x=log2FoldChange_1, y=log2FoldChange_2, 
                           col=diffexpressed, label=delabel, shape=saturated)) +
    geom_point(size=1.5) + 
    # geom_bin2d() + 
    geom_point(data = subset(de, diffexpressed != 0), 
               aes(x=log2FoldChange_1, y=log2FoldChange_2, 
                   col=diffexpressed, shape=saturated))
  
  if( blue_trend ) g <- g +
    geom_abline(intercept = blue_coef[1],
                slope = blue_coef[2],
                color = 'royalblue4') #'dodgerblue4')
  if( yellow_trend ) g <- g +
    geom_abline(intercept = yellow_coef[1],
                slope = 1/yellow_coef[2],
                color = 'darkgoldenrod')
  
  g <- g + theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1.5)) +
    geom_text_repel(max.overlaps=max.overlaps, size = 3, colour = "black", seed=seed) +
    scale_x_continuous(expand = expansion(mult = .015)) +
    scale_y_continuous(expand = expansion(mult = c(.01, .015))) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    geom_hline(yintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("no-degs (",n_gene_class[1],")"), 
                           paste0("x-degs (",n_gene_class[2],")"), 
                           paste0("y-degs (",n_gene_class[3],")"),
                           paste0("same (",n_gene_class[4],")"),
                           paste0("opposite (",n_gene_class[5],")")
                         )[n_gene_class>0], 
                         type=color_scale[n_gene_class>0])

  # if( global_cor ) {
  #   annotations <- data.frame(
  #     xpos = c(Inf),
  #     ypos =  c(Inf),
  #     annotateText = as.character(round(cor(de$log2FoldChange_1, de$log2FoldChange_2, method='s'),2)),
  #     hjustvar = c(1) ,
  #     vjustvar = c(1))
  #   g <- g +
  #     geom_text(data=annotations,aes(x=xpos,y=ypos,hjust=hjustvar,vjust=vjustvar,label=annotateText))
  # }
  # 
  # if(is.finite(l2fc_thr_max)) {
  #   g <- g + xlim(-l2fc_thr_max, l2fc_thr_max) + ylim(-l2fc_thr_max, l2fc_thr_max)
  # }
  return(g)
  
}


gene_class <- function(resLFC1, resLFC2, alpha = .1, l2fc_thr_min=0) {
  merged_df <- merge(data.frame(resLFC1[,c('log2FoldChange','padj')]),
                     data.frame(resLFC2[,c('log2FoldChange','padj')]), by = 0, all=TRUE)
  resLFC1 <- merged_df[,2:3]; rownames(resLFC1) <- merged_df[,1]; 
  colnames(resLFC1) <- sub('\\..*','',colnames(resLFC1))
  resLFC2 <- merged_df[,4:5]; rownames(resLFC2) <- merged_df[,1]; 
  colnames(resLFC2) <- sub('\\..*','',colnames(resLFC2))
  resLFC1$padj[is.na(resLFC1$padj)] <- 1; resLFC1$log2FoldChange[is.na(resLFC1$log2FoldChange)] <- 0
  resLFC2$padj[is.na(resLFC2$padj)] <- 1; resLFC2$log2FoldChange[is.na(resLFC2$log2FoldChange)] <- 0
  
  signif1 <- resLFC1$padj < alpha & abs(resLFC1$log2FoldChange) > l2fc_thr_min
  signif2 <- resLFC2$padj < alpha & abs(resLFC2$log2FoldChange) > l2fc_thr_min
  logFCsign1 <- resLFC1$log2FoldChange > 0
  logFCsign2 <- resLFC2$log2FoldChange > 0
  
  sign1pos <- signif1 & !signif2 & logFCsign1 # 1
  sign2pos <- !signif1 & signif2 & logFCsign2 # 2
  sign1neg <- signif1 & !signif2 & !logFCsign1 # 3
  sign2neg <- !signif1 & signif2 & !logFCsign2 # 4
  pospos   <- signif1 & signif2 & logFCsign1 & logFCsign2 # 5
  posneg   <- signif1 & signif2 & logFCsign1 & !logFCsign2 # 6
  negneg   <- signif1 & signif2 & !logFCsign1 & !logFCsign2 # 7
  negpos   <- signif1 & signif2 & !logFCsign1 & logFCsign2 # 8
  
  gene_class <- rep(0,nrow(resLFC1)) # both no degs
  gene_class[sign1pos|sign1neg] <- 1 # x-degs
  gene_class[sign2pos|sign2neg] <- 2 # y-degs
  gene_class[pospos|negneg] <- 3 # concordant
  gene_class[posneg|negpos] <- 4 # discordant
  
  gene_class <- factor(gene_class, levels=0:4)
  return(data.frame(class=gene_class, row.names = rownames(resLFC1)))
  
}

# gene_class <- function(resLFC1, resLFC2, alpha = .1, l2fc_thr_min=0) {
#   resLFC1$padj[is.na(resLFC1$padj)] <- 1
#   resLFC2$padj[is.na(resLFC2$padj)] <- 1
#   signif_genes <- 1*(resLFC1$padj < alpha) + 1*(resLFC2$padj < alpha)
#   signif1 <- resLFC1$padj < alpha & abs(resLFC1$log2FoldChange) > l2fc_thr_min
#   signif2 <- resLFC2$padj < alpha & abs(resLFC2$log2FoldChange) > l2fc_thr_min
#   logFCsign1 <- resLFC1$log2FoldChange > 0
#   logFCsign2 <- resLFC2$log2FoldChange > 0
#   
#   sign1pos <- signif1 & !signif2 & logFCsign1 # 1
#   sign2pos <- !signif1 & signif2 & logFCsign2 # 2
#   sign1neg <- signif1 & !signif2 & !logFCsign1 # 3
#   sign2neg <- !signif1 & signif2 & !logFCsign2 # 4
#   pospos   <- signif1 & signif2 & logFCsign1 & logFCsign2 # 5
#   posneg   <- signif1 & signif2 & logFCsign1 & !logFCsign2 # 6
#   negneg   <- signif1 & signif2 & !logFCsign1 & !logFCsign2 # 7
#   negpos   <- signif1 & signif2 & !logFCsign1 & logFCsign2 # 8
#   
#   gene_class <- rep(0,nrow(resLFC1))
#   gene_class[sign1pos|sign1neg] <- 1
#   gene_class[sign2pos|sign2neg] <- 2
#   gene_class[pospos|negneg] <- 3
#   gene_class[posneg|negpos] <- 4
# 
#   gene_class <- factor(gene_class, levels=0:4)
#   return(data.frame(class=gene_class, row.names = rownames(resLFC1)))
#   
#   # palette <- c('black','blue','gold','forestgreen','red')
#   # cols <- palette[1+gene_class]
#   # 
#   # return(data.frame(class=gene_class, col=cols, row.names = rownames(resLFC1)))
#   
# }

# comparisons <- c(
#   "condition_SOD_45_vs_WT_45",
#   "condition_SOD_90_vs_WT_90",
#   "condition_SOD_120_vs_WT_120",
#   "condition_SOD_135_vs_WT_135"
# )
# genotypic_res <- lapply(comparisons, function(name) makeComparison(name, dds_endo_TRAP))

classicFpkm <- function(dds) {
  stopifnot('libsize' %in% names(colData(dds)))
  stopifnot('basepairs' %in% names(rowData(dds)))
  t(t(counts(dds)/rowData(dds)$basepairs)/dds$libsize)*10^9
}

makeComparison <- function(name, obj, shrink=FALSE) {
  ref <- strsplit(name, split = '_vs_')[[1]][2]
  obj$condition <- relevel(obj$condition, ref = ref)
  obj <- DESeq(obj)
  if( shrink ) lfcShrink(obj, coef = name) else results(obj, name = name)
}

makeStrictComparisonOld <- function(name, obj,
                                 counts_filter = -1,
                                 fpkms_filter = -1,
                                 mean_counts_filter = -1,
                                 mean_fpkms_filter = -1,
                                 f = any,
                                 robustSF = FALSE,
                                 shrink=FALSE) {
  fac <- strsplit(name, split = '_')[[1]][1]
  cnd <- strsplit(sub(paste0(fac, '_'), '', name), split = '_vs_')[[1]][1]
  ref <- strsplit(name, split = '_vs_')[[1]][2]
  obj <- obj[,obj[[fac]] %in% c(cnd, ref)]
  obj[[fac]] <- factor(as.character(obj[[fac]]), levels = c(ref, cnd))
  if(robustSF) {
    obj <- estimateSizeFactors(obj)
  } else {
    if(is.null(sizeFactors(obj))) {
      print('sizeFactors not set, computing now for the given dds (obj)')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  ref_fpkms <- fpkm(obj)[,obj[[fac]] == ref, drop=F]
  cnd_fpkms <- fpkm(obj)[,obj[[fac]] == cnd, drop=F]
  mean_fpkms <- cbind(rowMeans(ref_fpkms), rowMeans(cnd_fpkms))
  ref_counts <- counts(obj)[,obj[[fac]] == ref, drop=F]
  cnd_counts <- counts(obj)[,obj[[fac]] == cnd, drop=F]
  mean_counts <- cbind(rowMeans(ref_counts), rowMeans(cnd_counts))
  idx <-
    (apply(ref_counts>counts_filter,1,f) |
       apply(cnd_counts>counts_filter,1,f)) &
    (apply(ref_fpkms>fpkms_filter,1,f) |
       apply(cnd_fpkms>fpkms_filter,1,f)) &
    apply(mean_fpkms>mean_fpkms_filter,1,f) &
    apply(mean_counts>mean_counts_filter,1,f)
  obj <- obj[idx,]
  # in the robust case, re-estimate sizeFactors after the filtering
  if(robustSF) obj <- estimateSizeFactors(obj)
  obj <- DESeq(obj)
  if( shrink ) lfcShrink(obj, coef = name) else results(obj, name = name)
}

makeStrictComparison <- function(name, obj, 
                                 counts_filter = -1, 
                                 fpkms_filter = -1, 
                                 mean_counts_filter = -1, 
                                 mean_fpkms_filter = -1, 
                                 f = any, 
                                 shrink=FALSE) {
  fac <- strsplit(name, split = '_')[[1]][1]
  cnd <- strsplit(sub(paste0(fac, '_'), '', name), split = '_vs_')[[1]][1]
  ref <- strsplit(name, split = '_vs_')[[1]][2]
  obj <- obj[,obj[[fac]] %in% c(cnd, ref)]
  obj[[fac]] <- factor(as.character(obj[[fac]]), levels = c(ref, cnd))
  ref_fpkms <- classicFpkm(obj)[,obj[[fac]] == ref, drop=F]
  cnd_fpkms <- classicFpkm(obj)[,obj[[fac]] == cnd, drop=F]
  mean_fpkms <- cbind(rowMeans(ref_fpkms), rowMeans(cnd_fpkms))
  ref_counts <- counts(obj)[,obj[[fac]] == ref, drop=F]
  cnd_counts <- counts(obj)[,obj[[fac]] == cnd, drop=F]
  mean_counts <- cbind(rowMeans(ref_counts), rowMeans(cnd_counts))
  idx <- 
    (apply(ref_counts>counts_filter,1,f) | 
       apply(cnd_counts>counts_filter,1,f)) & 
    (apply(ref_fpkms>fpkms_filter,1,f) | 
       apply(cnd_fpkms>fpkms_filter,1,f)) &
    apply(mean_fpkms>mean_fpkms_filter,1,f) & 
    apply(mean_counts>mean_counts_filter,1,f)
  obj <- obj[idx,]
  # estimate sizeFactors before testing filtering
  obj <- estimateSizeFactors(obj)
  obj <- DESeq(obj)
  if( shrink ) lfcShrink(obj, coef = name) else results(obj, name = name)
}

makeStableComparison <- function(name, obj, 
                                 counts_filter = -1, 
                                 fpkms_filter = -1, 
                                 mean_counts_filter = -1, 
                                 mean_fpkms_filter = -1, 
                                 f = any, 
                                 robustSF = FALSE,
                                 shrink=FALSE) {
  if(robustSF) {
    obj <- estimateSizeFactors(obj)
  } else {
    if(is.null(sizeFactors(obj))) {
      print('sizeFactors not set, computing now for the given dds (obj)')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  fac <- strsplit(name, split = '_')[[1]][1]
  cnd <- strsplit(sub(paste0(fac, '_'), '', name), split = '_vs_')[[1]][1]
  ref <- strsplit(name, split = '_vs_')[[1]][2]
  #
  ref_fpkms <- fpkm(obj)[,obj[[fac]] == ref, drop=F]
  cnd_fpkms <- fpkm(obj)[,obj[[fac]] == cnd, drop=F]
  mean_fpkms <- cbind(rowMeans(ref_fpkms), rowMeans(cnd_fpkms))
  ref_counts <- counts(obj)[,obj[[fac]] == ref, drop=F]
  cnd_counts <- counts(obj)[,obj[[fac]] == cnd, drop=F]
  mean_counts <- cbind(rowMeans(ref_counts), rowMeans(cnd_counts))
  #
  idx <- 
    (apply(ref_counts>counts_filter,1,f) | 
       apply(cnd_counts>counts_filter,1,f)) & 
    (apply(ref_fpkms>fpkms_filter,1,f) | 
       apply(cnd_fpkms>fpkms_filter,1,f)) &
    apply(mean_fpkms>mean_fpkms_filter,1,f) & 
    apply(mean_counts>mean_counts_filter,1,f)
  # subset object (after fpkm calculation not to alter them)
  obj <- obj[idx,obj[[fac]] %in% c(cnd, ref)]
  obj[[fac]] <- factor(as.character(obj[[fac]]), levels = c(ref, cnd))
  # in the robust case, re-estimate sizeFactors after the filtering
  if(robustSF) obj <- estimateSizeFactors(obj)
  obj <- DESeq(obj)
  if( shrink ) lfcShrink(obj, coef = name) else results(obj, name = name)
}


library(ggplot2)
library(ggrepel)
plotPCA.san <- function (object, intgroup = "condition", ntop = 500, 
                         labels = NA, max.overlaps = 10,
                         fix.coord = TRUE, returnData = FALSE) 
{
  rv <- rowVars(assay(object))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, 
                                                     length(rv)))]
  pca <- prcomp(t(assay(object)[select, ]))
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop = FALSE])
  group <- if (length(intgroup) > 1) {
    factor(apply(intgroup.df, 1, paste, collapse = " : "))
  }
  else {
    colData(object)[[intgroup]]
  }
  
  if(all(is.na(labels))) {
    labels <- rownames(colData(object))
  } 
  if(is.logical(labels)) {
    if(!labels) labels <- rep('', ncol(object)) else labels <- rownames(colData(object))
  }
  stopifnot(length(labels)==length(rownames(colData(object))))
  
  ## Select the PCAs and percentVar that you like instead of 1 and 2
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], group = group, 
                  intgroup.df, name = labels)
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  } else {
    p <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = "group", label = "name")) + 
      geom_point(size = 3) + 
      xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) + 
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      geom_text_repel(size=3, max.overlaps = max.overlaps) 
    if(fix.coord) p <- p + coord_fixed()
    return(p)
  }
  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

enrichPlot <- function(filename, sheet, order.by='n', FUN='max', 
                       show.labels=TRUE) {
  enrich <- read.xlsx(filename, sheet = sheet)
  enrich$n <- as.numeric(sapply(strsplit(enrich$Overlap,'\\/'),'[',1))
  enrich$ratio <- enrich$n/as.numeric(sapply(strsplit(enrich$Overlap,'\\/'),'[',2))
  colnames(enrich)[colnames(enrich) == "Adjusted.P-value"] <- 'p_adjusted'
  colnames(enrich)[colnames(enrich) == "Summary.Term"] <- 'summary_term'
  colnames(enrich)[grep('selected.genes', colnames(enrich))] <- 'summary_genes'
  enrich$p_adjusted <- as.numeric(enrich$p_adjusted)
  enrich$summary_term <- replace_na(enrich$summary_term)
  enrich$summary_genes <- if(show.labels) replace_na(enrich$summary_genes) else NA
  enrich_top1 <- enrich[,c('summary_term','Term','n','ratio','p_adjusted','summary_genes')]
  enrich_top1 <- do.call('rbind', lapply(split(enrich_top1, enrich_top1$summary_term),
                                         function(x) {
                                           y <- x[,order.by]
                                           ind <- which(y == do.call(FUN, as.list(y)))[1]
                                           x[ind,]
                                         }))
  enrich_top1 <- enrich_top1[order(enrich_top1$n),]
  enrich_top1$seq <- 1:nrow(enrich_top1)
  gg <- ggplot(enrich_top1, 
               aes(x=n, 
                   y=seq, 
                   size = ratio,
                   label = if(show.labels) summary_genes else NULL,
                   color = -log10(p_adjusted))) + geom_point(alpha=0.7) + scale_x_continuous(trans='log2')
  if(show.labels) gg <- gg + geom_text_repel(aes(x=n, y=seq), size=2.5, hjust=-.2, colour='black')
  gg <- gg + scale_y_continuous(breaks=enrich_top1$seq, labels = enrich_top1$summary_term)
  gg <- gg + scale_color_gradient(low="blue", high="red") + scale_size_area()
  gg
}

# splitted_genes <- sapply(split(enrich$Genes, enrich$summary_term), 
#                        function(x) strsplit(x, split=';|,'))
# splitted_genes <- splitted_genes[unique(enrich$summary_term)]
# sorted_Genes <- lapply(splitted_genes, function(x) {
#   N <- length(x)
#   tabled_genes <- factor(table(unlist(x)), levels=1:N)
#   rev(sapply(split(names(tabled_genes), tabled_genes), paste, collapse=';'))
# })
# enrich$sorted_N <- unlist(sapply(sorted_Genes, function(x) length(x):1))
# enrich$sorted_Genes <- unlist(sorted_Genes)
# write.csv(enrich[,c('summary_term','sorted_N','sorted_Genes')], file='')


replace_na <- function(x) {
  for(i in 1:length(x))
    if(is.na(x[i])) x[i] <- x[i-1]
  return(x)
}

volcanoPlot <- function(de, 
                        l2fc_thr=c(0.6, Inf), 
                        padj_thr=c(1e-1, -Inf),
                        expr_filter=NULL,
                        gene_labels = NULL,
                        geom_point_size = 1.5, 
                        max.overlaps = 10,
                        seed = 1,
                        color_scale = c("dodgerblue4", "grey50", "firebrick")
) {
  
  de <- data.frame(de)
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  padj_thr_max <- padj_thr[2]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min] <- "UP"
  de$diffexpressed[de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c("DOWN","NO","UP"))
  if(!is.null(expr_filter)) de$diffexpressed[!expr_filter] <- "NO"
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  satPval <- de$padj < padj_thr_max
  de$padj[satPval] <- padj_thr_max
  de$saturated <- satFCup | satFCdn | satPval
  de$saturated[is.na(de$saturated)] <- FALSE
  
  de$delabel <- NA
  if(is.null(gene_labels)) {
    de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  n_gene_class <- table(de$diffexpressed)
  
  assign('volcanoPlotde', de, envir = globalenv())
  
  if(!is.finite(l2fc_thr_max)) l2fc_thr_max <- max(abs(de$log2FoldChange), na.rm=TRUE)
  if(!is.finite(padj_thr_max)) padj_thr_max <- min(de$padj, na.rm=TRUE)
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), 
                      col=diffexpressed, label=delabel, shape=saturated)) +
    geom_point(size=geom_point_size) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1.5)) +
    geom_text_repel(max.overlaps=max.overlaps, size = 3, colour = "black", seed=seed) +
    scale_x_continuous(expand = expansion(mult = .015), 
                       limits = c(-l2fc_thr_max, l2fc_thr_max)) +
    scale_y_continuous(expand = expansion(mult = c(.01, .015)),
                       limits = c(0, -log10(padj_thr_max))) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thr_min), linetype="dashed") + 
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")")
                         )[n_gene_class>0], 
                         type=color_scale[n_gene_class>0])
}

volcanoPlot2 <- function(de, 
                        l2fc_thr=c(0.6, Inf), 
                        padj_thr=c(1e-1, -Inf),
                        expr_filter=NULL,
                        gene_labels = NULL,
                        geom_point_size = 1.5, 
                        max.overlaps = 10,
                        seed = 1,
                        color_scale = c("dodgerblue4", "grey50", "firebrick","gold")
) {
  
  de <- data.frame(de)
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  padj_thr_max <- padj_thr[2]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min] <- "UP"
  de$diffexpressed[de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c("DOWN","NO","UP","LABEL"))
  if(!is.null(expr_filter)) de$diffexpressed[!expr_filter] <- "NO"
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  satPval <- de$padj < padj_thr_max
  de$padj[satPval] <- padj_thr_max
  de$saturated <- satFCup | satFCdn | satPval
  de$saturated[is.na(de$saturated)] <- FALSE
  de$size <- geom_point_size
  
  de$delabel <- NA
  if(!is.null(gene_labels)) {
    de[gene_labels,'diffexpressed'] <- 'LABEL'
    de[gene_labels,'size'] <- geom_point_size * 2
    de[gene_labels,'delabel'] <- gene_labels
  }
  n_gene_class <- table(de$diffexpressed)
  
  assign('volcanoPlotde', de, envir = globalenv())
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), label=delabel,
                      col=diffexpressed, shape=saturated, size=size)) +
    geom_point() + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1.5)) +
    geom_text_repel(max.overlaps=max.overlaps, size = 4, colour = "black", seed=seed) +
    scale_x_continuous(expand = expansion(mult = .015)) +
    scale_y_continuous(expand = expansion(mult = c(.01, .015))) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thr_min), linetype="dashed") + 
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")"),
                           paste0("Labeled (",n_gene_class[4],")")
                         )[n_gene_class>0], 
                         type=color_scale[n_gene_class>0]) + 
    guides(size = FALSE)
}

volcanoPlot3 <- function(de, 
                         l2fc_thr=c(0.6, Inf), 
                         padj_thr=c(1e-1, -Inf),
                         expr_filter=NULL,
                         gene_labels = NULL,
                         geom_point_size = 1.5, 
                         max.overlaps = 10,
                         seed = 1,
                         color_scale = c("dodgerblue4", "grey50", "firebrick",
                                         "lightblue","orange")
) {
  
  de <- data.frame(de)
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  padj_thr_max <- padj_thr[2]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min] <- "UP"
  de$diffexpressed[de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c("DOWN","NO","UP","DOWN-SOMA","UP-SOMA"))
  if(!is.null(expr_filter)) de$diffexpressed[!expr_filter] <- "NO"
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  satPval <- de$padj < padj_thr_max
  de$padj[satPval] <- padj_thr_max
  de$saturated <- satFCup | satFCdn | satPval
  de$saturated[is.na(de$saturated)] <- FALSE
  # de$size <- geom_point_size
  
  de$delabel <- NA
  if(!is.null(gene_labels)) {
    gene_labels <- intersect(gene_labels, rownames(de)[de$diffexpressed != 'NO'])
    de[gene_labels,'diffexpressed'] <- paste(de[gene_labels,'diffexpressed'],'SOMA',sep='-')
    # de[gene_labels,'size'] <- geom_point_size
    # de[gene_labels,'delabel'] <- gene_labels
  }
  n_gene_class <- table(de$diffexpressed)
  
  assign('volcanoPlotde', de, envir = globalenv())
  
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=-log10(padj), label=delabel,
                      col=diffexpressed, shape=saturated)) +
    geom_point(size=geom_point_size) + 
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_rect(colour = "black", size=1.5)) +
    geom_text_repel(max.overlaps=max.overlaps, size = 4, colour = "black", seed=seed) +
    scale_x_continuous(expand = expansion(mult = .015)) +
    scale_y_continuous(expand = expansion(mult = c(.01, .015))) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), linetype="dashed") +
    geom_hline(yintercept=-log10(padj_thr_min), linetype="dashed") + 
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")"),
                           paste0("Depleted - Soma DEGs (",n_gene_class[4],")"), 
                           paste0("Enriched - Soma DEGs (",n_gene_class[5],")")
                         )[n_gene_class>0], 
                         type=color_scale[n_gene_class>0])
}


mergeTwoRes <- function(res1, res2, cols, replaceNAby=NA, suffixes = c(".x",".y")) {
  l2fcmat <- merge(data.frame(res1[,unlist(sapply(cols, grep, colnames(res1))),drop=F]),
                   data.frame(res2[,unlist(sapply(cols, grep, colnames(res2))),drop=F]),
                   by = 0, all=TRUE, suffixes = suffixes)
  rownames(l2fcmat) <- l2fcmat$Row.names
  l2fcmat$Row.names <- NULL
  l2fcmat[is.na(l2fcmat)] <- replaceNAby
  return(l2fcmat)
}


degsMatOld <- function(genotypic_res, padj_thr_min=.1, log2fc_thr_min = 0) {
  l2fcmat <- sapply(genotypic_res, function(x) x$log2FoldChange)
  padjmat <- sapply(genotypic_res, function(x) x$padj)
  degsmat <- 1*(padjmat<padj_thr_min & l2fcmat>log2fc_thr_min) - 
    1*(padjmat<padj_thr_min & l2fcmat<(-log2fc_thr_min))
  degsmat[is.na(degsmat)] <- 0
  rownames(degsmat) <- rownames(genotypic_res[[1]])
  return(degsmat)
}

reportMat <- function(genotypic_res, l2fcNA=NA, padjNA=NA) {

  l2fcmat <- mergeTwoRes(genotypic_res[[1]], genotypic_res[[2]], 
                         'log2FoldChange', l2fcNA, suffixes = c('.1','.2'))
  if(length(genotypic_res)>2)
    for(i in 3:length(genotypic_res))
      l2fcmat <- mergeTwoRes(l2fcmat, genotypic_res[[i]], 'log2FoldChange', l2fcNA,
                             suffixes = c('',paste0('.',i)))
  colnames(l2fcmat) <- names(genotypic_res)
  
  padjmat <- mergeTwoRes(genotypic_res[[1]], genotypic_res[[2]], 
                         'padj', padjNA, suffixes = c('.1','.2'))
  if(length(genotypic_res)>2)
    for(i in 3:length(genotypic_res))
      padjmat <- mergeTwoRes(padjmat, genotypic_res[[i]], 'padj', padjNA,
                             suffixes = c('',paste0('.',i)))
  colnames(padjmat) <- names(genotypic_res)

  return(list(l2fcmat=l2fcmat, padjmat=padjmat))
}

degsMat <- function(genotypic_res, padj_thr_min=.1, log2fc_thr_min = 0) {
  
  repmat <- reportMat(genotypic_res, l2fcNA=0, padjNA=1)
  padjmat <- repmat$padjmat
  l2fcmat <- repmat$l2fcmat
  
  degsmat <- 1*(padjmat<padj_thr_min & l2fcmat>log2fc_thr_min) - 
    1*(padjmat<padj_thr_min & l2fcmat<(-log2fc_thr_min))
  degsmat[is.na(degsmat)] <- 0
  rownames(degsmat) <- rownames(l2fcmat)
  return(degsmat)
  
}

mypheatcts <- function(dds, genotypic_res, col_ann = c('sex','genotype','age'),
                       row_ann = c(), annot_colors = NA, annont_clustering = 0,
                       padj_thr_min=.1, log2fc_thr_min = 0, idx=NULL, 
                       markers=NULL, markers_padj_thr=1e-3, markers_l2fc_thr=1, 
                       cutree_rows=NA, ...) { #, ec_genes=NULL, ...) {
  
  # l2fcmat <- sapply(genotypic_res, function(x) x$log2FoldChange)
  # padjmat <- sapply(genotypic_res, function(x) x$padj)
  # degsmat <- 1*(padjmat<padj_thr_min & l2fcmat>log2fc_thr_min) - 1*(padjmat<padj_thr_min & l2fcmat<(-log2fc_thr_min))
  # rownames(l2fcmat) <- rownames(padjmat) <- rownames(degsmat) <- rownames(dds)
  
  degsmat <- degsMat(genotypic_res, padj_thr_min=padj_thr_min, 
                     log2fc_thr_min = log2fc_thr_min)
  dds <- dds[rownames(degsmat),]
  
  if( is.null(idx) ) {
    idx <- apply(degsmat!=0,1,any)
  }
  if(is.logical(idx)) idx <- rownames(dds)[which(idx)]
  
  stopifnot(require(RColorBrewer))
  stopifnot(require(pheatmap))

  # get the column annotation of selected fields from the dds
  ann_col <- lapply(col_ann, function(x) dds[[x]])
  ann_col <- data.frame(ann_col)
  rownames(ann_col) <- colnames(dds)
  colnames(ann_col) <- col_ann

  # order the dds using the column annotation
  ord <- do.call(order, ann_col)
  ann_col <- ann_col[ord,,drop=FALSE]
  dds <- dds[,ord]
  
  ann_row <- data.frame(apply(degsmat[idx,], 2, as.character), row.names = idx)
  # colnames(ann_row) <- substr(names(genotypic_res),1,10)
  ann_row <- ann_row[,rev(colnames(ann_row))]
  if(!is.null(markers)) {
    ind <- markers$padj < markers_padj_thr & markers$log2FoldChange > markers_l2fc_thr
    enriched <- rownames(markers)[which(ind)]
    ind <- markers$padj < markers_padj_thr & markers$log2FoldChange < (-markers_l2fc_thr)
    depleted <- rownames(markers)[which(ind)]
    ann_row$markers <- '0'
    ann_row$markers[rownames(ann_row) %in% enriched] <- '1'
    ann_row$markers[rownames(ann_row) %in% depleted] <- '-1'
  }
  ann_colors = lapply(colnames(ann_row), function(x) {
    cols <- rev(brewer.pal(n = 3, name = "RdYlBu"))
    names(cols) <- c('-1','0','1')
    return(cols)
  }); names(ann_colors) <- colnames(ann_row)
  if(length(row_ann>0)) {
    ann_row_add <- rowData(dds)[,rev(row_ann)]
    ann_row_add <- data.frame(ann_row_add)
    ann_row <- cbind(ann_row, ann_row_add[rownames(ann_row),])
  }
  if(!is.na(annot_colors[1])) {
    stopifnot(is.list(annot_colors))
    # remove already assigned colors
    ann_colors <- ann_colors[setdiff(names(ann_colors), names(annot_colors))]
    # annot_colors <- annot_colors[setdiff(names(annot_colors), names(ann_colors))]
    ann_colors <- c(ann_colors, annot_colors)
  }
  
  countsmat <- counts(dds, normalized = TRUE)[idx,]
  countsmat[countsmat==0] <- min(countsmat[countsmat>0])/exp(1)
  colnames(countsmat) <- rownames(ann_col)
  countsmat <- log(countsmat)
  countsmat <- (countsmat - apply(countsmat, 1, mean))/apply(countsmat, 1, sd)
  # gaps_col <- max(which(ann_col[,1] == ann_col[1,1]))
  gaps_col <- which(diff(as.numeric(as.factor(ann_col[,1]))) == 1)
  
  if(annont_clustering>0) {
    clustmat <- countsmat
    for(i in 1:annont_clustering)
      clustmat <- cbind(ann_row,clustmat)
    cluster_rows <- hclust(dist(clustmat))
  } else {
    cluster_rows <- TRUE
  }
  
  out <- pheatmap(countsmat, cluster_cols = FALSE, cluster_rows = cluster_rows,
                  cutree_rows = cutree_rows,
                  annotation_col = ann_col, annotation_row = ann_row,
                  annotation_colors = ann_colors, gaps_col = gaps_col, ...)
  # attach report to pheatmap output (giving appropriate colnames to the numeric matrix)
  
  condition_names <- sub(' ','',apply(ann_col,1,paste,collapse='_'))
  colnames(countsmat) <- paste(condition_names, 
                               colnames(countsmat), 
                               sep='_')
  report <- cbind(ann_row, countsmat)
  if(!is.na(cutree_rows)) {
    clst <- cutree(out$tree_row, k=cutree_rows)
    clst <- factor(clst, unique(clst[out$tree_row$order]))
    clst <- as.numeric(clst)
    report <- cbind(clusterID=clst, report)
    out$clusters_gene_list <- split(rownames(report), report$clusterID)
  }
  out$report <- report[out$tree_row$order,]
  out$countmeans <- t(apply(countsmat, 1, function(x) 
    tapply(x, condition_names, mean)))[out$tree_row$order,]
  out$clusterID <- out$report$clusterID
  out <- out
  
}

pheatcts <- function(ddsMN, res, l2fc_thr=c(0,Inf), padj_thr=.1, ...) {
  stopifnot(require(pheatmap))
  stopifnot(require(RColorBrewer))
  ann_col <- data.frame(genotype=as.character(ddsMN$genotype),
                        row.names = colnames(ddsMN)
  )
  idx <- res$padj<padj_thr & abs(res$log2FoldChange)>l2fc_thr[1]
  idx[is.na(idx)] <- FALSE
  ann_row <- data.frame(
    degs = with(subset(res, idx), 
                1*(padj<padj_thr & log2FoldChange>l2fc_thr[1]) - 1*(padj<padj_thr & log2FoldChange<l2fc_thr[1])),
    row.names = rownames(ddsMN[idx,]))
  countsmat <- counts(ddsMN[idx,], normalized = TRUE)
  countsmat[countsmat==0] <- min(countsmat[countsmat>0])/exp(1)
  countsmat_means <- apply(countsmat,1,mean)
  countsmat <- log2(countsmat/countsmat_means)
  cluster_rows <- hclust(dist(t(countsmat)))
  countsmat[countsmat>(+l2fc_thr[2])] <- +l2fc_thr[2]
  countsmat[countsmat<(-l2fc_thr[2])] <- -l2fc_thr[2]
  ann_colors = lapply(colnames(ann_row), function(x) {
    cols <- c("dodgerblue4", "grey50", "firebrick") # rev(brewer.pal(n = 3, name = "RdYlBu"))
    names(cols) <- c('-1','0','1')
    return(cols)
  }); names(ann_colors) <- colnames(ann_row)
  out <- pheatmap(countsmat, cluster_cols = FALSE, # scale="row", 
                  annotation_col = ann_col, annotation_row = ann_row,
                  annotation_colors = ann_colors, 
                  breaks = if(is.finite(l2fc_thr[2])) seq(-l2fc_thr[2], l2fc_thr[2], length.out=101) else NA,
                  # color = colorRampPalette(c("dodgerblue4", "grey50", "firebrick"))(100), ...)
                  color = colorRampPalette(c("dodgerblue4", "firebrick"))(100), ...)
  # color = colorRampPalette(c("blue", "red"))(100), ...)
}

kmpheatcts <- function(dds, genotypic_res, col_ann = c('sex','genotype','age'),
                       row_ann = c(), annot_colors = NA, annont_clustering = 0,
                       padj_thr_min=.1, log2fc_thr_min = 0, idx=NULL, 
                       markers=NULL, markers_padj_thr=1e-3, markers_l2fc_thr=1, 
                       cutree_rows=1, ...) { #, ec_genes=NULL, ...) {
  
  # l2fcmat <- sapply(genotypic_res, function(x) x$log2FoldChange)
  # padjmat <- sapply(genotypic_res, function(x) x$padj)
  # degsmat <- 1*(padjmat<padj_thr_min & l2fcmat>log2fc_thr_min) - 1*(padjmat<padj_thr_min & l2fcmat<(-log2fc_thr_min))
  # rownames(l2fcmat) <- rownames(padjmat) <- rownames(degsmat) <- rownames(dds)
  
  degsmat <- degsMat(genotypic_res, padj_thr_min=padj_thr_min, 
                     log2fc_thr_min = log2fc_thr_min)
  dds <- dds[rownames(degsmat),]
  
  if( is.null(idx) ) {
    idx <- apply(degsmat!=0,1,any)
  }
  if(is.logical(idx)) idx <- rownames(dds)[which(idx)]
  
  stopifnot(require(RColorBrewer))
  stopifnot(require(pheatmap))
  
  # get the column annotation of selected fields from the dds
  ann_col <- lapply(col_ann, function(x) dds[[x]])
  ann_col <- data.frame(ann_col)
  rownames(ann_col) <- colnames(dds)
  colnames(ann_col) <- col_ann
  
  # order the dds using the column annotation
  ord <- do.call(order, ann_col)
  ann_col <- ann_col[ord,,drop=FALSE]
  dds <- dds[,ord]
  
  ann_row <- data.frame(apply(degsmat[idx,], 2, as.character), row.names = idx)
  # colnames(ann_row) <- substr(names(genotypic_res),1,10)
  ann_row <- ann_row[,rev(colnames(ann_row))]
  if(!is.null(markers)) {
    ind <- markers$padj < markers_padj_thr & markers$log2FoldChange > markers_l2fc_thr
    enriched <- rownames(markers)[which(ind)]
    ind <- markers$padj < markers_padj_thr & markers$log2FoldChange < (-markers_l2fc_thr)
    depleted <- rownames(markers)[which(ind)]
    ann_row$markers <- '0'
    ann_row$markers[rownames(ann_row) %in% enriched] <- '1'
    ann_row$markers[rownames(ann_row) %in% depleted] <- '-1'
  }
  ann_colors = lapply(colnames(ann_row), function(x) {
    cols <- rev(brewer.pal(n = 3, name = "RdYlBu"))
    names(cols) <- c('-1','0','1')
    return(cols)
  }); names(ann_colors) <- colnames(ann_row)
  if(length(row_ann>0)) {
    ann_row_add <- rowData(dds)[,rev(row_ann)]
    ann_row_add <- data.frame(ann_row_add)
    ann_row <- cbind(ann_row, ann_row_add[rownames(ann_row),])
  }
  if(!is.na(annot_colors[1])) {
    stopifnot(is.list(annot_colors))
    # remove already assigned colors
    annot_colors <- annot_colors[setdiff(names(annot_colors), names(ann_colors))]
    ann_colors <- c(ann_colors, annot_colors)
  }
  
  countsmat <- counts(dds, normalized = TRUE)[idx,]
  countsmat[countsmat==0] <- min(countsmat[countsmat>0])/exp(1)
  colnames(countsmat) <- rownames(ann_col)
  countsmat <- log(countsmat)
  countsmat <- (countsmat - apply(countsmat, 1, mean))/apply(countsmat, 1, sd)
  # gaps_col <- max(which(ann_col[,1] == ann_col[1,1]))
  gaps_col <- which(diff(as.numeric(as.factor(ann_col[,1]))) == 1)
  
  if(annont_clustering>0) {
    clustmat <- countsmat
    for(i in 1:annont_clustering)
      clustmat <- cbind(ann_row,clustmat)
  } else {
    clustmat <- countsmat
  }
  
  km <- kmeans(clustmat, cutree_rows, iter.max = 100, nstart = 5)
  ord <- order(km$cluster)
  gaps_row <- which(diff(km$cluster[ord]) == 1) + 1
  out <- pheatmap(countsmat[ord,], cluster_cols = FALSE, cluster_rows = FALSE,
                  gaps_row = gaps_row,
                  annotation_col = ann_col, annotation_row = ann_row,
                  annotation_colors = ann_colors, gaps_col = gaps_col, ...)
  # attach report to pheatmap output (giving appropriate colnames to the numeric matrix)
  
  condition_names <- sub(' ','',apply(ann_col,1,paste,collapse='_'))
  colnames(countsmat) <- paste(condition_names, 
                               colnames(countsmat), 
                               sep='_')
  report <- cbind(ann_row, countsmat)
  clst <- km$cluster
  clst <- factor(clst, unique(clst[ord]))
  clst <- as.numeric(clst)
  report <- cbind(clusterID=clst, report)
  out$clusters_gene_list <- split(rownames(report), report$clusterID)
  out$report <- report[ord,]
  out$countmeans <- t(apply(countsmat, 1, function(x) 
    tapply(x, condition_names, mean)))[ord,]
  out$clusterID <- out$report$clusterID
  out <- out
  
}


fpkmPlotWhich <- function(dds, res, which,
                          l2fc_thr=c(0.6, 5), 
                          padj_thr=1e-1,
                          point_size=1,
                          background_alpha=.5,
                          gene_labels = '',
                          max.overlaps = 10,
                          robustSF = FALSE
) {
  
  dds <- dds[rownames(res),]
  de <- data.frame(res)
  if(robustSF) {
    dds <- estimateSizeFactors(dds)
  } else {
    if(is.null(sizeFactors(dds))) {
      print('sizeFactors not set, computing now for the given dds')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  de$fpkm <- apply(fpkm(dds),1,mean)
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  de$saturated <- satFCup | satFCdn
  de$saturated[is.na(de$saturated)] <- FALSE
  
  de$diffexpressed_tmp <- "NO"
  de$diffexpressed_tmp[which(de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min)] <- "UP"
  de$diffexpressed_tmp[which(de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min)] <- "DOWN"
  
  de$delabel <- NA
  # de$delabel[de$diffexpressed_tmp=='NO'] <- NA
  if(is.null(gene_labels)) {
    de[which,]$delabel <- rownames(de[which,])
    #de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  
  de$diffexpressed <- 'other'
  de[which,]$diffexpressed <- de[which,]$diffexpressed_tmp
  de$diffexpressed <- factor(de$diffexpressed, levels=c('DOWN','NO','UP','other'))
  
  n_gene_class <- table(de$diffexpressed)
  
  de$alpha <- background_alpha
  de[which,]$alpha <- 1
  
  de_which <- subset(de, diffexpressed != 'other')
  de_which$diffexpressed <- factor(de_which$diffexpressed, levels=c('DOWN','NO','UP'))
  assign('fpkmPlotWhichde', de_which, envir = globalenv())
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=log10(fpkm), 
                      color=diffexpressed, label=delabel, 
                      shape=saturated, alpha=alpha)) +
    geom_point(size=point_size) + 
    # stat_density_2d() + 
    # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
    geom_point(data = de_which, aes(x = log2FoldChange, y = log10(fpkm), 
                                    color = diffexpressed, 
                                    shape=saturated)) + 
    theme_minimal() +
    geom_text_repel(data = de_which, max.overlaps=max.overlaps) +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")"),
                           "other"
                         )[n_gene_class>0], 
                         type=c("DOWN"="blue", "NO"="black", "UP"="red", "other"="grey70")[n_gene_class>0]) + 
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), col="red")
  # geom_hline(yintercept=-log10(padj_thr_min), col="red")
  
}

fpkmPlot <- function(dds, res, 
                     l2fc_thr=c(0.6, 5), 
                     padj_thr=1e-1,
                     max.overlaps = 10,
                     gene_labels = NULL,
                     dotsize = 1,
                     robustSF = FALSE
) {
  
  dds <- dds[rownames(res),]
  de <- data.frame(res)
  if(robustSF) {
    dds <- estimateSizeFactors(dds)
  } else {
    if(is.null(sizeFactors(dds))) {
      print('sizeFactors not set, computing now for the given dds')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  de$fpkm <- apply(fpkm(dds),1,mean)
  # de <- de[!is.na(de$padj),]
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[which(de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min)] <- "UP"
  de$diffexpressed[which(de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min)] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c('DOWN', 'NO', 'UP'))
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  de$saturated <- satFCup | satFCdn
  de$saturated[is.na(de$saturated)] <- FALSE
  
  de$delabel <- NA
  if(is.null(gene_labels)) {
    de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  
  n_gene_class <- table(de$diffexpressed)
  assign('fpkmPlotde', de, envir = globalenv())
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=log10(fpkm), 
                      color=diffexpressed, label=delabel, shape=saturated)) +
    geom_point(size=dotsize) + 
    theme_minimal() +
    geom_text_repel(max.overlaps=max.overlaps) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), col="red") +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")")
                         )[n_gene_class>0], 
                         type=c("blue", "black", "red")[n_gene_class>0])
}

fpkmPlotAll <- function(ddsCnd, ddsRef, res, 
                        l2fc_thr=c(0.6, 5), 
                        padj_thr=1e-1,
                        max.overlaps = 10,
                        gene_labels = '',
                        dotsize = 1,
                        robustSF = FALSE
) {
  
  stopifnot(identical(rownames(ddsCnd),rownames(ddsRef)))
  
  de <- data.frame(
    baseMean = NA,
    log2FoldChange = log2(rowMeans(counts(ddsCnd, normalized=TRUE))/rowMeans(counts(ddsRef, normalized=TRUE))),
    lfcSE=NA, stat=NA, pvalue=NA, padj=NA,
    row.names = rownames(ddsCnd)
  )
  de[rownames(res),] <- data.frame(res)
  if(robustSF) {
    ddsCnd <- estimateSizeFactors(ddsCnd)
  } else {
    if(is.null(sizeFactors(ddsCnd))) {
      print('sizeFactors not set, computing now for the given ddsCnd')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  de$fpkm <- apply(fpkm(ddsCnd),1,mean)
  
  # filter for zero counts in both conditions
  filt <- rowMeans(counts(ddsCnd)) > 10 | rowMeans(counts(ddsRef)) > 10
  de <- de[filt,]

  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  
  de$diffexpressed <- "NO"
  de$diffexpressed[which(de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min)] <- "UP"
  de$diffexpressed[which(de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min)] <- "DOWN"
  de$diffexpressed <- factor(de$diffexpressed, levels = c('DOWN', 'NO', 'UP'))
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  de$saturated <- satFCup | satFCdn
  de$saturated[is.na(de$saturated)] <- FALSE
  
  de$delabel <- NA
  if(is.null(gene_labels)) {
    de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }
  
  n_gene_class <- table(de$diffexpressed)
  
  assign('fpkmPlotAllde', de, envir = globalenv())
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=log10(fpkm), 
                      color=diffexpressed, label=delabel, shape=saturated)) +
    geom_point(size=dotsize) + 
    theme_minimal() +
    geom_text_repel(max.overlaps=max.overlaps) +
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), col="red") +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")")
                         )[n_gene_class>0], 
                         type=c("blue", "black", "red")[n_gene_class>0])
}

fpkmPlotAllWhich <- function(ddsCnd, ddsRef, res, which,
                          l2fc_thr=c(0.6, 5), 
                          padj_thr=1e-1,
                          point_size=1,
                          gene_labels = '',
                          background_alpha=.5,
                          max.overlaps = 10,
                          robustSF = FALSE
) {
  
  stopifnot(identical(rownames(ddsCnd),rownames(ddsRef)))
  
  de <- data.frame(
    baseMean = NA,
    log2FoldChange = log2(rowMeans(counts(ddsCnd, normalized=TRUE))/rowMeans(counts(ddsRef, normalized=TRUE))),
    lfcSE=NA, stat=NA, pvalue=NA, padj=NA,
    row.names = rownames(ddsCnd)
  )
  de[rownames(res),] <- data.frame(res)
  if(robustSF) {
    ddsCnd <- estimateSizeFactors(ddsCnd)
  } else {
    if(is.null(sizeFactors(ddsCnd))) {
      print('sizeFactors not set, computing now for the given ddsCnd')
    } else {
      print('Using pre-computed sizeFactors, if you want to re-evaluate them put robustSF to TRUE')
    }
  }
  de$fpkm <- apply(fpkm(ddsCnd),1,mean)
  
  # filter for zero counts in both conditions
  filt <- rowMeans(counts(ddsCnd)) > 10 | rowMeans(counts(ddsRef)) > 10 |
    rownames(ddsCnd) %in% which
  de <- de[filt,]
  
  l2fc_thr_min <- l2fc_thr[1]
  l2fc_thr_max <- l2fc_thr[2]
  padj_thr_min <- padj_thr[1]
  
  satFCup <- de$log2FoldChange > l2fc_thr_max
  de$log2FoldChange[satFCup] <- l2fc_thr_max
  satFCdn <- de$log2FoldChange < -l2fc_thr_max
  de$log2FoldChange[satFCdn] <- -l2fc_thr_max
  de$saturated <- satFCup | satFCdn
  de$saturated[is.na(de$saturated)] <- FALSE
  
  de$diffexpressed_tmp <- "NO"
  de$diffexpressed_tmp[which(de$log2FoldChange > l2fc_thr_min & de$padj < padj_thr_min)] <- "UP"
  de$diffexpressed_tmp[which(de$log2FoldChange < -l2fc_thr_min & de$padj < padj_thr_min)] <- "DOWN"
  
  de$delabel <- NA
  # de$delabel[de$diffexpressed_tmp=='NO'] <- NA
  if(is.null(gene_labels)) {
    de[which,]$delabel <- rownames(de[which,])
    #de$delabel[de$diffexpressed != "NO"] <- rownames(de)[de$diffexpressed != "NO"]
  } else {
    idx <- rownames(de) %in% gene_labels
    de$delabel[idx] <- rownames(de)[idx]
  }

  de$diffexpressed <- 'other'
  de[which,]$diffexpressed <- de[which,]$diffexpressed_tmp
  de$diffexpressed <- factor(de$diffexpressed, levels=c('DOWN','NO','UP','other'))
  de$alpha <- background_alpha
  de[which,]$alpha <- 1
  
  de_which <- subset(de, diffexpressed != 'other')
  de_which$diffexpressed <- factor(de_which$diffexpressed, levels=c('DOWN','NO','UP'))
  assign('fpkmPlotAllWhichde', de_which, envir = globalenv())

  n_gene_class <- table(de$diffexpressed)
  
  library(ggrepel)
  # plot adding up all layers we have seen so far
  ggplot(data=de, aes(x=log2FoldChange, y=log10(fpkm), 
                      color=diffexpressed, label=delabel, 
                      shape=saturated, alpha=alpha)) +
    geom_point(size=point_size) + 
    # scale_color_discrete(name = "Gene class", 
    #                      labels = c(
    #                        paste0("Depleted (",n_gene_class[1],")"), 
    #                        paste0("Neutral (",n_gene_class[2],")"), 
    #                        paste0("Enriched (",n_gene_class[3],")"),
    #                        "other"
    #                      )[n_gene_class>0], 
    #                      type=c("DOWN"="blue", "NO"="black", "UP"="red", "other"="grey70")[n_gene_class>0]) + 
    # stat_density_2d() + 
    # stat_density_2d(aes(fill = ..level..), geom="polygon") + 
    geom_point(data = de_which, 
               aes(x = log2FoldChange, 
                   y = log10(fpkm), 
                   color = diffexpressed, 
                   shape=saturated)) + 
    theme_minimal() +
    geom_text_repel(data = de_which, 
                    max.overlaps=max.overlaps) +
    scale_color_discrete(name = "Gene class", 
                         labels = c(
                           paste0("Depleted (",n_gene_class[1],")"), 
                           paste0("Neutral (",n_gene_class[2],")"), 
                           paste0("Enriched (",n_gene_class[3],")"),
                           'other'
                           )[n_gene_class>0], 
                         type=c("DOWN"="blue", "NO"="black", 
                                "UP"="red", "other"="grey70")[n_gene_class>0]) + 
    geom_vline(xintercept=c(-l2fc_thr_min, l2fc_thr_min), col="red") +
    guides(alpha = 'none')
  # geom_hline(yintercept=-log10(padj_thr_min), col="red")
  
}
