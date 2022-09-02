cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  rm(Hat)
  gc()
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}

plotCounts1 <- function(ddsmat, mygene, group, gene_name) {
  ind <- which(row.names(counts(ddsmat)) == mygene)
  mydata <- counts(ddsmat, normalized = TRUE)[ind,]
  mydata <- as.data.frame(mydata)
  mydata$group <- colData(ddsmat)[match(row.names(mydata), ddsmat$Sample), group] 
  gene_idx <- which(row.names(rowData(ddsmat)) == mygene)
  gname <- rowData(ddsmat)[gene_idx, gene_name]
  ggplot(mydata, aes(x = group, y = mydata)) + geom_beeswarm(cex = 3) + ggtitle(paste(mygene, " : ", gname))
}

checkY <- function(idx, region) {
  coldata <- coldata.ord[idx,]
  countdata <- countdata[,idx]
  
  geneid <- row.names(countdata)
  ens.mm.v97 <- AnnotationHub()[["AH73905"]]
  chr <- mapIds(ens.mm.v97, keys = geneid, keytype = "GENEID", column = "SEQNAME")
  all(names(chr) == row.names(countdata))
  y <- countdata[which(chr == "Y"), ]
  #y2 <- countdata.drg.t12[which(row.names(countdata.drg.t12) %in% ygenes$V1), ]
  s <- colSums(y)
  smelt <- melt(s)
  smelt$group <- coldata[match(row.names(smelt), coldata$Sample), "group"]
  smelt$gender <- coldata[match(row.names(smelt), coldata$Sample), "Gender"]
  smelt$sample <- row.names(smelt)
  #pdf(here(region, "troubleshooting", "Counts_on_Y.pdf"))
  ggplot(smelt, aes(x = sample, y = value, col = gender)) + geom_point(aes(shape = group), size = 3) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) + ggtitle(paste0(region, " Total counts on the Y chromosome"))
  #dev.off()
}