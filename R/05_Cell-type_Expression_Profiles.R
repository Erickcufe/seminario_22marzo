load("data/ADgenes.rda")
load("data/geneAnno2.rda")
targetname <- "AD"
targetgene <- ADhgnc
cellexp <- read.table("data/DER-20_Single_cell_expression_processed_TPM_backup.tsv",header=T,fill=T)
cellexp[1121,1] <- cellexp[1120,1]
cellexp <- cellexp[-1120,]
rownames(cellexp) <- cellexp[,1]
cellexp <- cellexp[,-1]
datExpr <- scale(cellexp,center=T, scale=F)
datExpr <- datExpr[,789:ncol(datExpr)]

# Extract cellular expression profiles of AD risk genes
exprdat <- apply(datExpr[match(targetgene, rownames(datExpr)),],2,mean,na.rm=T)
dat <- data.frame(Group=targetname, cell=names(exprdat), Expr=exprdat)
dat$celltype <- unlist(lapply(strsplit(dat$cell, split="[.]"),'[[',1))
dat <- dat[-grep("Ex|In",dat$celltype),]
dat$celltype <- gsub("Dev","Fetal",dat$celltype)
dat$celltype <- factor(dat$celltype, levels=c("Neurons","Astrocytes","Microglia","Endothelial",
                                             "Oligodendrocytes","OPC","Fetal"))
pdf(file ="images/singlecell_expression_ADgenes.pdf")
ggplot(dat,aes(x=celltype, y=Expr, fill=celltype)) +
  ylab("Normalized expression") + xlab("") +
  geom_violin() +
  theme(axis.text.x=element_text(angle = 90, hjust=1)) +
  theme(legend.position="none") +
ggtitle(paste0("Cellular expression profiles of AD risk genes"))
dev.off()
