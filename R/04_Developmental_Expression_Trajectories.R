library(reshape)
library(ggplot2)
library(GenomicRanges)
library(biomaRt)
library("WGCNA")

# Process expression and meta data
datExpr <- read.csv("data/gene_array_matrix_csv/expression_matrix.csv", header = FALSE)
datExpr <- datExpr[,-1]
datMeta <- read.csv("data/gene_array_matrix_csv/columns_metadata.csv")
datProbes <- read.csv("data/gene_array_matrix_csv/rows_metadata.csv")
datExpr <- datExpr[datProbes$ensembl_gene_id!="",]
datProbes <- datProbes[datProbes$ensembl_gene_id!="",]
datExpr.cr <- collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id,
                          rowID= rownames(datExpr))
datExpr <- datExpr.cr$datETcollapsed
gename <- data.frame(datExpr.cr$group2row)
rownames(datExpr) <- gename$group

# Specify developmental stages.
datMeta$Unit <- "Postnatal"
idx <- grep("pcw", datMeta$age)
datMeta$Unit[idx] <- "Prenatal"
idx <- grep("yrs", datMeta$age)
datMeta$Unit[idx] <- "Postnatal"
datMeta$Unit <- factor(datMeta$Unit, levels=c("Prenatal", "Postnatal"))

# Select Cortical regions
datMeta$Region <- "SubCTX"
r <- c("A1C", "STC", "ITC", "TCx", "OFC", "DFC",
       "VFC", "MFC", "M1C", "S1C", "IPC", "M1C-S1C",
       "PCx", "V1C", "Ocx")
datMeta$Region[datMeta$structure_acronym %in% r] <- "CTX"
datExpr <- datExpr[,which(datMeta$Region=="CTX")]
datMeta <- datMeta[which(datMeta$Region=="CTX"),]
save(datExpr, datMeta, file="data/devExpr.rda")

# Extract developmental expression profiles of AD risk genes.
load("data/ADgenes.rda")
exprdat <- apply(datExpr[match(ADgenes, rownames(datExpr)),], 2, mean, na.rm=T)
dat <- data.frame(Region=datMeta$Region, Unit=datMeta$Unit, Expr=exprdat)

# Compare prenatal versus postnatal expression levels of AD risk genes.
pdf(file="images/developmental_expression.pdf")
ggplot(dat,aes(x=Unit, y=Expr, fill=Unit, alpha=Unit)) +
  ylab("Normalized expression") + geom_boxplot(outlier.size = NA) +
  ggtitle("Brain Expression") + xlab("") +
  scale_alpha_manual(values=c(0.2, 1)) + theme_classic() +
  theme(legend.position="na")
dev.off()



