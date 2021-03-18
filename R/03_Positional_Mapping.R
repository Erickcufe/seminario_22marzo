library(GenomicRanges)
load("data/AD_credibleSNP.rda")

# Load promoter and exonic region and generate a GRange object.
exon <- read.table("data/Gencode19_exon.bed")
exonranges <- GRanges(exon[,1],IRanges(exon[,2], exon[,3]), gene = exon[,4])
promoter <- read.table("data/Gencode19_promoter.bed")
promoterranges <- GRanges(promoter[,1], IRanges(promoter[,2], promoter[,3]), gene = promoter[,4])

# Overlap credible SNPs with exonic regions.
olap <- findOverlaps(credranges, exonranges)
credexon <- credranges[queryHits(olap)]
mcols(credexon) <- cbind(mcols(credexon), mcols(exonranges[subjectHits(olap)]))

# Overlap credible SNPs with promoter regions.
olap <- findOverlaps(credranges, promoterranges)
credpromoter <- credranges[queryHits(olap)]
mcols(credpromoter) <- cbind(mcols(credpromoter), mcols(promoterranges[subjectHits(olap)]))

# Link SNPs to their putative target genes using chromatin interactions.
# Load Hi-C dataset and generate a GRange object.
hic <- read.table("data/Promoter-anchored_chromatin_loops.bed", skip=1)
colnames(hic) <- c("chr", "TSS_start", "TSS_end", "Enhancer_start", "Enhancer_end")
hicranges <- GRanges(hic$chr, IRanges(hic$TSS_start, hic$TSS_end), enhancer=hic$Enhancer_start)
olap <- findOverlaps(hicranges, promoterranges)
hicpromoter <- hicranges[queryHits(olap)]
mcols(hicpromoter) <- cbind(mcols(hicpromoter), mcols(promoterranges[subjectHits(olap)]))
hicenhancer <- GRanges(seqnames(hicpromoter), IRanges(hicpromoter$enhancer, hicpromoter$enhancer+10000),
                       gene=hicpromoter$gene)

# Overlap credible SNPs with Hi-C GRange object.
olap <- findOverlaps(credranges, hicenhancer)
credhic <- credranges[queryHits(olap)]
mcols(credhic) <- cbind(mcols(credhic), mcols(hicenhancer[subjectHits(olap)]))

# Compile AD candidate genes defined by positional mapping and chromatin interaction profiles.
### The resulting candidate genes for AD:
ADgenes <- Reduce(union, list(credhic$gene, credexon$gene, credpromoter$gene))
### to convert Ensembl Gene ID to HGNC symbol
load("data/geneAnno2.rda")
ADhgnc <- geneAnno1[match(ADgenes, geneAnno1$ensembl_gene_id), "hgnc_symbol"]
ADhgnc <- ADhgnc[ADhgnc!=""]
save(ADgenes, ADhgnc, file="data/ADgenes.rda")
write.table(ADhgnc, file="data/ADgenes.txt", row.names=F, col.names=F, quote=F, sep="\t")


