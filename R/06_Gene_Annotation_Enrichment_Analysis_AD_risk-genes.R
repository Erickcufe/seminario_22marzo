library(ggpubr)

pdf("images/GO_enrichment.pdf",width=15,height=8)
plot_barplot = function(dbname,name,color){
  input = read.delim(paste0(dbname,".txt"),header=T)
  input = input[,c(-1,-10,-11)]
  input = unique(input)
  input$FDR = p.adjust(exp(input$logP))
  input_sig = input[input$FDR < 0.1,]
  input_sig$FDR = -log10(input_sig$FDR)
  input_sig = input_sig[order(input_sig$FDR),]
  p = ggbarplot(input_sig, x = "Term", y = "FDR", fill = color, color = "white", sort.val = "asc", ylab = expression(-log[10](italic(FDR))), xlab = paste0(name," Terms"), rotate = TRUE, label = paste0(input_sig$Target.Genes.in.Term,"/",input_sig$Genes.in.Term), font.label = list(color = "white", size = 9), lab.vjust = 0.5, lab.hjust = 1)
  p = p+geom_hline(yintercept = -log10(0.05), linetype = 2, color = "lightgray")
  return(p)
}
p1 <- plot_barplot("data/homer_analysis_enrichment/biological_process","GO Biological Process","#00AFBB")
p2 <- plot_barplot("data/homer_analysis_enrichment/kegg","KEGG","#E7B800")
p3 <- plot_barplot("data/homer_analysis_enrichment/reactome","Reactome","#FC4E07")
ggarrange(p1, p2, p3, labels = c("A", "B", "C"), ncol = 2, nrow = 2)
dev.off()

plot_barplot("data/homer_analysis_enrichment/cosmic","Cosmic","#FC4E07")

ver <- read.delim("data/homer_analysis_enrichment/cosmic.txt")
