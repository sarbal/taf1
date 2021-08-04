 
# Pre-filtering 
load("deg.wilcox.adj.Rdata")
deg = deg.wilcox.new
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
genes.ranks.pre = t(genesets.fcranks.all[filt,])

# Combat adjusted results 
load("deg.wilcox.combat.adj.Rdata") 
deg = deg.wilcox.combat.adj
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.comb = t(genesets.fcranks.all[filt,])


# DESeq2 batch aware
load("deg.deseq.batch.Rdata")
deg = deg.deseq.batch
genesets.fcranks.all = sapply(7:63, function(i) rank(deg[[i]][,2] ) )    
genesets.padj.all = sapply(7:63, function(i) rank( log10(deg[[i]][,6] ) )   )  
genesets.both.all = sapply(7:63, function(i) rank(rank( log10(deg[[i]][,6] ) ) +  rank(deg[[i]][,2] ) )   )  

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.deseqBfc = t(genesets.fcranks.all[filt,])
genes.ranks.deseqBp = t(genesets.padj.all[filt,])
genes.ranks.deseqBboth = t(genesets.both.all[filt,])


# DESeq2 filtered ?
load("deg.deseq.filt.Rdata")
deg = deg.deseq.filt
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
genesets.padj.all = sapply(1:63, function(i) rank( log10(deg[[i]][,6] ) )   )  
genesets.both.all = sapply(1:63, function(i) rank(rank( log10(deg[[i]][,6] ) ) +  rank(deg[[i]][,2] ) )   )  

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.deseqfc = t(genesets.fcranks.all[filt,])
genes.ranks.deseqp = t(genesets.padj.all[filt,])
genes.ranks.deseqboth = t(genesets.both.all[filt,])


# DESeq2 default 
load("deg.deseq.Rdata")
deg = deg.deseq 
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
genesets.padj.all = sapply(1:63, function(i) rank( log10(deg[[i]][,6] ) )   )  
genesets.both.all = sapply(1:63, function(i) rank(rank( log10(deg[[i]][,6] ) ) +  rank(deg[[i]][,2] ) )   )  

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.deseq1fc = t(genesets.fcranks.all[filt,])
genes.ranks.deseq1p = t(genesets.padj.all[filt,])
genes.ranks.deseq1both = t(genesets.both.all[filt,])


