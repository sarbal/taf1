load("gene_annotations_v22.Rdata")
f.e = !is.na(attr$entrezID )
load("genes.Rdata")
m = match( genes[,1], attr$entrezID)
f.g = !is.na(m)
f.a = m[f.g]
 

m =  match(rownames(deg[[1]]), attr$ensemblID[f.a])
f.m = !is.na(m)
f.am = m[f.m]


N = length(f.am)



load(file="res.prec.deseq.test2.min5.Rdata")
down100.f5 = sapply(1:6, function(i)  match(names(res.prec[[i]][[2]][100,][res.prec[[i]][[2]][100,]>0] ), attr$name[f.a][f.am]) )
up100.f5 = sapply(1:6, function(i)  match(names(res.prec[[i]][[1]][100,][res.prec[[i]][[1]][100,]>0] ), attr$name[f.a][f.am]) )
down200.f5  = sapply(1:6, function(i)  match(names(res.prec[[i]][[2]][200,][res.prec[[i]][[2]][200,]>0] ), attr$name[f.a][f.am]) )
up200.f5 = sapply(1:6, function(i)  match(names(res.prec[[i]][[1]][200,][res.prec[[i]][[1]][200,]>0] ), attr$name[f.a][f.am]) )


load(file="res.prec.deseq.test2.min6.Rdata")
down100.f6 = sapply(1:6, function(i)  match(names(res.prec[[i]][[2]][100,][res.prec[[i]][[2]][100,]>0] ), attr$name[f.a][f.am]) )
up100.f6 = sapply(1:6, function(i)  match(names(res.prec[[i]][[1]][100,][res.prec[[i]][[1]][100,]>0] ), attr$name[f.a][f.am]) )
down200.f6  = sapply(1:6, function(i)  match(names(res.prec[[i]][[2]][200,][res.prec[[i]][[2]][200,]>0] ), attr$name[f.a][f.am]) )
up200.f6 = sapply(1:6, function(i)  match(names(res.prec[[i]][[1]][200,][res.prec[[i]][[1]][200,]>0] ), attr$name[f.a][f.am]) )



cc = rep(100,63)
cc2 = rep(200,63)
od100f6.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(down100.f6[[j]] ))) )
od100f5.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(down100.f5[[j]] ))) )
od200f6.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(down200.f6[[j]] ))) )
od200f5.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(down200.f5[[j]] ))) )
ou100f6.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(up100.f6[[j]] ))) )
ou100f5.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(up100.f5[[j]] ))) )
ou200f6.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(up200.f6[[j]] ))) )
ou200f5.n = sapply(mat.list, function(i) max(sapply(i, function(j)  length(up200.f5[[j]] ))) )



max.ns=
cbind(	cc,
        cc,
	od100f6.n,
	od100f5.n,

    cc2,
	cc2,
	od200f6.n,
	od200f5.n,

	cc,
	cc,
	ou100f6.n,
	ou100f5.n,

    cc2,
	cc2,
	ou200f6.n,
	ou200f5.n
)

 

btests = list()
btests.adj = list()

for(x in 1:200){

	temp.mat     = matrix(NA, ncol=102, nrow=100)
        temp.mat.adj = matrix(NA, ncol=102, nrow=100)
	p = x/N
	b.tests = lapply(1:100, function(i) pbinom(-1:i, i, p, lower.tail=F))
	b.tests.adj = lapply(1:100, function(i) p.adjust(pbinom(-1:i, i, p, lower.tail=F), n=N ))
        for(i in 1:100) {
		temp.mat[i,1:length(-1:i)] = b.tests[[i]]
		temp.mat.adj[i,1:length(-1:i)] = b.tests.adj[[i]]
	}

        btests[[x]] = temp.mat
        btests.adj[[x]] = temp.mat.adj

}


test.all = lapply(1:16, function(j) t(sapply(1:length(mat.list), function(i) -log10(btests[[ max.ns[i,j]]][1:6,1:7][n.s[i],][-1]) ) ) )
test.masks = lapply(1:16, function(j) test.all[[j]][,1:6] >= -log10(0.05/N) )
test.masked = lapply(1:16, function(j) t(counts.mats[[j]] ))



# Pre-filtering 
load("deg.wilcox.adj.Rdata")
deg = deg.wilcox.adjusted
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
genesets.padj.all = sapply(1:63, function(i) rank( log10(deg[[i]][,6] ) )   )  
genesets.both.all = sapply(1:63, function(i) rank(rank( log10(deg[[i]][,6] ) ) +  rank(deg[[i]][,2] ) )   )  

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.pre = t(genesets.fcranks.all[filt,])

# Combat adjusted results 
load("deg.wilcox.combat.adj.Rdata") 
deg = deg.wilcox.combat.adj
genesets.fcranks.all = sapply(1:63, function(i) rank(deg[[i]][,2] ) )    
genesets.padj.all = sapply(1:63, function(i) rank( log10(deg[[i]][,6] ) )   )  
genesets.both.all = sapply(1:63, function(i) rank(rank( log10(deg[[i]][,6] ) ) +  rank(deg[[i]][,2] ) )   )  

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.comb = t(genesets.fcranks.all[filt,])

# Combat  results 
load("deg.wilcox.combat.Rdata") 
deg = deg.wilcox.combat
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
genes = rownames(deg[[63]])
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
genes = rownames(deg[[63]])

xist =  which(  rownames(deg[[63]])  == "ENSG00000229807")
cacna =  which(  rownames(deg[[63]])  == "ENSG00000100346")
igfbp =   which(  rownames(deg[[63]])  == "ENSG00000146674")
filt = c(xist, cacna,igfbp)
genes.ranks.deseq1fc = t(genesets.fcranks.all[filt,])
genes.ranks.deseq1p = t(genesets.padj.all[filt,])
genes.ranks.deseq1both = t(genesets.both.all[filt,])



xist =  which(  commongenes  == "ENSG00000229807")
cacna =  which( commongenes == "ENSG00000100346")
igfbp =   which( commongenes == "ENSG00000146674")
filt = c(xist, cacna,igfbp)






filesdir = dir()[9:14]
for(i in filesdir) { load(i); genesets.all[[i]] = genesets.both.all }
for(i in filesdir[1:3]) { load(i); rownames(genesets.both.all) = genes; genesets.all[[i]] = genesets.both.all }
freqgenes = plyr::count(unlist(lapply(1:6, function(i) rownames(genesets.all[[i]] ) ) ) )
commongenes = freqgenes[freqgenes[,2] == 6 ,1]

m = lapply(1:6, function(i)  match(commongenes, rownames(genesets.all[[i]]  )  ) )
for(i in 1:6){ genesets.all[[i]] = genesets.all[[i]][m[[i]],] }

genesets.rank = lapply(1:6, function(i) apply(genesets.all[[i]], 2, rank ) )

xist =  which(  commongenes  == "ENSG00000229807")
cacna =  which( commongenes == "ENSG00000100346")
igfbp =   which( commongenes == "ENSG00000146674")
taf1 =  which( commongenes == "ENSG00000147133")
filt = c(xist, cacna,igfbp, taf1)


for(i in filesdir) { load(i); genesets.all[[i]] = genesets.fcranks.all }
for(i in filesdir[1:3]) { load(i); rownames(genesets.fcranks.all) = genes; genesets.all[[i]] = genesets.fcranks.all }
for(i in 1:6){ genesets.all[[i]] = genesets.all[[i]][m[[i]],] }
genesets.rank2 = lapply(1:6, function(i) apply(genesets.all[[i]], 2, rank ) )

genesets.rank.filt = list() 
 
for(i in c(1:3,5:6)) { 
genesets.rank.filt[[i]] =  cbind(t(genesets.rank[[i]][filt,]) , t(genesets.rank2[[i]][filt,])  )   
} 

fill = rep(NA, 6 ) 
fill2 = cbind(fill, fill,fill,fill)
i = 4 
genesets.rank.filt[[i]] =  cbind(rbind(fill2,t(genesets.rank[[i]][filt,])) , t(genesets.rank2[[i]][filt, ])  )   

i = 1 
genesets.rank.filt[[i]] =  cbind(rbind(fill2,t(genesets.rank[[i]][filt,])) , rbind(fill2,t(genesets.rank2[[i]][filt,])) )   


load("generic/ranks2.Rdata") 
load("generic/ranks.filt2.Rdata") 
 

 genesets.rank.filt[[8]] =  t(ranks.candidates.filt2)
 genesets.rank.filt[[7]]=  t(ranks.candidates2)



save( mat.list, n.s,genesets.rank2, filt, filesdir, genesets.rank, file="genesets.rank.all.Rdata")
save( mat.list, n.s, filesdir, genesets.rank.filt, file="genesets.rank.filt.Rdata")

 


j = 2 + 3 
i = 3 
plot( log10(genesets.rank.filt[[i]][,j]) ~ n.s , ylim=c(0,5)) 


 
beanplot( genesets.rank.filt[[i]][,j] ~ n.s   , pch=19, xlab="Number of probands", ylab="Rank of CACNA1I", col=list("plum"), side="f", bw="nrd0")
  

n.s2 = c(n.s,n.s+0.05, n.s+0.1,n.s+0.15,n.s+0.2,n.s+0.25,n.s+0.3, n.s+0.35)
n.s3 = unlist(lapply(1:8, function(i) rep(i,63) ) )


temp = c(genesets.rank.filt[[1]][,j], 
	genesets.rank.filt[[2]][,j],
	genesets.rank.filt[[3]][,j],
	genesets.rank.filt[[4]][,j], 
	genesets.rank.filt[[5]][,j], 
	genesets.rank.filt[[6]][,j], 
	genesets.rank.filt[[7]][,j], 
	genesets.rank.filt[[8]][,j]) 
labels = c("DESeq2 batch", "DESeq2 filt", "DESeq2", "Wilcox + adj", "Wilcox + adj + combat", "Wilcox + combat", "Recurrence", "Recurrence + filt") 
temp[is.na(temp)] = jitter(rep(13000, 12) )

beanplot( temp ~ n.s2   , pch=19, xlab="Number of probands", ylab="Rank of CACNA1I", 
	col=list("purple", "grey", "plum", "magenta", "pink", "darkmagenta", "red"), side="f", bw="nrd0")

 
beeswarm( log10(temp) ~ n.s2   , pch=19, xlab="Number of probands", ylab="Rank of CACNA1I", pwcol=n.s2)

beeswarm( log10(temp) ~ n.s2   , pch=19, xlab="Number of probands", ylab="Rank of CACNA1I", pwcol=n.s3)

j = 2   
par(mfrow=c(3,3))
for(i in 1:8){ 
beeswarm( log10(genesets.rank.filt[[i]][,j]) ~ n.s   , main=labels[i], ylim=c(0,5), pch=19, xlab="Number of probands",
 ylab="Rank of CACNA1I" )
} 


 j = 1   
par(mfrow=c(3,3))
for(i in 1:8){ 
beeswarm( log10(genesets.rank.filt[[i]][,j]) ~ n.s   , main=labels[i], ylim=c(0,5), pch=19, xlab="Number of probands",
 ylab="Rank of XIST" )
} 


  j = 2   
templist = list(genesets.rank.filt[[1]][,j], 
	genesets.rank.filt[[2]][,j],
	genesets.rank.filt[[3]][,j],
	genesets.rank.filt[[4]][,j], 
	genesets.rank.filt[[5]][,j], 
	genesets.rank.filt[[6]][,j], 
	genesets.rank.filt[[7]][,j], 
	genesets.rank.filt[[8]][,j]) 
n.s4 = lapply(1:8, function(i) n.s)


 j = 2   
par(mfrow=c(3,3))
for(i in 1:8){ 
	templist = lapply(2:6, function(k) log10(genesets.rank.filt[[i]][n.s==k,j]) )
	vioplot.2( templist , plotpoints=T, main=labels[i], colin=cols6  )
} 




pdf("output.power.pdf", height=12, width=12)
cols_range = cividis(10)
 j = 2   
par(mfrow=c(3,3))
for(i in 1:8){ 
	templist = lapply(2:6, function(k) log10(genesets.rank.filt[[i]][n.s==k,j]) )
	names(templist) = 2:6 
	vioplot.2( templist, xlab="Number of probands",ylab="Rank of CACNA1I" , 
		colin=cols_range, plotpoints=T , main=labels[i], yrange=c(0,5 ) ) 
} 



 j = 3   
par(mfrow=c(3,3))
for(i in 1:8){ 
	templist = lapply(2:6, function(k) log10(genesets.rank.filt[[i]][n.s==k,j]) )
	names(templist) = 2:6 
	vioplot.2( templist, xlab="Number of probands",ylab="Rank of IGFBP3" , 
		colin=cols_range, plotpoints=T , main=labels[i], yrange=c(0,5 ) ) 
} 

 j = 1   
par(mfrow=c(3,3))
for(i in 1:8){ 
	templist = lapply(2:6, function(k) log10(genesets.rank.filt[[i]][n.s==k,j]) )
	names(templist) = 2:6 
	vioplot.2( templist, xlab="Number of probands",ylab="Rank of XIST" , 
		colin=cols_range, plotpoints=T , main=labels[i], yrange=c(0,5 ) ) 
} 


 j = 4   
par(mfrow=c(3,3))
for(i in 1:8){ 
	templist = lapply(2:6, function(k) log10(genesets.rank.filt[[i]][n.s==k,j]) )
	names(templist) = 2:6
	vioplot.2( templist, xlab="Number of probands",ylab="Rank of TAF1" , 
		colin=cols_range, plotpoints=T , main=labels[i], yrange=c(0,5 ) ) 
} 


dev.off() 





