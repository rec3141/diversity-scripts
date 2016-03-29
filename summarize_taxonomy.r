library(RColorBrewer)
library(gplots)

ab <- read.table("arctic_18S.taxonomy.summary.txt",header=T,check.names=F,sep="\t")
ac <- read.table("arctic_16S.taxonomy.summary.txt",header=T,check.names=F,sep="\t")

	
pdf(file="arctic_16S_18S.taxa.pdf",width=60,height=12)

for (taxl in 2:6) {
	ad <- as.matrix(subset(ab,taxlevel==taxl)[,6:ncol(ab)])
	colnames(ad) <- colnames(ab)[6:ncol(ab)]
	rownames(ad) <- subset(ab,taxlevel==taxl,select="taxon")[,1]

	ae <- as.matrix(subset(ac,taxlevel==taxl)[,6:ncol(ac)])
	colnames(ae) <- colnames(ac)[6:ncol(ac)]
	rownames(ae) <- subset(ac,taxlevel==taxl,select="taxon")[,1]

	af <- rbind(ad,ae)
	
	af <- af[rowSums(af)>sum(af)/1e7,]
	heatmap.2(t(af),scale="none",margins=c(20,15),col=brewer.pal(11,"RdBu"),trace="none",main="raw counts")
	heatmap.2((((t(prop.table(af,2))))),scale="none",margins=c(20,15),col=brewer.pal(11,"RdBu"),trace="none",main="sample %")
	plot(cmdscale(dist((((t(prop.table(af,2))))))))
	heatmap.2(sqrt(sqrt((t(prop.table(af,2))))),scale="none",margins=c(20,15),col=brewer.pal(11,"RdBu"),trace="none",main="sqrt sqrt sample %",cexCol=0.5,cexRow=0.5)
}
dev.off()






#	heatmap.2(t(ad),scale="col",margins=c(20,15),col=brewer.pal(11,"RdBu"),trace="none",main="taxon stdev normalized")
#	heatmap.2(sqrt(sqrt(t(ad))),scale="none",margins=c(20,15),col=brewer.pal(11,"RdBu"),trace="none",main="sqrt sqrt sample %")


