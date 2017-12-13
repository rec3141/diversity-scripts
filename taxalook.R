## Rscript for network analysis based on OTUs

# would be cool to see seasonal cycle as an animation...
# static locations for members, bubble size changes with proportion over time
# connections change color for positive/negative correlations with time, width changes with strength
# you should be able to see networks form and break over time
# could do this for a transect? for the 24 h experiment using gene expression?
# high-resolution depth profile!

library(gplots) #heatmap.2
library(ggplot2)
library(corrplot)
library(MASS)
library(data.table) #fread
library(dendextend) #cutree
library(spatstat) #im

## network analysis

#install_github("zdk123/SpiecEasi")
#biocLite('phyloseq')
library(SpiecEasi)
library(Matrix)
library(phyloseq)
library(igraph)

#sharedfile = "allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.pick.shared" #this one has fewer, not sure what the diff is

bac_sharedfile = "allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared"
bac_taxfile = "allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy"

euk_sharedfile = "allsamples18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.shared"
euk_taxfile = "allsamples18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.0.02.cons.taxonomy"

bac.shared.in <- fread(input=bac_sharedfile,sep="\t",header=T, data.table=F)
euk.shared.in <- fread(input=euk_sharedfile,sep="\t",header=T, data.table=F)

bac.tax.in <- fread(input=bac_taxfile,sep=";",header=F,skip=1, data.table=F)
euk.tax.in <- fread(input=euk_taxfile,sep=";",header=F,skip=1, data.table=F)

bactmp <- as.data.frame(matrix(unlist(strsplit(as.character(bac.tax.in[,1]), split="\t",fixed=T)), ncol=3, byrow=T))
euktmp <- as.data.frame(matrix(unlist(strsplit(as.character(euk.tax.in[,1]), split="\t",fixed=T)), ncol=3, byrow=T))

bac.tax.in <- cbind(bactmp,bac.tax.in[,-1])
euk.tax.in <- cbind(euktmp,euk.tax.in[,-1])

bac.shared.meta <- bac.shared.in[,1:3]
euk.shared.meta <- euk.shared.in[,1:3]

bac.shared.abund <- bac.shared.in[,4:ncol(bac.shared.in)]
euk.shared.abund <- euk.shared.in[,4:ncol(euk.shared.in)]

bac.shared.abund <- t(bac.shared.abund)
euk.shared.abund <- t(euk.shared.abund)

colnames(bac.shared.abund) <- bac.shared.meta[,2]
colnames(euk.shared.abund) <- euk.shared.meta[,2]

bac.tax.meta <- bac.tax.in[,1:2]
euk.tax.meta <- euk.tax.in[,1:2]

bac.tax.tax <- bac.tax.in[,3:ncol(bac.tax.in)]
euk.tax.tax <- euk.tax.in[,3:ncol(euk.tax.in)]

rownames(bac.tax.tax) <- paste0("bac_",bac.tax.meta[,1])
rownames(euk.tax.tax) <- paste0("euk_",euk.tax.meta[,1])

all.tax.tax <- rbind(bac.tax.tax,euk.tax.tax)

rownames(bac.shared.abund) <- paste0("bac_",rownames(bac.shared.abund))
rownames(euk.shared.abund) <- paste0("euk_",rownames(euk.shared.abund))

#now taxa and abundances are in rows named by OTU

# need to do some scaling here to separate chl/mito in bacteria, and meta/protist in euks
chl.rows <- grep("Chloroplast",all.tax.tax[rownames(bac.shared.abund),22])
mit.rows <- grep("Mitochondria",all.tax.tax[rownames(bac.shared.abund),22])
prok.rows <- setdiff(1:nrow(bac.shared.abund),c(chl.rows,mit.rows))

met.rows <- grep("Metazoa",all.tax.tax[rownames(euk.shared.abund),5])
prot.rows <- setdiff(1:nrow(euk.shared.abund),met.rows)

org.shared.abund <- bac.shared.abund[c(chl.rows,mit.rows),]
prok.shared.abund <- bac.shared.abund[prok.rows,]

met.shared.abund <- euk.shared.abund[met.rows,]
prot.shared.abund <- euk.shared.abund[prot.rows,]

#select samples with more than minseqs bact and protist seqs
#cut samples/columns that don't have more than 1000 seqs
minseqs=1000

prok.tmp.abund <- prok.shared.abund[,which(colSums(prok.shared.abund) > minseqs)]
prot.tmp.abund <- prot.shared.abund[,which(colSums(prot.shared.abund) > minseqs)]

#select for SEWARD LINE
prok.SL <- grep("SL",colnames(prok.tmp.abund))
prot.SL <- grep("SL",colnames(prot.tmp.abund)) 

prok.tmp.abund <- prok.tmp.abund[,prok.SL]
prot.tmp.abund <- prot.tmp.abund[,prot.SL]

#select samples present in both reduced sets
#on this run some samples had multiple 16S but only one or zero 18S
#here we copy the 18S if it's present but delete the pair if it's not

prok.samples <- sub("bac_","",colnames(prok.tmp.abund))
prok.samples <- sub("_rep1","",prok.samples)
prok.samples <- sub("_rep2","",prok.samples)
prok.samples <- sub("_rep3","",prok.samples)

prot.samples <- sub("euk_","",colnames(prot.tmp.abund))

bac.in.euk <- match(prok.samples,prot.samples)

prok.small.abund <- prok.tmp.abund[,!is.na(bac.in.euk)]
prot.small.abund <- prot.tmp.abund[,bac.in.euk]
prot.small.abund <- prot.small.abund[,!is.na(bac.in.euk)]

org.small.abund <- org.shared.abund[, colnames(prok.small.abund)]
met.small.abund <- met.shared.abund[, colnames(prot.small.abund)]


#select OTUs with more than minprop proportion averaged over all samples
minprop = 0.1/minseqs #

#make table of proportions/relative abundance
prok.small.rel <- prop.table(as.matrix(prok.small.abund), margin=2)
prot.small.rel <- prop.table(as.matrix(prot.small.abund), margin=2)
org.small.rel <- prop.table(as.matrix(org.small.abund), margin=2)
met.small.rel <- prop.table(as.matrix(met.small.abund), margin=2)

# normalize counts to smallest remaining prok or prot sample
minsize <- min(c(colSums(prok.small.abund),colSums(prot.small.abund)))
prok.small.new <- round(prok.small.rel * minsize)
prot.small.new <- round(prot.small.rel * minsize)
org.small.new <- round(org.small.rel * minsize * colSums(org.small.abund) / colSums(prok.small.abund) )
met.small.new <- round(met.small.rel * minsize * colSums(met.small.abund) / colSums(prot.small.abund) )
met.small.new[is.na(met.small.new)] <- 0

#keeps rows with at least minprop in any single sample
prok.small.resample <- prok.small.new[which(apply(prok.small.new,1,function(x) any(x>0))),]
prot.small.resample <- prot.small.new[which(apply(prot.small.new,1,function(x) any(x>0))),]
org.small.resample <- org.small.new[which(apply(org.small.new,1,function(x) any(x>0))),]
met.small.resample <- met.small.new[which(apply(met.small.new,1,function(x) any(x>0))),]

prok.otus <- sub("bac_","prok_",rownames(prok.small.resample))
prot.otus <- sub("euk_","prot_",rownames(prot.small.resample))
org.otus <- sub("bac_","org_",rownames(org.small.resample))
met.otus <- sub("euk_","met_",rownames(met.small.resample))


#group the bac and euk together and calculate combined proportions and correlations

all.small.combined <- rbind(prok.small.resample, org.small.resample, prot.small.resample, met.small.resample)
all.small.otus <- rownames(all.small.combined)
all.small.taxanames <- all.tax.tax[rownames(all.small.combined),22]
rownames(all.small.combined) <- paste(c(prok.otus, org.otus, prot.otus, met.otus), all.small.taxanames)

#singleton rows are all singletons
max(all.small.combined[rowSums(all.small.combined)<2,]) # 1

#throw out singleton rows
all.small.resample <- all.small.combined[rowSums(all.small.combined)>0,]

#table of relative proportions
all.small.rel <- prop.table(as.matrix(all.small.resample),margin=2)

#manhattan transformation of relative abundance
all.small.manh <- sqrt(sqrt(all.small.rel))

#correlation matrices of transformed 
all.small.manh.colcor <- cor(all.small.manh)
all.small.manh.rowcor <- cor(t(all.small.manh))

#corrected correlations /stackoverflow
all.small.manh.colnorm <- sqrt(1 - all.small.manh.colcor^2)
all.small.manh.rownorm <- sqrt(1 - all.small.manh.rowcor^2)

#as.dist avoids another transformation, 
pdf(file="heat-samples.pdf",width=36,height=36)
hm1a <- heatmap.2(all.small.manh.colnorm, dist=as.dist, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

#distance calculation on normalized correlations, branching is the same but nicer branch lengths
pdf(file="heat-samples-dist.pdf",width=36,height=36)
hm1 <- heatmap.2(all.small.manh.colnorm, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

pdf(file="heat-otus.pdf",width=96,height=96)
hm2 <- heatmap.2(all.small.manh.rownorm,dist=as.dist, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()



# ND,D1,D2 all cluster together -- so let's just use ND (non-diluted)
all.small.nd.rel <- all.small.rel[,grep("ND",colnames(all.small.rel))]
all.small.nd <- all.small.nd.rel[rowSums(all.small.nd.rel)>0,]

#simpler names
newnames <- sub("_1","1",colnames(all.small.nd))
newnames <- sub("_2","2",newnames)
newnames <- data.frame(do.call(rbind,strsplit(newnames,"_")))
newnames <- paste(newnames[,4],newnames[,5])

colnames(all.small.nd) <- newnames

#manhattan transformation of relative abundance
all.small.nd.manh <- sqrt(sqrt(all.small.nd))

#correlation matrices of transformed 
all.small.nd.manh.colcor <- cor(all.small.nd.manh)
all.small.nd.manh.rowcor <- cor(t(all.small.nd.manh))

#corrected correlations /stackoverflow
all.small.nd.manh.colnorm <- sqrt(1 - all.small.nd.manh.colcor^2)
all.small.nd.manh.rownorm <- sqrt(1 - all.small.nd.manh.rowcor^2)


#as.dist avoids another transformation, 
pdf(file="heat-samples-nd.pdf",width=36,height=36)
hm1 <- heatmap.2(all.small.nd.manh.colnorm, dist=as.dist, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

#distance calculation on normalized correlations, branching is the same but nicer branch lengths
pdf(file="heat-samples-nd-dist.pdf",width=36,height=36)
hm1d <- heatmap.2(all.small.nd.manh.colnorm, margin=c(40,40),trace="none",col=bluered,breaks=100,cexRow=4,cexCol=4)
colcut <- cutree(hm1d$colDendrogram, k=4)
colcolors <- topo.colors(4)[colcut]
hm1d <- heatmap.2(all.small.nd.manh.colnorm, margin=c(40,40),trace="none",col=bluered,breaks=100,cexRow=4,cexCol=4,ColSideColors=colcolors)
dev.off()

#still too many rows to plot, let's use just those that have at least 0.5% in any sample
all.small.one <- all.small.nd[which(apply(all.small.nd,1,function(x) any(x>0.01))),]

#manhattan transformation of relative abundance
all.small.one.manh <- sqrt(sqrt(all.small.one))

#correlation matrices of transformed 
all.small.one.manh.colcor <- cor(all.small.one.manh)
all.small.one.manh.rowcor <- cor(t(all.small.one.manh))

#corrected correlations /stackoverflow
all.small.one.manh.colnorm <- sqrt(1 - all.small.one.manh.colcor^2)
all.small.one.manh.rownorm <- sqrt(1 - all.small.one.manh.rowcor^2)

#now for rows/OTUs
pdf(file="heat-otus-one.pdf",width=96,height=96)
hm2 <- heatmap.2(all.small.one.manh.rownorm,dist=as.dist, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

#now with extra distance calc
pdf(file="heat-otus-one-dist.pdf",width=96,height=96)
hm2d <- heatmap.2(all.small.one.manh.rownorm, margin=c(20,20),trace="none",col=bluered,breaks=100)
rowcut <- cutree(hm2d$colDendrogram, k=8)
rowcolors <- topo.colors(8)[rowcut]
hm2d <- heatmap.2(all.small.one.manh.rownorm, margin=c(20,20),trace="none",col=bluered,breaks=100,RowSideColors=rowcolors)
dev.off()


# heat map with relative abundance, dendrogram of samples and taxa from correlation meta-analysis
# use the same dendrogram but plot abundance
pdf(file="heat-relabund.pdf",width=96,height=96)
lmat = rbind(c(1,2),c(3,4))
lwid = c(25,1)
lhei = c(25,1)
#hm4 <- heatmap.2( log10((all.small.one*1000+1)), Rowv=hm2d$colDendrogram, Colv=hm1d$colDendrogram, margin=c(10,10),trace="none",col=bluered,breaks=100,scale="none", lmat=lmat, lwid = lwid, lhei = lhei, RowSideColors=rowcolors, ColSideColors=colcolors)
hm4 <- heatmap.2( sqrt(sqrt(all.small.one)), Rowv=hm2d$colDendrogram, Colv=hm1d$colDendrogram, margin=c(60,100),trace="none",col=bluered,breaks=100,scale="none", RowSideColors=rowcolors, ColSideColors=colcolors, cexRow=3, cexCol=8)
#dev.off()

# heat map with relative abundance, dendrogram of samples and taxa from correlation meta-analysis
# use the same dendrogram but plot abundance

#first heat map, use to get dendrogram colors / could just use hclust
hm5a <- heatmap.2( log10(all.small.one*1000000+1), margin=c(10,10),trace="none",col=bluered,breaks=100,scale="none")
dev.off()

colcut <- cutree(hm5a$colDendrogram, k=3, order_clusters_as_data=F)
#colcut <- cutree(hm5a$colDendrogram, h=round(max(get_branches_heights(hm5$colDendrogram))/1.5))

#bothcolors <- colorpanel(3,low="white",mid="lightpink",high="skyblue1")
#collist <- c("#FFFFFF","#FFB6C1","#87CEFF")
collist <- c("white","lightpink","lightblue")
colcolors <- as.vector(rbind(collist,rev(collist)))[colcut]

c("rare"="grey","rare deep"="lightblue1","rare surface"="lightpink1","endemic deep"="lightblue2","endemic surface"="lightpink2","cosmopolitan deep"="lightblue3","cosmopolitan surface"="lightpink3")

rowcut <- cutree(hm5a$rowDendrogram, k=7, order_clusters_as_data=F)
#rowcut <- cutree(hm5a$rowDendrogram, h=round(max(get_branches_heights(hm5$rowDendrogram))/1.8))
#rowcolors <- colorpanel(max(rowcut),low="sienna1",mid="white",high="slateblue")
#rowcolors <- as.vector(rbind(rowcolors,rev(rowcolors)))[rowcut]
#rowlist <- c("#FFFFFF","#87CEFF","#FFB6C1","#87CEFF","#FFB6C1","#87CEFF","#FFB6C1")
rowlist <- c("grey","lightblue1","lightpink1","lightblue2","lightpink2","lightblue3","lightpink3")
rowcolors <- rowlist[rowcut]

all.small.oneb <- log10(all.small.one[names(rowcut),names(colcut)]*1000000+1)

#second heat map, use to get dendrograms / could just use hclust
hm5b <- heatmap.2( all.small.oneb, margin=c(10,10),trace="none",col=bluered,breaks=100,scale="none")
dev.off()

#some hijinks to get columns in order
colweights <- colcut
colweights[colweights==1] <- 2.5
coldend <- reorder(hm5b$colDendrogram,colweights,mean)

#final heat map to print
pdf(file="heat-relabund.pdf",width=96,height=96)
hm5d <- heatmap.2( all.small.oneb, Rowv=hm5b$rowDendrogram, Colv=coldend, trace="none", breaks=20,scale="none", RowSideColors=rowcolors, ColSideColors=colcolors, cexRow=3, cexCol=8, key.par=list(cex=6), key.title="relative abundance (per million)", key.xlab = "log10", key.ylab = NA, keysize=1, denscol="black", 
	col=colorpanel(19,low="gray66",mid="white",high="blue"),
	lmat=rbind(
		c(6,0,0,0,0),
		c(0,0,0,2,0),
		c(0,0,0,5,0),
		c(0,1,4,3,1),
		c(0,0,0,2,0)),
    lhei=c(3,.1,1,9,1),
    lwid=c(3,.1,1,9,1.5),
    offsetRow=5, offsetCol=5, margin=c(0,0))
dev.off()

# Figure 1 heatmap of relabundance of samples with OTUs
#DONE

break

#MDS

# Figure 2 ODV plot of relative abundances of each of the major clades of heatmap

name.rows <- c("rare","rare deep","rare surface","endemic deep","endemic surface","cosmopolitan deep","cosmopolitan surface")
name.cols <- c("mixed","deep","surface")

all.small.agg.row <- aggregate(all.small.oneb, by=list(rowcut), FUN=mean)
all.small.agg.row <- t(all.small.agg.row[,2:ncol(all.small.agg.row)])
colnames(all.small.agg.row) <- name.rows

all.small.agg.col <- aggregate(t(all.small.oneb), by=list(colcut), FUN=mean)
all.small.agg.col <- t(all.small.agg.col[,2:ncol(all.small.agg.col)])
colnames(all.small.agg.col) <- name.cols

all.small.agg.both <- aggregate(all.small.agg.row, by=list(colcut), FUN=mean)
all.small.agg.both <- t(all.small.agg.both[,2:ncol(all.small.agg.both)])
colnames(all.small.agg.both) <- name.cols


# OTU MDS
all.small.mds.otus <- as.data.frame(cmdscale(dist(all.small.oneb)))
colnames(all.small.mds.otus) <- c("x","y")
all.small.mds.otus$size <- rowMeans(all.small.agg.col)
all.small.mds.otus$color <- rep.int(x=name.rows,times=unname(rle(rowcut)$lengths))
#qplot(x,y,data=all.small.mds.otus, color=color, size=size)

g <- ggplot(all.small.mds.otus, aes(x, y, color=color, size=size)) + 
	geom_point(
		aes(fill=color, size=size), 
		colour="black", shape=21, stroke = .5) +
	labs(x="MDS1",y="MDS2") +
	theme(
		legend.title=element_blank(), 
		legend.key=element_rect(fill=NA), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = 'grey89'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
    ) +
    guides(fill="legend",colour = guide_legend(override.aes = list(size=4))) +
    scale_fill_manual(guide="legend",name='', values=c("rare"="grey","rare deep"="lightblue1","rare surface"="lightpink1","endemic deep"="lightblue2","endemic surface"="lightpink2","cosmopolitan deep"="lightblue3","cosmopolitan surface"="lightpink3"))
g
ggsave(file="mds-otus.pdf", width=8, height=6)

# MDS Samples
all.small.mds.sample <- as.data.frame(cmdscale(dist(t(all.small.oneb))))
colnames(all.small.mds.sample) <- c("x","y")
all.small.mds.sample$size <- rowMeans(all.small.agg.row)
all.small.mds.sample$color <- rep.int(x=name.cols,times=unname(rle(colcut)$lengths))
#qplot(x,y,data=all.small.mds.sample, color=color, size=size)

h <- ggplot(all.small.mds.sample, aes(x, y, color=color, size=size)) + 
	geom_point(
		aes(fill=color, size=size), 
		colour="black", shape=21, stroke = .5) +
	labs(x="MDS1",y="MDS2") +
	theme(
		legend.title=element_blank(), 
		legend.key=element_rect(fill=NA), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = 'grey89'),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank()
    ) +
    guides(fill="legend",colour = guide_legend(override.aes = list(size=4))) + scale_fill_manual(guide="legend",name='', values=c("mixed"="white","deep"="lightblue","surface"="lightpink"))
    h
ggsave(file="mds-samples.pdf", width=8, height=6)


cbind(rle(all.small.mds.otus$color)$values,rle(all.small.mds.otus$color)$lengths)
cbind(rle(all.small.mds.sample$color)$values,rle(all.small.mds.sample$color)$lengths)


## ODV export
all.small.export <- rbind(all.small.one,t(all.small.mds.sample[colnames(all.small.one),]),t(10^all.small.agg.row[colnames(all.small.one),]))
all.small.export <- cbind(all.small.export,rbind(all.small.mds.otus[rownames(all.small.one),],0,0,0,0,0,0,0,0,0,0,0))
all.small.export <- t(all.small.export)
write.table(file="all.small.export.txt",all.small.export,quote=F,sep="\t")




# Figure 3 network diagram showing major OTUs, colors=taxonomy, shape=heatmap clade

#NETWORK DIAGRAMS

all.small.se <- round(t(all.small.one*1000000)) #matrix of normalized counts
all.small.se.taxanames <- sub("^met_","euk_",colnames(all.small.se))
all.small.se.taxanames <- sub("^prot_","euk_",all.small.se.taxanames)
all.small.se.taxanames <- sub("^prok_","bac_",all.small.se.taxanames)
all.small.se.taxanames <- sub("^org_","bac_",all.small.se.taxanames)

all.small.se.otus <- do.call(rbind,strsplit(all.small.se.taxanames," "))[,1]

#make nice colors for our blobs
domain.levels <- all.tax.tax[all.small.se.otus,2] #Bacteria, Archaea
kingdom.levels <- all.tax.tax[all.small.se.otus,5] #Metazoa
class.levels <- all.tax.tax[all.small.se.otus,13] #Chloroplast, Alpha, Beta, Gamma, Flavo, Syndiniales
family.levels <- all.tax.tax[all.small.se.otus,20] #Mitochondria

vcolor <- rep("#427df4",times=ncol(all.small.se)) #blue default Bacteria
vcolor[grep("Archaea",domain.levels)] <- "#f9f03eCC" #yellow
vcolor[grep("Eukaryota",domain.levels)] <- "#f44265" #red
vcolor[grep("Metazoa",kingdom.levels)] <- "#f463f9CC" #purple
vcolor[grep("Chloroplast",class.levels)] <- "#0fe047CC" #green
vcolor[grep("Alphaproteobacteria",class.levels)] <- "#0033ff" #blue
vcolor[grep("Betaproteobacteria",class.levels)] <- "#1947ff" #blue
vcolor[grep("Gammaproteobacteria",class.levels)] <- "#1dcdf9" #blue
vcolor[grep("Flavobacteriia",class.levels)] <- "#0fe047CC" #blue
vcolor[grep("Syndiniales",class.levels)] <- "#ba1bf9" #purple
vcolor[grep("Mitochondria",family.levels)] <- "#42e5f4" #purple


#### run sparCC algorithm for correlation analysis

all.small.sparcc <- sparcc(all.small.se,iter=100,inner_iter=20,th=0.1) #normal
#all.small.sparcc <- sparcc(all.small.se,iter=1,inner_iter=1,th=0.1) #quick

all.small.sparcc.cor <- all.small.sparcc$Cor
rownames(all.small.sparcc.cor) <- colnames(all.small.se)
colnames(all.small.sparcc.cor) <- colnames(all.small.se)

#heatmap, use correlation matrix directly
pdf(file="heat-sparcc.pdf",width=36,height=36)
hm0 <- heatmap.2(all.small.sparcc.cor,distfun=function(c) as.dist(1 - c), margin=c(20,20),trace="none",col=bluered,breaks=20,cexRow=0.5, cexCol=0.5)
dev.off()


#changing the first block will affect the layout
#changing the second and third will affect which edges are shown

## Define arbitrary threshold for SparCC correlation matrix for the graph

#	lowcutoff <- rev(sort(abs(all.small.sparcc$Cor)))[round(length(all.small.sparcc$Cor)/2/5)] #arbitrary cutoff, i = 10 = uses top 10% correlations
	lowcutoff <- 0.8
	sparcc.graph <- abs(all.small.sparcc.cor) >= lowcutoff
	diag(sparcc.graph) <- 0
	sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)
	ig.sparcc <- adj2igraph(sparcc.graph)
	to.del <- which(degree(ig.sparcc)<1)
	to.keep <- which(degree(ig.sparcc)>0)
	ig.sparcc <- delete.vertices(ig.sparcc,to.del)
	
	#try a few different layouts
	fr.coord <- layout.fruchterman.reingold(ig.sparcc) #ok, shows separation
	mds.coord <- layout_with_mds(ig.sparcc) #olverlapping
	drl.coord <- layout_with_drl(ig.sparcc) #ok
	kk.coord <- layout_with_kk(ig.sparcc) #spread out
	lgl.coord <- layout_with_lgl(ig.sparcc) #spread out
	use.coord <- mds.coord

	xlim <- range(use.coord[,1])
	ylim <- range(use.coord[,2])
	vsize <- sqrt(colSums(all.small.se))
	vsize <- vsize/max(vsize)*100


#for(i in 100) {
	highcutoff <- round(lowcutoff,2)
	highcutoff <- .80

#	sparcc.graph.pos <- (all.small.sparcc$Cor >  (i-1)/100 ) & (all.small.sparcc$Cor <= i/100)	
	sparcc.graph.pos <- all.small.sparcc.cor >= highcutoff
	diag(sparcc.graph.pos) <- 0
	sparcc.graph.pos <- Matrix(sparcc.graph.pos, sparse=TRUE)
	pos.weights <- (sparcc.graph.pos & upper.tri(sparcc.graph.pos))*all.small.sparcc.cor
	ig.sparcc.pos <- adj2igraph(sparcc.graph.pos)
	E(ig.sparcc.pos)$weight <- pos.weights[pos.weights>0]^4
	E(ig.sparcc.pos)$weight <- E(ig.sparcc.pos)$weight/max(E(ig.sparcc.pos)$weight)
	E(ig.sparcc.pos)$color = rgb(0,0,1,E(ig.sparcc.pos)$weight)
	ig.sparcc.pos <- delete.vertices(ig.sparcc.pos,to.del)

#	sparcc.graph.neg <- (all.small.sparcc$Cor <= -1 * (i-1)/100 ) & (all.small.sparcc$Cor > -1 * i/100)	
	sparcc.graph.neg <- all.small.sparcc.cor <= -1 * highcutoff
	diag(sparcc.graph.neg) <- 0
	sparcc.graph.neg <- Matrix(sparcc.graph.neg, sparse=TRUE)
	neg.weights <- (sparcc.graph.neg & upper.tri(sparcc.graph.neg))*all.small.sparcc.cor
	ig.sparcc.neg <- adj2igraph(sparcc.graph.neg)
	E(ig.sparcc.neg)$weight <- abs(neg.weights[neg.weights<0])^4
	E(ig.sparcc.neg)$weight <- E(ig.sparcc.neg)$weight/max(E(ig.sparcc.neg)$weight)
	E(ig.sparcc.neg)$color = rgb(1,0,0,E(ig.sparcc.neg)$weight)
	ig.sparcc.neg <- delete.vertices(ig.sparcc.neg,to.del)


	pdf(file=paste0("graph-sparcc-",highcutoff,".pdf"),width=96,height=96)
	plot(ig.sparcc, layout=use.coord, xlim=xlim, ylim=ylim, vertex.shape="none", vertex.label=NA, edge.label=NA, edge.color="#a3a3a3CC", edge.width=0.1, main="sparcc", rescale=F)
	plot(ig.sparcc.pos, layout=use.coord, xlim=xlim, ylim=ylim, vertex.shape="none", vertex.label=NA, edge.label=NA, edge.width=10*E(ig.sparcc.pos)$weight, main="sparcc", rescale=F, add=T)
	plot(ig.sparcc.neg, layout=use.coord, xlim=xlim, ylim=ylim, vertex.shape="none", vertex.label=NA, edge.label=NA, edge.width=10*E(ig.sparcc.neg)$weight, main="sparcc", rescale=F, add=T)
	plot(ig.sparcc, layout=use.coord, xlim=xlim, ylim=ylim, vertex.shape="circle", vertex.size=vsize[to.keep], vertex.color=vcolor[to.keep], vertex.label.cex=4, vertex.label.color="black", 
		vertex.label=as.character(all.small.mds.otus[colnames(all.small.se),"color"])[to.keep], 
		all.small.mds.otus$color,		
		edge.lty="blank", main="sparcc", rescale=F, add=T)
	dev.off()
#		vertex.label=as.character(colnames(all.small.se))[to.keep], 

# }


#as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22]))

#### run SpiecEasi mb algorithm for correlation analysis
se.mb.small <- spiec.easi(all.small.se, method='mb', lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=50,ncores=4))
#took hours on 300x1100

#make graph
ig.mb <- adj2igraph(se.mb.small$refit)
fr.coord <- layout.fruchterman.reingold(ig.mb) #ok
mds.coord <- layout_with_mds(ig.mb) #overlapping
drl.coord <- layout_with_drl(ig.mb) #pretty good
kk.coord <- layout_with_kk(ig.mb) #spread out
lgl.coord <- layout_with_lgl(ig.mb) #spread out
use.coord <- lgl.coord

xlim <- range(use.coord[,1])
ylim <- range(use.coord[,2])
vsize <- sqrt(colSums(all.small.se))
vsize <- vsize/max(vsize)*50

pdf(file="graph-mb.pdf",width=96,height=96)
#plot(ig.mb, layout=use.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.mb, layout=use.coord, xlim=xlim, ylim=ylim, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=3, vertex.label.color="black", vertex.label=paste(as.character(colnames(all.small.se)),as.character(all.small.mds.otus[colnames(all.small.se),"color"])), edge.color="#a3a3a3CC", main="sparcc", rescale=F)
dev.off()
#as.character(colnames(all.small.se))



#### run SpiecEasi glasso algorithm for correlation analysis
se.gl.small <- spiec.easi(all.small.se, method='glasso', lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=50,ncores=4))
#took overnight on 300x1100

#make graph
#add transposed graphs to get symmetrical graph due to inconsistent algorithm
se.gl.small.refit <- tril(se.gl.small$refit) + t(tril(se.gl.small$refit))

#se.gl.small.icov <- tril(Matrix(se.gl.small$icov)) + t(tril(se.gl.small$icov)) #20 matrices

ig.gl <- adj2igraph(se.gl.small.refit)
fr.coord <- layout.fruchterman.reingold(ig.gl) #has structure
mds.coord <- layout_with_mds(ig.gl) #ok
drl.coord <- layout_with_drl(ig.gl) #strong separation
kk.coord <- layout_with_kk(ig.gl) #some structure
lgl.coord <- layout_with_lgl(ig.gl) #nothing
circle.coord <- layout_in_circle(ig.gl)
use.coord <- fr.coord

xlim <- range(use.coord[,1])
ylim <- range(use.coord[,2])
vsize <- sqrt(colSums(all.small.se))
vsize <- vsize/max(vsize)*10


pdf(file="graph-glasso.pdf",width=96,height=96)
plot(ig.gl, layout=use.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.9, vertex.label.color="black", vertex.label=colnames(all.small.se), main="glasso")
dev.off()

#as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22]))
#plot(ig.gl, layout=use.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
#par(mfrow=c(2,3))
#plot(ig.gl, layout=fr.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="fr")
#plot(ig.gl, layout=mds.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="mds")
#plot(ig.gl, layout=kk.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="kk")
#plot(ig.gl, layout=lgl.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="lgl")
#plot(ig.gl, layout=circle.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="circle")


#glasso correlation matrix

se.gl.small.cor <- triu(cov2cor(se.gl.small$opt.cov)*se.gl.small$refit, k=1)
se.gl.small.cor <- se.gl.small.cor + t(se.gl.small.cor)
ig.gl.cor <- adj2igraph(se.gl.small.cor)

drl.coord <- layout_with_drl(ig.gl.cor)
#Error in .Call("R_igraph_layout_drl", graph, seed, use.seed, options,  : 
#  At DensityGrid.cpp:230 : Exceeded density grid in DrL, Internal DrL error
use.coord <- drl.coord

xlim <- range(use.coord[,1])
ylim <- range(use.coord[,2])
vsize <- sqrt(colSums(all.small.se))
vsize <- vsize/max(vsize)*10
pdf(file="tmp-graph-gl.pdf",width=96,height=96)
plot(ig.gl.cor, layout=drl.coord, vertex.shape="circle", vertex.size=vsize, vertex.color=vcolor, vertex.label.cex=0.3, vertex.label.color="black", vertex.label=as.character(paste(colnames(all.small.se),all.tax.tax[colnames(all.small.se),22])), main="drl")
dev.off()


#from tutorial

elist.gl <- summary(triu(cov2cor(se.gl.small$opt.cov)*se.gl.small$refit, k=1))
elist.mb <- summary(symBeta(getOptBeta(se.mb.small), mode='maxabs'))
elist.sparcc <- summary(sparcc.graph*all.small.sparcc$Cor)

hist(elist.sparcc[,3], main="", xlab="edge weights", ylim=c(0,3000))
hist(elist.mb[,3], add=TRUE, col='forestgreen', ylim=c(0,3000))
hist(elist.gl[,3], add=TRUE, col='red', ylim=c(0,1000))





## more plotting

pdf(file="heat-samples-new.pdf",width=36,height=36)
hm1 <- heatmap.2(all.small.colcor,margin=c(20,20),trace="none",col=bluered,breaks=100)
#corrplot(cor(small.manh))
dev.off()

pdf(file="heat-otus-new.pdf",width=96,height=96)
hm2 <- heatmap.2(all.small.rowcor,margin=c(20,20),trace="none",col=bluered,breaks=100)
#corrplot(cor(small.manh))
dev.off()


#use the same dendrogram but label with sample containing highest relative abundance
#plot correlation matrix by label
#> dim(all.small.se)
#[1]  268 1122
#replace colnames (OTUs) with sample name with most abundance
all.replot <- all.small.resample
rownames(all.replot) <- paste(rownames(all.replot),colnames(all.replot)[unname(apply(all.replot,1,which.max))])

all.replot.rel <- prop.table(as.matrix(all.replot),margin=2)
all.replot.manh <- sqrt(sqrt(all.replot.rel))

all.replot.colcor <- cor(all.replot.manh)
all.replot.rowcor <- cor(t(all.replot.manh))

all.replot.colnorm <- sqrt(1 - all.replot.colcor^2)
all.replot.rownorm <- sqrt(1 - all.replot.rowcor^2) #stackoverflow

for(i in 1:nrow(all.replot.rowcor)) {for(j in 1:nrow(all.replot.rowcor)) {plot(all.replot.rowcor[i,],all.replot.rowcor[j,],xlab=rownames(all.replot.rowcor)[i],ylab=rownames(all.replot.rowcor)[j])}}


pdf(file="heat-samples-replot.pdf",width=48,height=48)
hm3 <- heatmap.2(all.replot.colnorm,dist=as.dist,margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

pdf(file="heat-otus-replot.pdf",width=96,height=96)
hm4 <- heatmap.2(all.replot.rownorm,dist=as.dist, margin=c(20,20),trace="none",col=bluered,breaks=100)
dev.off()

rownames(all.small.resample) <- paste(rownames(all.small.resample), all.tax.tax[rownames(all.small.resample),22])

# heat map with relative abundance, dendrogram of samples and taxa from correlation meta-analysis
# use the same dendrogram but plot abundance
pdf(file="heat-relabund.pdf",width=96,height=96)
lmat = rbind(c(1,2),c(3,4))
lwid = c(25,1)
lhei = c(25,1)
hm4 <- heatmap.2( log10((all.small.resample+1)), Rowv=hm2$colDendrogram, Colv=hm1$colDendrogram, margin=c(10,10),trace="none",col=bluered,breaks=100,scale="none", lmat=lmat, lwid = lwid, lhei = lhei)
dev.off()


#get pairs of taxa for mito, chloropast, syndiniales...
coln <- nrow(all.replot.colcor)
colr <- nrow(all.replot.rowcor)

rowcut <- cutree(hm2$colDendrogram, k=max(get_branches_heights(hm2$colDendrogram))/2)
colcut <- cutree(hm1$colDendrogram, k=8)

#write.table(file="pairs-rows.txt",rowcut, sep="\t")
#write.table(file="pairs-cols.txt",colcut, sep="\t")

rowphylo <- as.phylo(hm2$colDendrogram)
colphylo <- as.phylo(hm1$colDendrogram)

xchl <- grep("Chloroplast",rowphylo$tip.label)
xmit <- grep("Mitochondria",rowphylo$tip.label)
xsyn <- grep("Syndiniales",rowphylo$tip.label)

ychl <- unlist(lapply(xchl, function(x) getSisters(rowphylo, x)))
ymit <- unlist(lapply(xmit, function(x) getSisters(rowphylo, x)))
ysyn <- unlist(lapply(xsyn, function(x) getSisters(rowphylo, x)))

chl <- cbind(rowphylo$tip.label[xchl],rowphylo$tip.label[ychl])
mit <- cbind(rowphylo$tip.label[xmit],rowphylo$tip.label[ymit])
syn <- cbind(rowphylo$tip.label[xsyn],rowphylo$tip.label[ysyn])

#pdf(file="tmp-corr.pdf",width=36,height=36)
#small.colcorrplot <- corrplot(small.colcor, order="hclust",method="color",tl.cex=0.5)
#dev.off()

pdf("tree.pdf",width=12,height=48)
#plotTree(rowphylo,cex=0.01)
plot(rowphylo,cex=0.2)
dev.off()


#get.oturep(fasta=allsamples18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=allsamples18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, large=true, column=allsamples18S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, name=allsamples18S.trim.contigs.good.names)
#get.oturep(fasta=allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.opti_mcc.unique_list.list, large=true, column=allsamples16S.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, name=allsamples16S.trim.contigs.good.names)




