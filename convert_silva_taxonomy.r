# download the SILVA taxa mapping
# wget http://www.arb-silva.de/fileadmin/silva_databases/current/Exports/taxonomy/tax_slv_ssu_123.txt

# cut -f3 tax_slv_ssu_123.txt | sort | uniq -c
# what a mess

# reordered by taxonomic level
# 	8275 genus
# 	  26 subfamily
# 	1412 family
# 	  10 superfamily
# 	  21 suborder
# 	1138 order
# 	  16 superorder
# 	   6 infraclass
# 	  58 subclass
# 	 541 class
# 	   1 superclass
# 	   1 infraphylum
# 	  33 subphylum
# 	 239 phylum
# 	   5 superphylum
# 	   1 infrakingdom
# 	   4 subkingdom
# 	  15 kingdom
# 	   2 superkingdom
# 	   4 major_clade
# 	   3 domain
#      1       # out of all the labels, they leave Escherichia blank??

# start from scratch with the silva.nr_v123.align headers
# grep '>' silva.nr_v123.align | cut -f1,3 | cut -f2 -d'>' > silva.nr_v123.full

#!/usr/bin/R

map.in <- read.table("tax_slv_ssu_123.txt",header=F,sep="\t",stringsAsFactors=F)
map.in <- map.in[,c(1,3)]
colnames(map.in) <- c("taxlabel","taxlevel")

#fix Escherichia nonsense
map.in$taxlevel[which(map.in$taxlabel=="Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales;Enterobacteriaceae;Escherichia;")] = "genus"

taxlevels <- c("root","domain","major_clade","superkingdom","kingdom","subkingdom","infrakingdom","superphylum","phylum","subphylum","infraphylum","superclass","class","subclass","infraclass","superorder","order","suborder","superfamily","family","subfamily","genus")
taxabb <- c("ro","do","mc","pk","ki","bk","ik","pp","ph","bp","ip","pc","cl","bc","ic","po","or","bo","pf","fa","bf","ge")
tax.mat <- matrix(data="",nrow=nrow(map.in),ncol=length(taxlevels))
tax.mat[,1] <- "root"
colnames(tax.mat) <- taxlevels

#change this if you want just the default 6
outlevels <- taxlevels
#outlevels <- c("domain","phylum","class","order","family","genus")

#put all the taxa into a matrix

for (i in 1:nrow(map.in)) {
	taxname <- unlist(strsplit(as.character(map.in[i,1]), split=';'))
	print(taxname);

	while ( length(taxname) > 0) {
		#regex to look for exact match

		tax.exp <- paste(paste(taxname,collapse=";"),";",sep="")
		tax.match <- match(tax.exp,map.in$taxlabel)
		tax.mat[i,map.in[tax.match,2]] <- tail(taxname,1)
		taxname <- head(taxname,-1)
	}
}


for (i in 1:nrow(tax.mat)) {
	#this fills in the empty gaps by using the closest higher taxonomic level appended with an abbreviation for the current taxonomic level
	#if you don't want this behavior, cut it out
	for (j in 1:ncol(tax.mat)) {
		if(tax.mat[i,j] < 0) { tax.mat[i,j] <- paste(tmptax,taxabb[j],sep="_")}
		else { tmptax <- tax.mat[i,j]}
	}

	#this maps the new name to the input taxonomic levels
	map.in[i,"taxout"] <- paste(paste(tax.mat[i,outlevels],collapse=";"),";",sep="")
}

# replace spaces with underscores
map.in$taxout <- gsub(" ","_",map.in$taxout)

# bring in the old taxonomic levels from SILVA and remap them using the new levels

tax.in <- read.table("silva.nr_v123.full",header=F,stringsAsFactors=F,sep="\t")
colnames(tax.in) <- c("taxid","taxlabel")
tax.in$id <- 1:nrow(tax.in)

tax.write <- merge(tax.in,map.in,all.x=T,sort=F)
tax.write <- tax.write[order(tax.write$id),]

write.table(tax.write[,c("taxid","taxout")],file="silva.nr_v123.full.tax",sep="\t",row.names=F,quote=F,col.names=F)


  
  


