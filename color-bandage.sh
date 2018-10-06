#/bin/bash

# this program provides a CSV file that can be loaded into Bandage
# containing taxonomy labels and colors for each taxa

# input is one or more .gfa files
GFA="$@"

for file in $GFA; do
	echo "start $file"

	f1=`basename $file .gfa`
	mkdir $f1
	cd $f1
	ln -s ./../$file $file
	grep '^S'  $f1.gfa | cut -f2,3 | sed 's/^\([0-9]\)/>\1/' | tr $'\t' $'\n' > $f1.assembly.fasta

    #run sendsketch
    /work/cryomics/apps/bbmap/sendsketch.sh in=$f1.assembly.fasta address=nt mode=sequence records=1 out=$f1.sketch.tsv ow=t

	# match orfs to nodes
	awk -v OFS=$'\t' 'BEGIN{feature=0}{if($0 ~ /^Query/){feature=$2};if($0 ~ /^[0-9]/){print feature,$12}}' $f1.sketch.tsv | sort -k2,2 | tr $'\t' ',' > $f1.nodes.csv

	#get orf colors
    cut -f2 -d',' $f1.nodes.csv | uniq > $f1.names.tsv
	while read line; do echo "#"`md5 <(echo -e $line) | cut -f4 -d' '| cut -b1-6`;  done < $f1.names.tsv > $f1.colors.tsv
	paste -d',' $f1.colors.tsv $f1.names.tsv > $f1.taxcolors.csv
	
	#merge orfs and colors
	echo -e "Node,Taxonomy,Color" > $f1.bandage.csv
    join -1 2 -2 2 -t',' $f1.nodes.csv $f1.taxcolors.csv | awk -v OFS=',' -v FS=',' '{print $2,$1,$3}' >> $f1.bandage.csv
    
	echo "done $file"
	cd -
done;

