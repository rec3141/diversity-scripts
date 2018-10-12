#/bin/bash

# this program provides a CSV file that can be loaded into Bandage
# containing taxonomy labels and colors for each taxa
# or used to color vizbin tsne's in R

# input is one or more genome bins as .fa files
Bin="$@"

#pathway=sketch or pathway=prokka
pathway="prokka"

for file in $Bin; do
	echo "start $file"

	f1=`basename $file .fasta`

	sed -E 's/(^>.{37}).*/\1/' $file > $f1.assembly.fasta #prokka limits to 37 chars

    if [ "$pathway" == "sketch" ]; then
    #run sendsketch
    /work/cryomics/apps/bbmap/sendsketch.sh in=$file address=nt mode=sequence records=1 out=$f1.sketch.tsv ow=t

	# match orfs to nodes
	awk -v OFS=$'\t' 'BEGIN{feature=0}{if($0 ~ /^Query/){feature=$2};if($0 ~ /^[0-9]/){print feature,$12 " " $13}}' $f1.sketch.tsv | sort -k2,2 | tr $'\t' ',' > $f1.nodes.csv

    fi;

    if [ "$pathway" == "prokka" ]; then
	#run prokka

	/work/cryomics/apps/prokka-biolinux/prokka/bin/prokka --metagenome --outdir prokka_$f1 --prefix $f1 --locustag $f1 --force --fast --rawproduct $f1.assembly.fasta

	#run kaiju
	/work/cryomics/apps/kaiju-1.6.3/bin/kaiju -t /scratch/cryomics/reference_dbs/kaiju/2017/nodes.dmp -f /scratch/cryomics/reference_dbs/kaiju/2017/kaiju_db_nr_euk.fmi -p -i prokka_$f1/$f1.faa -a greedy -e 5 -o $f1.kaiju.tsv -z 12

	#get orf taxonomy
	/work/cryomics/apps/kaiju-1.6.3/bin/addTaxonNames -t /scratch/cryomics/reference_dbs/kaiju/2017/nodes.dmp -n /scratch/cryomics/reference_dbs/kaiju/2017/names.dmp -r superkingdom,class,genus -i $f1.kaiju.tsv -o $f1.kaiju_names.tsv

  #taxonomy list
  sort -k3 $f1.kaiju_names.tsv | cut -f3,4 | sort -k1,1 -u | tr $'\t' '|' > $f1.kaiju_names_sort.tsv

	# match orfs to nodes
	awk -v OFS='\t' 'BEGIN{feature=0}{if($0 ~ /^>/){feature=$2};if($0 ~ /locus_tag/){print feature,$2}}' prokka_$f1/$f1.tbl > $f1.kaiju_nodes.tsv

	#merge orfs and nodes, keeping only most common taxonomic hit
    join -1 2 -2 2 -t $'\t' $f1.kaiju_nodes.tsv $f1.kaiju_names.tsv | cut -f2,4 | sort | uniq -c | sed 's/^ *//' | sort -k2,2 -k1,1rn | awk '{print $3 "\t" $2}' |  grep -v '^0' | uniq -f1 | sort -k1,1 | tr $'\t' '|' > $f1.nodes.tsv
    join -1 1 -2 1 -t '|' $f1.nodes.tsv $f1.kaiju_names_sort.tsv | cut -f2- -d'|'| tr '|' ',' | sort -k2 -t',' >  $f1.nodes.csv

    fi;


	#get orf colors
  cut -f2 -d',' $f1.nodes.csv | uniq > $f1.names.tsv
	while read line; do echo -e "#"`md5 <(echo -e $line) | cut -f4 -d' '| cut -b1-6`;  done < $f1.names.tsv > $f1.colors.tsv
	paste -d',' $f1.colors.tsv $f1.names.tsv > $f1.taxcolors.csv

	#merge orfs and colors
	echo -e "Node,Taxonomy,Color" > $f1.bandage.csv
  join -1 2 -2 2 -t',' $f1.nodes.csv $f1.taxcolors.csv | awk -v OFS=',' -v FS=',' '{print $2,$1,$3}' >> $f1.bandage.csv

  #merge with vizbin
  #sort -k1,1 -t',' $f1.vizbin.txt > tmp1
  #sort -k1,1 -t',' $f1.bandage.csv > tmp2
  #join -1 1 -2 1 -t',' tmp1 tmp2 > $f1.rcolor.csv

	echo "done $file"
done;

