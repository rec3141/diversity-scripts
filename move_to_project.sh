cat sample_ids_project.csv | cut -f1-3 -d$'\t' | tr '\t' '_' > tmpa
cat sample_ids_project.csv | cut -f6 -d$'\t' > tmpb
paste tmpa tmpb > tmpc

for DIR in `cat tmpb | sort | uniq`; do
mkdir $DIR;
done

while read LINE;
do
	strarr=($LINE)
	mv ${strarr[0]}_R1.fq ${strarr[1]}
	mv ${strarr[0]}_R2.fq ${strarr[1]}
done < tmpc

mkdir empty
mv *R1.fq empty/
mv *R2.fq empty/
