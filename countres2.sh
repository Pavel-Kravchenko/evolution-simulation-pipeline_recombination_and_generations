echo "/Bash script run successfully/"
mask=$1
home=`pwd`

cd ./Species_files

a=*.fasta.aligned.fasta
for x in $a
do
 echo "/Removing directory/"
 rm -R ${x%.fasta.aligned.fasta}"_"$mask
 echo "/Python script run/" 
 python "$home/alignment_reader_3.4.py" $x ${x%.fasta.aligned.fasta} $mask
 echo "/Done/"
 cd ${x%.fasta.aligned.fasta}"_"$mask
 pwd
 echo "/R script run/"
 Rscript "$home/R_plot.R" --no-save --no-restore --args `pwd` $mask ${x%.fasta.aligned.fasta}
 echo "/Done/"
 echo "/Cleaning the derictory/"
 mv ../${x%.fasta.aligned.fasta}.dnd ./
 mv ../$x ./
 echo "/Done/"
 echo "/Drawing trees/"
 fdrawtree -intreefile ${x%.fasta.aligned.fasta}.dnd -plotfile ${x%.fasta.aligned.fasta}"_"$mask"_"$wind.fdrawtree -auto
 fdrawgram -intreefile ${x%.fasta.aligned.fasta}.dnd -plotfile ${x%.fasta.aligned.fasta}"_"$mask"_"$wind.c.fdrawgram -auto -style c -previewer n
 cd ..
 echo "/Done/"
done
echo "/Bash script end successfully/"

