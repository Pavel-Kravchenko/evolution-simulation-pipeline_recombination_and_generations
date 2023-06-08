echo "/Bash script run successfully/"
mask=0.25
home=`pwd`
cd ./Species_files

a=*.fasta*
for x in $a
do
 cd ${x%.fasta}"_"$mask
 echo "/R script run/"
 Rscript "$home/R_plot.R" --no-save --no-restore --args `pwd` $mask ${x%.fasta}
 echo "/Done/"
 cd ..
done
echo "/Bash script end successfully/"


