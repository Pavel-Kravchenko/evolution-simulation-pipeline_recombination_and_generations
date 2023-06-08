echo "/Bash script srarted/"
mask=0.25
home=`pwd`
species=128

#rm Files
#mkdir Files
cd ./Files

#for generation in 4 8 16 32 64 128 256 512 
#do
# for recombination in 0 1 10 25 50 80 100
# do
#  for mutation in  1 2 4 #8 16 32
#  do
#   python "$home/simulation_script.py" $species $generation $mutation $recombination
#  done
# done
#done


a=*.fasta
for x in $a
do
 echo "/Removing directory/"
 rm -R ${x%.fasta}"_"$mask
 echo "/Python script run/" 
 python "$home/alignment_reader_script.py" $x ${x%.fasta} $mask
 echo "/Done/"
 cd ${x%.fasta}"_"$mask
 pwd
 echo "/R script run/"
 Rscript "$home/R_plot.R" --no-save --no-restore --args `pwd` $mask ${x%.fasta}
 echo "/Done/"
 #echo "/Cleaning the derictory/"
 #mv ../$x ./
 #echo "/Done/"
 cd ..
done
echo "/Python heatmap script run/"
python "$home/heatmap_df_maker.py"
echo "/Done/"
cd ..
#tar -cvf archive.tar ./Files

echo "/Bash script ended successfully/"

