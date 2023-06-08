echo "/Bash script srarted/"
mask=0.25
home=`pwd`
species=500

#rm Files
mkdir Files2
cd ./Files2

for generation in 4 8 16 32 64 128 256 512 
do
 for recombination in 0 1 5 10 15 25 50 80 100
 do
  for mutation in 0 1 2 4 #8 16 32
  do
   python "$home/simulation_script.py" $species $generation $mutation $recombination
  done
 done
done
