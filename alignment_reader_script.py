#--------------biblio & modules-------------------------------------
import sys
from sys import argv
from Bio import AlignIO 
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
from time import time
import os

tic = time()
#--------------extracting file-----------------------------

if len(sys.argv) == 1:
    print("Invalid arguments")
    sys.exit (1)

in_file = open(argv[1], "r")   # It must be a alignment with equal lengths of sequences
print(argv[1], " is on air")

directory = str(argv[2]) + "_" + str(argv[3])  #Creating a directory
os.mkdir(directory)
os.chdir(directory) 
print(directory, " directory have been created")


#-----------modules-----------------------------------------

round_ = 4  # The round parameter! 
h = float(argv[3]) #The mask parameter!


def slice_module(slice_number, align_slice):   #It takes a slice bar os sequences and write in to the out_file
    align_slice = align_slice.upper()
    A = align_slice.count("A")
    T = align_slice.count("T")
    G = align_slice.count("G")
    C = align_slice.count("C")
    counter_A_pr = round(float( align_slice.count("A") / len(align_slice)), round_)
    counter_T_pr = round(float( align_slice.count("T") / len(align_slice)), round_)
    counter_G_pr = round(float( align_slice.count("G") / len(align_slice)), round_)
    counter_C_pr = round(float( align_slice.count("C") / len(align_slice)), round_)
    counter_N_pr = round(float( align_slice.count("N") / len(align_slice)), round_)
    counter__pr = round(float( align_slice.count("-") / len(align_slice)), round_)
    if ((counter_A_pr or counter_T_pr or counter_G_pr or counter_C_pr or counter_N_pr or counter__pr) == 1) or (counter_N_pr or counter__pr != 0):
        return 0  
    else:
        if A+T == len(align_slice) or A+G == len(align_slice) or A+C == len(align_slice) or T+G == len(align_slice) or T+C == len(align_slice) or G+C == len(align_slice):
            if min(counter_A_pr, counter_T_pr) > h or min(counter_A_pr, counter_G_pr) > h or min(counter_A_pr, counter_C_pr) > h or min(counter_T_pr, counter_G_pr) > h or min(counter_T_pr, counter_C_pr) > h or min(counter_G_pr, counter_C_pr) > h:
                out_slice_file.write(str(slice_number) + "\t" + str(align_slice) + "\n")
        return 1


def nuc_identification_module(line): #looking for nuc and nuc% in line
    array = list(line)
    result = {i: array.count(i) for i in array};
    return result


def duplet_compute_module(duplet_list, A, B):
    A_keys = []
    B_keys = []
    for key in A.keys():
        A_keys = A_keys + [key]
    for key in B.keys():
        B_keys = B_keys + [key]

    p1 = round(A[A_keys[0]]/(A[A_keys[0]] + A[A_keys[1]]), round_)
    p2 = round(A[A_keys[1]]/(A[A_keys[0]] + A[A_keys[1]]), round_)
    q1 = round(B[B_keys[0]]/(B[B_keys[0]] + B[B_keys[1]]), round_)
    q2 = round(B[B_keys[1]]/(B[B_keys[0]] + B[B_keys[1]]), round_)

    p1q1 = str(A_keys[0]) + str(B_keys[0])
    p1q2 = str(A_keys[0]) + str(B_keys[1])
    p2q1 = str(A_keys[1]) + str(B_keys[0])
    p2q2 = str(A_keys[1]) + str(B_keys[1])

    counter_p1q1_pr = round(float( duplet_list.count(p1q1) / len(duplet_list)), round_)
    counter_p1q2_pr = round(float( duplet_list.count(p1q2) / len(duplet_list)), round_)
    counter_p2q1_pr = round(float( duplet_list.count(p2q1) / len(duplet_list)), round_)
    counter_p2q2_pr = round(float( duplet_list.count(p2q2) / len(duplet_list)), round_)
    return counter_p1q1_pr, counter_p1q2_pr, counter_p2q1_pr, counter_p2q2_pr, p1, p2, q1, q2



#--------------import alignment----------------------------------

alignment = AlignIO.read(in_file, "fasta") #Importing all alignments from in file
in_file.close()
out_slice_file = open("out_slice_file_" + directory + ".txt", "w")
out_slice_file.write("№" + "\t" + "Slice" + "\n")

align_slice = ''

print("1. SNPs selection...")
genome_size = len(alignment[0, :])
for slice_number in range(genome_size):  # Taking a nucleotide from each bar/slice in alignment and from each sequence. For nuc in range of length of seq
    for seq in range(len(alignment)):  # for seq from all seqs
        align_slice = align_slice + str(alignment[seq, slice_number]) #Making slices
    a = slice_module(slice_number, align_slice) #Collecting it and putting in the slice_module
    align_slice = ''

out_slice_file.close()    

#----------------------------------------------------------------

slices_file = open("out_slice_file_" + directory + ".txt", "r")
line_id = slices_file.readline() # Take out the heat of the table
slice_id_dic = {}
t_calc = 0
																																																																																																																																																																													
while len(line_id) > 0:
    line_id = slices_file.readline().strip().split()  # Reading a slice in the file
    if len(line_id) != 0:
        slice_id_dic[line_id[0]] = line_id[1]  # Making a dict with seqs and keys
        t_calc += 1

out_slice_file.close()
snps = (len(slice_id_dic))
print("			- Done.")
																																																																											

#---------------------------------------------------------------------------------

print("2. SNPs pairs generation...")

out_2_slices_file = open("out_2_slices_file_" + directory + ".txt", "w")
c = 0
for key in list(slice_id_dic):  # For keys in dict
    one = slice_id_dic[key]  #retrieving a slice 
    for q in list(slice_id_dic):
        if key != q:
            summ = min(abs(int(q) - int(key)), genome_size - abs(int(key) - int(q))) #Calculating the distance
            out_2_slices_file.write(">" + "\t" + str(key) + "\t" + "-" + "\t" + str(q) + "\t" + str(summ) + "\n") #Writing into out_2_slices_file.txt - the temp file - all non snps slices
            out_2_slices_file.write(one + "\n")
            out_2_slices_file.write(slice_id_dic[q] + "\n")
    c += 1 
out_2_slices_file.close() 
print("			- Done.")


#------------------------------------------------------------------------

print("3. LD calculation...")
 
out_2_slices_file = open("out_2_slices_file_" + directory + ".txt", "r") 

LD_out_file = open("LD.txt", "w")
LD_out_file.write("Len" + "\t" + "LD" + "\n")

r2_out_file = open("r2.txt", "w")
r2_out_file.write("Len" + "\t" + "r2" + "\n")


line = out_2_slices_file.readline().strip().split()
calc = 1               


while len(line) != 0:
    if line[0] == ">": 
        l = line[4]
        l = int(l)
        
        line1 = out_2_slices_file.readline().strip()
        line2 = out_2_slices_file.readline().strip()
        A = nuc_identification_module(line1)
        B = nuc_identification_module(line2)

        duplet_list = []
        for i in range(len(line1)):
            duplet_list = duplet_list + [line1[i] + line2[i]]
        AB = duplet_compute_module(duplet_list, A, B) 
 
        LD = AB[0]*AB[3] - AB[1]*AB[2]
        if LD >= 0:
            LD_out = round(LD/(min(AB[4]*AB[7], AB[5]*AB[6])), round_)
        if LD < 0:
            LD_out = round(LD/(max( -AB[4]*AB[6], -AB[5]*AB[7])), round_)   
  
    LD_out_file.write(str(l) + "\t" + str(LD_out) + "\n")

    r = round((LD**2)/(AB[4]*AB[5]*AB[6]*AB[7]), round_) 
    r2_out_file.write(str(l) + "\t" + str(r) + "\n")
    line = out_2_slices_file.readline().strip().split()

out_2_slices_file.close()    
print("			- Done")
toc = time()

print("")
print("")
print("There were ", len(alignment), " sequences in the alignment. It has been found ", snps, "SNPs.")
time_ = round(toc - tic, round_)
if time_ > 60:
    print("Timing:", str(round((toc - tic)/60, round_)), "min.")
else:
    print("Timing:", str(round(toc - tic, round_)), "sec.")
print("")
print("")

