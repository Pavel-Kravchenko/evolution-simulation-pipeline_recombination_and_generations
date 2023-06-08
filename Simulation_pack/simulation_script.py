import random
import datetime
from sys import argv
from time import time

tic = time()
# constants
length_of_sequence = 100 #int(input("Type the length of sequence "))
fasta_line_const = 60
round_ = 4  # The round parameter! 

# vars
number_of_species = int(argv[1]) #int(input("Type the number of species in population "))
number_of_cycles = int(argv[2]) #int(input("Type the number of evolution cycles "))
mutations_per_cycle = int(argv[3]) #int(input("Type the number of mutations per cycle "))
is_recombination = "y" #str(input("Do you want to use the heteroplasmy model? [y]/[n] "))
recombination_rank = int(argv[4])

if is_recombination == "y":
    out_file_name = "simulation_heteroplasmy|gen_" + str(number_of_cycles) + "|rec_" + str(recombination_rank) + "|mut_" + str(mutations_per_cycle) + ".fasta"
    recombination_rank = round(recombination_rank*number_of_species//100)
    if (recombination_rank % 2) != 0:
        recombination_rank = recombination_rank -1  

if is_recombination == "n":
    out_file_name = "simulation|gen" + str(number_of_cycles) + "|species_" + str(number_of_species) + "|mut_" + str(mutations_per_cycle) + ".fasta"
print(out_file_name, " have been created")


def make_random_sequence(length, alphabet="ATGC"):
    seq_lst = [] # it's better to use here list comprehensions
    for i in range(length):
        seq_lst.append(random.choice(alphabet))
    seq = "".join(seq_lst) # it is better than appending symbols directly to seq, why? 
    return seq


def sequence_mutation_module(seq, mutations_per_cycle, alphabet="ATGC"): #  for unic mutations
    seq_list = list(seq)
    for i in range(mutations_per_cycle):
        ran_nuc = random.randint(1, len(seq)) - 1
        seq_list[ran_nuc] = random.choice(alphabet)
    seq = "".join(seq_list)
    return seq


def population_dic_maker(seq, number_of_species, mutations_per_cycle, iteration=0, in_log=0):
    calc = 0
    dic = {}
    for i in range(number_of_species):
        name = ">iter_" + str(iteration) + "_in|log_" + str(in_log) + "_sp_" + str(calc)
        calc += 1
        dic[name] = sequence_mutation_module(seq, mutations_per_cycle)
    return dic


def seq_selection(in_dic, number_of_survived_species):
    dic = {}
    random_seq = random.sample(list(in_dic), number_of_survived_species)
    for i in random_seq:
        dic[i] = in_dic[i]
        
    if is_recombination == "y":
        #print("heteroplasmy mode")
        random_seq = random.sample(list(dic), recombination_rank)
        list_seq = []
        for i in range(len(random_seq)):
            if i%2 == 0:
                q = random_seq[i]
                q_1 = random_seq[i + 1]
                seq1 = dic[q]
                seq2 = dic[q_1]
                #print(seq1)
                #print(seq2)
                ran_nuc = random.randint(1, len(seq1)) - 1
                if ran_nuc <= len(seq1)/2:
                    seq1_mut = seq2[:ran_nuc] + seq1[ran_nuc:]
                    seq2_mut = seq1[:ran_nuc] + seq2[ran_nuc:]
                if ran_nuc > len(seq1)/2:
                    seq1_mut = seq1[:ran_nuc] + seq2[ran_nuc:]
                    seq2_mut = seq2[:ran_nuc] + seq1[ran_nuc:]
                #print(seq1_mut)
                #print(seq2_mut)
                list_seq.append(seq1_mut)
                list_seq.append(seq2_mut)
            else:
                pass
        count = 0
        for j in random_seq:
            dic[j] = list_seq[count]
            count += 1
    return dic
    
    
def module_fasta_writer(out_file, line):  # writes in fasta format by 60 nums in each line
    start_point = 0
    end_point = fasta_line_const
    if len(line) > fasta_line_const:
        len_line = len(line)
        while len_line > 0:
            j = line[start_point:end_point]
            len_line = len_line - fasta_line_const
            out_file.write(j + "\n")
            start_point = start_point + fasta_line_const
            end_point = end_point + fasta_line_const
    else:
        out_file.write(line + "\n")
        
        
# -----------------------------------------------------------------------------------

seq = make_random_sequence(length_of_sequence)  # making a random seq

iteration = 1  # creating a population
print("iteration:", iteration)
dic_temp = population_dic_maker(seq, number_of_species, mutations_per_cycle, iteration)
iteration = 2
in_log = 1
for i in range(number_of_cycles):  # making evolution rates
    print("iteration:", iteration)
    species_dic = dic_temp.copy() 
    dic_temp.clear()
    for j in species_dic.keys():
        seq = species_dic[j]
        new_seqs = population_dic_maker(seq, number_of_species, mutations_per_cycle, iteration, in_log)
        dic_temp.update(new_seqs)
        in_log += 1
    dic_temp = seq_selection(dic_temp, number_of_species)
    iteration += 1

with open(out_file_name, 'w') as f:
    count = 0
    for i in dic_temp.keys():
        if count <= 10:
            f.write(str(i) + '\n')
            module_fasta_writer(f, dic_temp[i])
            count += 1 
toc = time()
print("")
print("")
time_ = round(toc - tic, round_)
if time_ > 60:
    print("Timing:", str(round((toc - tic)/60, round_)), "min.")
else:
    print("Timing:", str(round(toc - tic, round_)), "sec.")
print("")
print("")
