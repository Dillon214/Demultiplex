#!/usr/bin/env python3

import argparse
from os import sep, write
import gzip
import numpy

print("hello")

 
parser = argparse.ArgumentParser()
parser.add_argument("-seq1", default=1)
parser.add_argument("-seq2", default=1)
parser.add_argument("-index1", default=1)
parser.add_argument("-index2", default=1)
parser.add_argument("-dictionary")
parser.add_argument("-indexqscore", default=0)
parser.add_argument("-readqscore", default=0)
parser.add_argument("-rev_comp_indexes_in_output", default = False)

args = parser.parse_args()

seq1 = args.seq1
seq2 = args.seq2
index1 = args.index1
index2 = args.index2
dictfile = args.dictionary
readqthreshold = args.readqscore
indexqthreshold = args.indexqscore
rev_comp_indexes_in_output_files = args.rev_comp_indexes_in_output


inputfiles = [seq1, seq2, index1, index2]



def convert_phred(string: str):
    #I'm not sure what the difference between convert_phred() and qual_score() is meant to be, but this function either converts a single ascii character or all in a string depending on input
    if len(string) == 1:
        return ord(string) - 33
    else:
        return tuple(ord(x) - 33 for x in string)


def reverse_comp(inputseq:str):
    transformer = {"A":"T", "T":"A", "G":"C", "C":"G"}
    return ("".join((transformer[x] for x in inputseq)))[::-1]

def convert_phred2(string: str):
    
    
    return numpy.array(tuple(ord(x) - 33 for x in string))




index1_identity_dict = {}
index2_identity_dict = {}
index_MATCHES = {}

with open(dictfile, "r") as indexes:
    for i, line in enumerate(indexes):
        if i == 0:
            continue
        line = line.strip()
        lineelements = line.split(sep = "\t")
        index1_identity_dict[lineelements[-1]] = lineelements[1]
        index2_identity_dict[reverse_comp(lineelements[-1])] = lineelements[1]
        index_MATCHES[(lineelements[-1], reverse_comp(lineelements[-1]))] = lineelements[1]
#print(index1_identity_dict)
#print(index2_identity_dict)

        


indexnames = index1_identity_dict.values()
index1_seqs = index1_identity_dict.keys()
index2_seqs = index2_identity_dict.keys()
COUNTERS = {}
OUTPUT_DICT = {}
DESTINY_COUNTER = {}
for i, indexname in enumerate(indexnames):
    OUTPUT_DICT[indexname + "R1"] = open(indexname + "_R1.fastq","w") 
    OUTPUT_DICT[indexname + "R2"] = open(indexname + "_R2.fastq","w")
    DESTINY_COUNTER[indexname] = 0
OUTPUT_DICT["SWAPS_R1"] = open("SWAPS_R1.fastq", "w")
OUTPUT_DICT["SWAPS_R2"] = open("SWAPS_R2.fastq", "w")
DESTINY_COUNTER["SWAPS_"] = 0
OUTPUT_DICT["BAD_DATA_R1"] = open("BAD_DATA_R1.fastq", "w")
OUTPUT_DICT["BAD_DATA_R2"] = open("BAD_DATA_R2.fastq", "w")
DESTINY_COUNTER["BAD_DATA_"] = 0




#INDEX_PAIR_COUNTER_ARRAY = numpy.










inputfiles_unzipped = []
for inputfile in inputfiles:
    if inputfile[-1] == "z":
        inputfiles_unzipped.append(gzip.open(inputfile, "rt"))
    else:
        inputfiles_unzipped.append(open(inputfile, "r")) 


ALL_RECORDS_COUNTER = 0
ALL_UNIQUE_INDEXES = {}



print(index_MATCHES)

#stored_lines = {0:None, 1:None, 2:None, 3:None}
stored_lines = {}
for i,line in enumerate(zip(*inputfiles_unzipped)):
    line = tuple(i[:-1] for i in line)
    
    stored_lines[i%4] = line
    
    if i%4 == 3:
        ALL_RECORDS_COUNTER += 1

        indices = (stored_lines[1][2:])
        
        readqscores = convert_phred("".join(stored_lines[3][:2]))
        indexqscores = convert_phred("".join(stored_lines[3][2:]))
        read1 = "\n".join(((stored_lines[0][0] + "-".join(indices)),stored_lines[1][0],"+",stored_lines[3][0], "\n"))
        read2 = "\n".join(((stored_lines[0][1] + "-".join(indices)),stored_lines[1][1],"+",stored_lines[3][1], "\n"))
        
        

        if min(readqscores) > readqthreshold and min(indexqscores) > indexqthreshold and stored_lines[1][2] in index1_seqs and stored_lines[1][3] in index2_seqs:
            destiny = index_MATCHES.get(indices, "SWAPS_")
            
        else:
            destiny = "BAD_DATA_"

        
        
        DESTINY_COUNTER[destiny] += 1
        OUTPUT_DICT[destiny + "R1"].writelines(read1)
        OUTPUT_DICT[destiny + "R2"].writelines(read2)
        

        if destiny == "SWAPS_":
            try:
                ALL_UNIQUE_INDEXES[indices] += 1
            except:
                ALL_UNIQUE_INDEXES[indices] = 0

print(ALL_UNIQUE_INDEXES)
print(DESTINY_COUNTER)





            
        
 
        

        







        

            


            


#remaining possibilities
#index swapping
#not in indexes
#good

        #stored_lines.fromkeys(stored_lines, 0)





#TACGCTAC
#GTAGCGTA



        
        
#To run on test files, use python demultiplex.py -seq1 READ_1_INPUT.fastq -seq2 READ_2_INPUT.fastq -index1 INDEX_1_INPUT.fastq -index2 INDEX_2_INPUT.fastq -dictionary /projects/bgmp/shared/2017_sequencing/indexes.txt
#To run on REAL files, use 

#read 1 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#read 2 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
#index 1 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#index 2 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
