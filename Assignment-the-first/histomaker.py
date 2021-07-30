#!/usr/bin/env python3

import argparse
import numpy 
import matplotlib.pyplot as plt
import gzip
print("pythonworking")


#argparse. 
parser = argparse.ArgumentParser()
parser.add_argument("-seq1", default=1)
parser.add_argument("-seq2", default=1)
parser.add_argument("-index1", default=1)
parser.add_argument("-index2", default=1)
parser.add_argument("-dictionary")
args = parser.parse_args()

seq1 = args.seq1
seq2 = args.seq2
index1 = args.index1
index2 = args.index2
dictfile = args.dictionary

print(seq1)




inputfiles = [seq1, seq2, index1, index2]

def savefig(name, x, y):
    plt.figure(name)
    plt.bar(x, height=y)
    plt.xlabel("Position")
    plt.ylabel("Average Qscore")
    plt.title(name + " qscore dist")
    plt.savefig(name)
    

def convert_phred(string: str):
    #I'm not sure what the difference between convert_phred() and qual_score() is meant to be, but this function either converts a single ascii character or all in a string depending on input
    if len(string) == 1:
        return ord(string) - 33
    else:
        return [ord(x) - 33 for x in string]

index_dictionary = {}
Firstline = True


with open(dictfile, "r") as indexes:
    for line in indexes:
        if Firstline:
            Firstline = False
            continue

        line = line.replace("\n", "").split("\t")
        
        index_dictionary[line[4]] = line[1]

indexseqs = index_dictionary.keys()
indexnames = index_dictionary.values()

""" for indexname in indexnames:
    open(indexname + "_R1.fastq","w")
    open(indexname + "_R2.fastq","w") """
""" open("SWAPS.fastq", "w")
open("BAD_DATA.fastq", "w") """

error_default_names = ["read1", "read2", "index1", "index2"]
errorcounter = 0 

for inputfile in inputfiles:
    print("started one")
    Firstline2 = True
    sums = 1
    records = 0
    linelength = 0
    with gzip.open(inputfile, "rb") as input:
        for i, line in enumerate(input):

            
            
            if i%4 == 3:
                
                line = line.decode()
                line = line.strip()
                records += 1
                if Firstline2:
                    linelength = len(line)
                    sums = numpy.zeros(linelength)
                    Firstline2 = False
                
                sums += numpy.array(convert_phred(line))
    averages = sums/records
    x = [v + 1 for v in range(linelength)]
    inputfile = inputfile.replace("/", "_")
    inputfile = inputfile.replace(".", "")
    try:
        savefig(inputfile, x, averages)
    except:
        savefig(error_default_names[errorcounter], x, averages)
    errorcounter += 1
    

                

    
            
            
        
        


        
        
#To run on test files, use python demultiplex.py -seq1 READ_1_INPUT.fastq -seq2 READ_2_INPUT.fastq -index1 INDEX_1_INPUT.fastq -index2 INDEX_2_INPUT.fastq -dictionary /projects/bgmp/shared/2017_sequencing/indexes.txt
#To run on REAL files, use 

#read 1 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
#read 2 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz
#index 1 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
#index 2 = zcat /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz