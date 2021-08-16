#!/usr/bin/env python3

import argparse
from os import sep, write
import gzip
from datetime import datetime
import matplotlib.pyplot as plt




#Runs in 1 hour 40 minutes on 1 cpu under default setting


#Functions and assertions below. Program will not run unless asserions are passed. 

#Convert phred is called in the main loop. Takes in an iterable, returns the ord() - 33 or each value in said iterable in a tuple. 
def convert_phred(string):
    return tuple(ord(x) - 33 for x in string)
assert(convert_phred("JJJ") == (41,41,41))

#Reverse comp is NOT called in main loop, so inefficiency is more tolerated here by me. Simply returns the reverse compliment of a string of nucleotides. 
transformer = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}
def reverse_comp(inputseq:str):
    
    return ("".join((transformer[x] for x in inputseq)))[::-1]
assert(reverse_comp("AAAT") == "ATTT")

print("assertions passed")







#Argparse. I chose to keep the read and index inputs seperate (as opposed to all in one) so that the user can easily see ehat inputs are what. 
#Options for qscore cutoffs and one for preference of viewing index 2 (either same as index or revcomp). Electing for this option incurs a timecost. 
parser = argparse.ArgumentParser()
parser.add_argument("-read1", default=1, help="path to read1 file")
parser.add_argument("-read2", default=1, help="path to read2 file")
parser.add_argument("-index1", default=1, help="path to index1 file")
parser.add_argument("-index2", default=1, help="path to index1 file")
parser.add_argument("-dictionary", help="path to index identity file")
parser.add_argument("-indexqscore", default=0, help="optional index qscore cutoff. default is 0")
parser.add_argument("-readqscore", default=0, help="optional read qscore cutoff. default is 0")
parser.add_argument("-indexmatch", default = False, help="input anything to ")


args = parser.parse_args()

seq1 = args.read1
seq2 = args.read2
index1 = args.index1
index2 = args.index2
dictfile = args.dictionary
readqthreshold = args.readqscore
indexqthreshold = args.indexqscore
rev_comp_indexes_in_output_files = args.indexmatch


inputfiles = [seq1, seq2, index1, index2]



#Below 4 lines open a file for stats on data/program, adds it to a list of open files, and writes the time to that file.
ALL_OPEN_FILES = []
statoutput = open("OUTPUT_INFORMATION.md", "w")
ALL_OPEN_FILES.append(statoutput)
statoutput.writelines("Started at " + str(datetime.now()) + "\n\n")







#ATTENTION 
#I elected to generally focus on the main loop's optimization over anything outside it, so be warned.




#The mess of code below is not worth fully explaining detail, but suffice it to say that it parses through our index identity file and makes our program faster.
#Instead of checking whether index 2 is the reverse complement of index 1 in our loop, we instead generate a combinatorial dictionary.
#This dictionary is structured like so ... key: (validindex1, reverse compliment of valid index 1) value: index IDENTITY (2A, 3H, etc)
#If (actual index1, actual index 2) is in this dictionary, we know a couple of things...
#We know that the reads associated with this index-pair fall into the bucket defined by the value
#We know that index swapping did NOT occur 

#So, if we can prove beforehand that the indexes we're looking at aren't BAD DATA, we can prove a case of index swapping.
#We do this by simply checking whether index 1 is in index1_identity_dict keys and index 2 is in index2_identity_dict keys
#Long story short, this is the code we use to offload some stress from our actual main loop. We use information we do know to infer what we can, and find what we can't know here in the loop. 

index1_identity_dict = {}
index2_identity_dict = {}
index_MATCHES = {}
code_indexes = {}
with open(dictfile, "r") as indexes:
    ALL_OPEN_FILES.append(indexes)
    for i, line in enumerate(indexes):
        if i == 0:
            continue
        line = line.strip()
        lineelements = line.split(sep = "\t")
        index1_identity_dict[lineelements[-1]] = lineelements[1]
        index2_identity_dict[reverse_comp(lineelements[-1])] = lineelements[1]

        if rev_comp_indexes_in_output_files:
            index_MATCHES[(lineelements[-1], lineelements[-1])] = lineelements[1]
            code_indexes[lineelements[1]] = (lineelements[-1], lineelements[-1])
        else:
            index_MATCHES[(lineelements[-1], reverse_comp(lineelements[-1]))] = lineelements[1]
            code_indexes[lineelements[1]] = (lineelements[-1], reverse_comp(lineelements[-1]))


        

#This code below automatically opens up R1 and R2 files for our buckets and stores them in a ditionary.
#We alse set up "destiny counter", which simply records how records are sorted into which buckets.
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
    ALL_OPEN_FILES += (OUTPUT_DICT[indexname + "R1"], OUTPUT_DICT[indexname + "R2"])
OUTPUT_DICT["SWAPS_R1"] = open("SWAPS_R1.fastq", "w")
OUTPUT_DICT["SWAPS_R2"] = open("SWAPS_R2.fastq", "w")
DESTINY_COUNTER["SWAPS_"] = 0
OUTPUT_DICT["BAD_DATA_R1"] = open("BAD_DATA_R1.fastq", "w")
OUTPUT_DICT["BAD_DATA_R2"] = open("BAD_DATA_R2.fastq", "w")
DESTINY_COUNTER["BAD_DATA_"] = 0

ALL_OPEN_FILES += (OUTPUT_DICT["SWAPS_R1"], OUTPUT_DICT["SWAPS_R2"], OUTPUT_DICT["BAD_DATA_R1"], OUTPUT_DICT["BAD_DATA_R2"])












#This code below makes it so unzipped diles get left alone, and zipped files get unzipped. 
inputfiles_unzipped = []
for inputfile in inputfiles:
    if inputfile[-1] == "z":
        inputfiles_unzipped.append(gzip.open(inputfile, "rt"))
        ALL_OPEN_FILES.append(gzip.open(inputfile, "rt"))
    else:
        inputfiles_unzipped.append(open(inputfile, "r")) 
        ALL_OPEN_FILES.append(open(inputfile, "r"))



#All records counter keeps track of total records
#All swap indexes keeps track of each unique index swapping event and their occurances. We track this in-loop. 
ALL_RECORDS_COUNTER = 0
ALL_SWAP_INDEXES = {}






#This next part is where we begin our loop. 
#It took some trying/experiementation to determine a good approach for temporarily storing data as we iterate (numpy & dictionary were both busts) but I settled on simple variables which concatenate once per line.
#These variables all change once every 4 lines/record. 
#Destiny is the same as bucket


index1 = ""
index2 = ""
read1 = ""
read2 = ""
destiny = ""

#zip is a really neat little way to iterate which I just discovered during this project. It allows multiple iterators, such as our 4 files, to be iterated through simultaneously. Each iteration produces a tuple of lines.

#The below if statement is for user output customization. While difficult to read, calling this one if outside the loop prevents us from calling it inside, and that's preferable at over a billion repeats. 
if not rev_comp_indexes_in_output_files:
    for i,line in enumerate(zip(*inputfiles_unzipped)):
    
        #This is the part of our script visited most often, so it must be brief. This one grows 4 fastq records from our 4 files one line at a time, then operates on them once complete. Error underlines are liars.
        read1 += line[0]
        read2 += line[1]
        index1 += line[2]
        index2 += line[3]
        
        
        if i%4 == 3:
            #With 4 fastq records, one from each file, prepared, now we can act on them. 

            ALL_RECORDS_COUNTER += 1

            #We split them up for easily accessing qscores, reads, whatever we need. 
            read1split = read1.split("\n")
            read2split = read2.split("\n")
            index1split = index1.split("\n")
            index2split = index2.split("\n")

            #Here we have a tuple of indexes which we use to access values in index_MATCHES, which is how we find out our bucket. NO REVERSE COMPLIMENT NECESSARY (unless the user asks for it)
            indexes = (index1split[1], index2split[1])

            #Here we create our index pair for adding to header. 
            index_for_reads = "".join((" ", indexes[0], "-", indexes[1], "\n"))
            
            #Big check for BAD DATA time. Here we check that index quality scores exceed thresholds and that both indexes are valid (though not necessarily a valid PAIR, keep that in mind)
            #Failing a condition here results in destiny/bucket becoming BAD DATA
            #Becuase python makes each check individually, we don't generate the lowest qscores before this step, as we can avoid generating the second in the case that the first is failed. 
            if min(convert_phred(set("".join((read1split[3], read2split[3]))))) > readqthreshold and min(convert_phred(set("".join((index1split[3], index2split[3]))))) > indexqthreshold and indexes[0] in index1_seqs and indexes[1] in index2_seqs:
                #This distinguishes between remaining 2 buckets: either valid pair (resulting in assigning that pair's correcponding index identity bucket) or index swapping. 
                destiny = index_MATCHES.get(indexes, "SWAPS_")
            else:
                destiny = "BAD_DATA_"
            
            #Slotting our index pair in. 
            read1 = read1.replace("\n", index_for_reads, 1)
            read2 = read2.replace("\n", index_for_reads, 1)

            #Tally destiny counter. 
            DESTINY_COUNTER[destiny] += 1
            
            #Writing tp files. I learned that writelines and write have a huge difference in testing this section. 
            OUTPUT_DICT[destiny + "R1"].write(read1)
            OUTPUT_DICT[destiny + "R2"].write(read2)

            #We only catalog occurances of unique pairs for swaps here. We already know occurances of vali indexes from destiny counter, so we only record what is unknown. 
            #Note, we could even skip this "except" by doing some additional pregeneration, but the cost it incurs naturally decreases as the program works anyways, so I won't bother.  
            #We use try/except to avoid calling if/else. Excepts naturally decrease over time, so time cost of it is marginal. 
            if destiny == "SWAPS_":
                try:
                    ALL_SWAP_INDEXES[indexes] += 1
                except:
                    ALL_SWAP_INDEXES[indexes] = 1

            #Reset temporary storage
            index1, index2, read2, read1 = "", "", "", ""
else:
    for i,line in enumerate(zip(*inputfiles_unzipped)):
        read1 += line[0]
        read2 += line[1]
        index1 += line[2]
        index2 += line[3]
        if i%4 == 3:
            ALL_RECORDS_COUNTER += 1
            read1split = read1.split("\n")
            read2split = read2.split("\n")
            index1split = index1.split("\n")
            index2split = index2.split("\n")
            indexes = (index1split[1], reverse_comp(index2split[1])) 
            
            index_for_reads = "".join((" ", indexes[0], "-", indexes[1], "\n"))
            if min(convert_phred(set("".join((read1split[3], read2split[3]))))) > readqthreshold and min(convert_phred(set("".join((index1split[3], index2split[3]))))) > indexqthreshold and indexes[0] in index1_seqs and indexes[1] in index1_seqs:  
                destiny = index_MATCHES.get(indexes, "SWAPS_")
            else:
                destiny = "BAD_DATA_"
            read1 = read1.replace("\n", index_for_reads, 1)
            read2 = read2.replace("\n", index_for_reads, 1)
            DESTINY_COUNTER[destiny] += 1
            OUTPUT_DICT[destiny + "R1"].write(read1)
            OUTPUT_DICT[destiny + "R2"].write(read2)
            if destiny == "SWAPS_":
                try:
                    ALL_SWAP_INDEXES[indexes] += 1
                except:
                    ALL_SWAP_INDEXES[indexes] = 1
            index1, index2, read2, read1 = "", "", "", ""



        

#Below we report stats. I chose to report percentages of records goinf into each bucket, percent of total valid records, occurances of valid index pairs, and occurances of index swapping pairs.
#All stats both printed and recorded to a file. 

#ATTENTION: I messed up a '/n' on one of these, leading to a section of awkward markdown. Fizxed here, but I didn't have the hour +40 minutes to both fix and rerun. SORRY!

ALL_VALID_INDEXES = {}
Totalvalid = 0

print("---RESULTS---\n\n")
statoutput.writelines("---RESULTS---\n\n")
print(DESTINY_COUNTER)


for index in DESTINY_COUNTER.keys():
    if index in code_indexes.keys():
        ALL_VALID_INDEXES[code_indexes[index]] = DESTINY_COUNTER[index]
        Totalvalid += DESTINY_COUNTER[index]


print("\nTOTAL READS SORTED = " + str(ALL_RECORDS_COUNTER*2) + "\n")
statoutput.writelines("\nTOTAL READS SORTED = " + str(ALL_RECORDS_COUNTER*2) + "\n")
for destiny in DESTINY_COUNTER.keys():
    print(destiny + " bucket received " + str((DESTINY_COUNTER[destiny]/ALL_RECORDS_COUNTER)*100) + " percent of all reads\n")
    statoutput.writelines(destiny + " bucket received " + str((DESTINY_COUNTER[destiny]/ALL_RECORDS_COUNTER)*100) + " percent of all reads\n\n")

print("\n" + str((Totalvalid/ALL_RECORDS_COUNTER)*100) + " percent of all entries had a valid index pair\n")
statoutput.writelines("\n" + str((Totalvalid/ALL_RECORDS_COUNTER)*100) + " percent of all entries had a valid index pair\n")



print("\n\nTHE OCCURANCES OF RECOGNIZED INDEX COMBINATIONS THAT WERE NOT DEEMED 'BAD DATA' ARE...\n")
statoutput.writelines("\n\nTHE OCCURANCES OF RECOGNIZED INDEX COMBINATIONS THAT WERE NOT DEEMED 'BAD DATA' ARE...\n")
for valid_indexpair in ALL_VALID_INDEXES.keys():
    print(str(valid_indexpair) + " (" + index_MATCHES[valid_indexpair] + ") occured " + str(ALL_VALID_INDEXES[valid_indexpair]) + " times\n")
    statoutput.writelines(str(valid_indexpair) + " (" + index_MATCHES[valid_indexpair] + ") occured " + str(ALL_VALID_INDEXES[valid_indexpair]) + " times\n\n")

print("\n\nTHE OCCURANCES OF INDEX COMBINATIONS RESULTING FROM INDEX SWAPPING AND NOT DEEMED 'BAD DATA' ARE...\n")
statoutput.writelines("\n\nTHE OCCURANCES OF INDEX COMBINATIONS RESULTING FROM INDEX SWAPPING AND NOT DEEMED 'BAD DATA' ARE...\n\n")
for swap_indexpair in ALL_SWAP_INDEXES.keys():
    print(str(swap_indexpair) +   "occured " + str(ALL_SWAP_INDEXES[swap_indexpair]) + " times\n")
    statoutput.writelines(str(swap_indexpair) +   "occured " + str(ALL_SWAP_INDEXES[swap_indexpair]) + " times\n\n")




#Record ending time.
statoutput.writelines("\n\nEnded at " + str(datetime.now()))

#Close files.    
for file in ALL_OPEN_FILES:
    file.close()

#Making a graph for viewing bucket data
plt.bar(x = DESTINY_COUNTER.keys(), height= DESTINY_COUNTER.values())
plt.xlabel("Bucket")
plt.ylabel("Records")
plt.xticks(rotation = 90, size = 10)
plt.tight_layout()
plt.savefig("BUCKET_GRAPH")

#As an aside, I'm very glad we were given bonus work time. I struggled for a while with a bad approach because I didn't feel I had the time to change directions and experiment, but I got that chance and learned alot from it. 

