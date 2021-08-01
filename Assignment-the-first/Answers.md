# Assignment the First

## Part 1
1. Be sure to upload your Python script.

| File name | label |
|---|---|
| 1294_S1_L008_R1_001.fastq.gz | Read 1 |
| 1294_S1_L008_R2_001.fastq.gz | Index 1 |
| 1294_S1_L008_R3_001.fastq.gz | Index 2 |
| 1294_S1_L008_R4_001.fastq.gz | Read 2 |

2. Per-base NT distribution
    1. Use markdown to insert your 4 histograms here.
    ![_projects_bgmp_shared_2017_sequencing_1294_S1_L008_R1_001fastqgz](https://user-images.githubusercontent.com/73907611/127631626-0cf892f9-ac1e-4561-9bb0-5e3c4d2eedd8.png)
    ![_projects_bgmp_shared_2017_sequencing_1294_S1_L008_R2_001fastqgz](https://user-images.githubusercontent.com/73907611/127631648-516f1884-5f87-4def-91a5-7a768c7dfa30.png)
    ![_projects_bgmp_shared_2017_sequencing_1294_S1_L008_R3_001fastqgz](https://user-images.githubusercontent.com/73907611/127631662-45919bd5-78ad-4524-bd7f-a94e17f39291.png)
    ![_projects_bgmp_shared_2017_sequencing_1294_S1_L008_R4_001fastqgz](https://user-images.githubusercontent.com/73907611/127631677-b639eec6-421e-4413-8111-4dfd5333d75b.png)

    3. I believe that 30 is a good threshold for read quality scores. The read sequences by far the greatest unknown factor in any sequencing exercise, so it
       is imperitive that the information we gather about them is accurate and not misrepresentative.  
       On the pther hand, Index sequences are a little less mysterious. We already hold them to a strict standard in having to match our index identity matrix,
       so the only danger coming from low quality score is the off chance that a matching and seemingly normal index is, in reality, not. For indexes, I think
       a more relaxed threshold of 20 is decent. 
    4. I used "for FILE in /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz; do echo $(zcat $FILE | sed -n "2~4p" | grep -c "N"); done"
       All up, there were indexes 3976613 containing an N in each file.
    
## Part 2
1. Define the problem
    Illumina has produced our sequencing data, but it's all unsorted! To help us sort our reads into the proper index buckets, we have 4 files to work with. 2 contain indexes 
    1 and 2 for each sequence, and another 2 contain reads 1 and 2 for each sequence. We have to use the information found in the index files, and our index identity library, in
    order to sort our data.
    
    There are some additional problems to consider. 
    Firstly, index swapping is a common error in Illumina library prep, and has certinly affected our data. We need to bucket reads with swapped indexes separately.
    Secodly, there are almost certainly problems with our data. Low quality scores, occurances of unknown bases (N), and indexes that don't match any of our index identities.
    These ones must all be sorted into a low-quality bucket. 
    
2. Describe output
   The output of this program will include the number of occurances for each pair of matching indexes that we observe in our data, including index-swapped ones, the number of
   records bearning an unknown index, the number of records which include an 'N', and the number of records which do not meet a quality threshold. 
   
   Additionally, the program will produce 52 files. 2 for read 1 and read 2 of each index bucket (2 * 24), 2 for read 1 and read 2 of our index swapping bucket, and 2 for read
   1 and read 2 of our low-quality data bucket.
   
3. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
4. Pseudocode
   
   Set up env and variables with #!/usr/bin/env python3 and argparse inputs for read1, read2, index1, index2, and index identity file. 
   Also accept optional values for read_threshold and index_threshold
   
   Place input files in a list called inputfiles, arrange in the order seen above.
   
   Declare functions convert_phred, rev_comp
   
   
   convert_phred(inputPhred: string) -> Tuple:
        Takes a Phred score string, returns a Tuple of numerical, converted scores
   
   Convert phred("III") -> (40,40,40)
   
   rev_comp(inputSEQ: string) -> string:
        Takes a DNA nucleotide string, returns its reverse compliment 
   
   rev_comp(AAT) -> ATT
   
   
   Create empty index identity dictionary, populate by parsing through index identity file. Sequences will be keys, index buckets (B1, B9, etc) will be values. 
   Create a second dictionary from the above one by using rev_conp() on the keys. Alternitively, create this one in parallel to the first while parsing.
   
   These dictionaries are called index_1_identities and index_2_identities
   
   For each element in the aboce dictionary.values(), create/open 2 files:
        open(key + "_R1.fastq", "w") and open(key + "_R2.fastq", "w")
   Also create and open an additional 4 files (R1 and R2 each) for bad_data and index_swap buckets.  
   
   As we loop through the above logic, set unique counters for each bucket to 0 in a dictionary. Add new entries for index hopping index pairs. 
   
   
   Create a dictionary called stored_lines.
   keys will be 0:3, values will all be set to None
   
   
   
   Open all inputs (with gz.open file as seq1, file2 as index1, etc.)
   Begin iterating throug each line of zip(read1, read2, index1, index2):
   
   
   For each i, line in enumerate(zipped files):
        line = (i.strip() for i in line)
        stored_lines[i%4] = line
        if i%4 == 0:
        
           First we need to define some variables we know we'll need, so we use the dictionary we just made to piece them together.
            
            indices = stored_lines[1][2:]
            readqscores = "".join(stored_lines[3][0], stored_lines[3][1])
            indexqscores = "".join(stored_lines[3][2], stored_lines[3][3])
  
            read1 = "\n".join(((stored_lines[0][0] + "-".join(indices)),stored_lines[1][0],"+",stored_lines[3][0]))
            read2 = "\n".join(((stored_lines[0][1] + "-".join(indices)),stored_lines[1][1],"+",stored_lines[3][1]))
            
            check first if the qscores (reads and indexes) contain anything that falls below our thresholds defined in our inputs by calling phred_score()
            or if there are any Ns in the sequences:
                   if so, write read 1 and read 2 we generated above into our Bad data bucket files and incremennt that bucket's counter by 1
                   
            Next check if index 1 is in index_1_identities.keys() and index 2 is in index 2 identities.keys():
                if so, check whether index_1_identities[index1] == index_2_identities[index2]
                    if so, write read1 and read2 into the bucket files under index_1_identities[index1], and incrememnt that bucket's counter
                if not, write read1 and read2 into the swapping bucket files, and increment the index pair's unique counter by 1
            if not, write read 1 and read 2 we generated above into our Bad data bucket files and incremennt that bucket's counter by 1
            
            restore all values in our stored_lines dictionary to None using dict.fromkeys
        
                   
                    
    Once we reach the end of our massive loop, we can return the tallies we counted for each index pair/bucket
    
    My thoughts: I will avoid calling convert phred in loop by instead converting quality thresholds to ascii.
    
    Accepted edits: Beagan suggests using sets for quickening logic. I think this is a good change, and will be implementing.
   
   
5. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement



   convert_phred(inputPhred: string) -> Tuple:
        Takes a Phred score string, returns a Tuple of numerical, converted scores
        return(phred_score)
   Convert phred("III") -> (40,40,40)
   
   rev_comp(inputSEQ: string) -> string:
        Takes a DNA nucleotide string, returns its reverse compliment 
        return(reverse compliment)
   rev_comp(AAT) -> ATT
   
   
   
   
