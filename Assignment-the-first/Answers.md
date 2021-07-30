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
    2. ```Your answer here```
    3. ```Your answer here```
    
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
   
4. Upload your [4 input FASTQ files](../TEST-input_FASTQ) and your [>=6 expected output FASTQ files](../TEST-output_FASTQ).
5. Pseudocode
   
   Set up env and variables with #!/usr/bin/env python3 and argparse inputs for read1, read2, index1, index2, and index identity file. 
   
   Place input files in a list called inputfiles, arrange in the order seen above.
   
   Declare functions convert_phred and rev_comp.
   
   Create empty index identity dictionary, populate by parsing through index identity file. Sequences will be keys, index buxkets (B1, B9, etc) will be values. 
   
   For each element in the aboce dictionary.values(), create/open 2 files:
        open(key + "_R1.fastq", "w") and open(key + "_R2.fastq", "w")
   Also create and open an additional 4 files (R1 and R2 each) for bad_data and index_swap buckets.  
   
   As we loop through the above logic, set unique counters for each bucket to 0.
   
   
   Create a dictionary called stored_lines.
   keys will be 0:3, values will all be set to None
   
   Open all inputs (with gz.open file as seq1, file2 as index1, etc.)
   Begin iterating throug each line of zip(read1, read2, index1, index2):
   
   
   For each i, line in enumerate(zipped files):
        stored_lines[i%4] = tuple(line)
        if i%4 == 0:
            indices = stored_lines[1][2:3]
            qscores = 
   
   
   
   
7. High level functions. For each function, be sure to include:
    1. Description/doc string
    2. Function headers (name and parameters)
    3. Test examples for individual functions
    4. Return statement
