#!/usr/bin/env python

import matplotlib.pyplot as plt
import gzip
import argparse

parser = argparse.ArgumentParser()

# Include input file, sequence length, read number
parser.add_argument('-l', type=int, required=True)
parser.add_argument('-f', type=str, required=True)
parser.add_argument('-r', type=str, required=True)
args = parser.parse_args()

# variables to store command-line inputs
input_file = args.f
sequence_length = args.l
read_number = args.r

def init_list(array, value=0.0):
    '''This function takes an empty list and will populate it with
    the value passed in "value". If no value is passed, initializes list
    with 101 values of 0.0.'''

    # Use a value of 8 for the index files (Read2 and Read3) and 101 for biological reads (Read1 and Read3)
    array = [value] * sequence_length
    return array

mean_scores = []
mean_scores = init_list(mean_scores)

def convert_phred(letter):
    """Converts a single character into a phred score"""
    phredScore = ord(letter) - 33
    return phredScore

# this is the absolute path for the files.
read1 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
read2 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
read3 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
read4 = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
file = input_file

def populate_list(file):
    """Update with your own docstring"""
    
    # create an empy list
    mean_scores = []
    mean_scores = init_list(mean_scores)
    
    # open file and read every 4th line (quality score)
    with gzip.open(file, 'tr') as fh:
        LN = 0
        
        for line in fh:
            LN += 1
            # line = line.decode("utf-8")
            line = line.strip()
            if LN % 4 == 0: # quality score line found
                for counter, char in enumerate(line): # read through each character of the line
                    mean_scores[counter] += convert_phred(char)
    return mean_scores, LN

mean_scores, NR = populate_list(file)

for counter, num in enumerate(mean_scores):
    mean = num / (NR/4) # round this?
    mean_scores[counter] = mean

base_pair_num = []
mean_quality_score = []

for count, score in enumerate(mean_scores):
    base_pair_num.append(count)
    mean_quality_score.append(score)

file_name = "read" + read_number +".png"

plt.rcParams['figure.figsize'] = (11, 5)
plt.bar(base_pair_num, mean_quality_score, color="#444444", alpha = 1)
plt.xlabel("Base Pair Number")
plt.ylabel("Mean Quality Score")
plt.title("Mean Quality Score for Read " + read_number + " Base Pair Positions")
plt.savefig(file_name)