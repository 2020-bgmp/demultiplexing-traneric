#!/usr/bin/env python

from itertools import islice
from itertools import permutations
import gzip
import argparse

# Create parser object.
parser = argparse.ArgumentParser()

# Include input file, sequence length, read number
parser.add_argument('-c', type=int, required=True)
parser.add_argument('-f', type=str, required=True)
parser.add_argument('-r1', type=str, required=True)
parser.add_argument('-r2', type=str, required=True)
parser.add_argument('-i1', type=str, required=True)
parser.add_argument('-i2', type=str, required=True)
args = parser.parse_args()

# variables to store command-line inputs
cutoff_score = args.c
indexes_path = args.f
read1_path = args.r1
read2_path = args.r2
index1_path = args.i1
index2_path = args.i2

def create_file_dictionaries(indexes_path):
    '''
    This function returns a list with the names of the indexes, a dictionary with index names as keys and read1 file names as values,
    and a dictionary with index names as keys and read2 file names as values.
    '''
    with open(indexes_path, 'r') as f:
        index_list = []
        final_index_list = []
        forward_files = {}
        reverse_files = {}

        for index, line in enumerate(f):
            if index == 0:
                continue
            line_list = line.split()
            line_list = line_list[3:5]
            index_list.append(line_list)

        for index in index_list:
            identifier = index[0]
            sequence = index[1]

            forward_files[sequence] = identifier + "_read1.fq"
            reverse_files[sequence] = identifier + "_read2.fq"
            final_index_list.append(index[1])
    return final_index_list, forward_files, reverse_files

def reverse_complement(dna_string):
    '''Takes a string of DNA characters and determine the reverse complement.'''

    reverse_complement = ""

    for base in dna_string:

        if (base == 'A'):
            reverse_complement += 'T'
        elif (base == 'T'):
            reverse_complement += 'A'
        elif (base == 'G'):
            reverse_complement += 'C'
        else:
            reverse_complement += 'G'
    
    reverse_complement = reverse_complement[::-1]
            
    return reverse_complement

def is_known_index(index1):
    '''Takes the index from read2.txt and determines if they exist within the 24 libraries in the 24_index_list'''

    if 'N' in index1:
        return False

    if index1 in index_list:
        return True
    else:
        return False

def is_index_pair(index1, index2):
    '''Takes two indexes and determines if the reverse complement of one index matches the other index.'''

    index2_rev_comp = reverse_complement(index2)

    if index2_rev_comp == index1:
        return True
    else:
        return False

def convert_phred(letter):
    """Converts a single character into a phred score"""
    phredScore = ord(letter) - 33
    return phredScore

def mean_quality_score(quality_line):
    '''Calculates the mean quality score of a given quality sequence'''
    char_count = 0
    quality_score = 0

    for char in quality_line:
        quality_score += convert_phred(char)
        char_count += 1
    quality_score = quality_score / char_count
    return quality_score

def is_low_quality(quality_line, cutoff_score):
    '''Checks to see if the read is above (good quality) or below (low quality) the specified threshhold'''
    quality_cutoff = cutoff_score

    if mean_quality_score(quality_line) < 20:
        return True
    return False

def permutations_dictionary(index_list):
    '''Creates a permutation dictionary from the indexes.'''
    permutations_list = list(permutations(index_list, 2))
    permutations_dict = {}

    for tup in permutations_list:
        permutations_dict[tup] = 0
    return permutations_dict

def create_index_count_dict(index_list):
    '''Creates a dictionary to store the index counts.'''
    index_count_dictionary = {}

    for index in index_list:
        index_count_dictionary[index] = 0
    return index_count_dictionary

def create_output_files(index_count_dictionary, num_matched_indexes, num_unknown_index, num_index_hopping, permutations_dict):
    '''
    This function creates a statistics report and a report for index hop counts.
    '''

    # create demultiplex report
    with open("report.txt", "a") as f:

        total_indexes = num_matched_indexes + num_unknown_index + num_index_hopping
        match_percentage = str(round(((num_matched_indexes / total_indexes) * 100), 2))
        unknown_percentage = str(round(((num_unknown_index / total_indexes) * 100), 2))
        ihop_percentage = str(round(((num_index_hopping / total_indexes) * 100), 2))

        f.write("Matched Indexes: " + str(num_matched_indexes) + " (" + match_percentage + "%)" + "\n")
        f.write("Unknown or Low Quality Indexes: " + str(num_unknown_index) + " (" + unknown_percentage + "%)" + "\n")
        f.write("Index Hopping: " + str(num_index_hopping) + " (" + ihop_percentage + "%)" + "\n\n")

        # write out occurences of each index.
        f.write("Index Percentage\n")
        for key in index_count_dictionary:
            percentage = str(round(((index_count_dictionary[key] / num_matched_indexes) * 100), 2))
            line = key + " " + percentage + "\n"
            f.write(line)
    
    # create a report for index hops
    with open("index_hop_counts.txt", "a") as f:

        for key in permutations_dict:
            value = str(permutations_dict[key])
            key = str(key)
            f.write(key + " " + value + "\n")

# Create file name dictionaries and list
index_list, forward_files, reverse_files = create_file_dictionaries(indexes_path)

# Create a dictionary to track the occurence of each matching index
index_count_dictionary = create_index_count_dict(index_list)

# create dictionaries to hold index-hop count and index-pair occurences
permutations_dict = permutations_dictionary(index_list)

# Create variables to hold the total number of matched indexes, index hopping, and unknown or low quality indexes.
num_matched_indexes = 0
num_index_hopping = 0
num_unknown_index = 0

# open files
read1_file = gzip.open(read1_path, 'rt')
read2_file = gzip.open(read2_path, 'rt')
index1_file = gzip.open(index1_path, 'rt')
index2_file = gzip.open(index2_path, 'rt')

RECORD_LENGTH = 4

while True:

    # slice out the first four lines (record) of each file
    num_lines = 4

    read1_record = list(islice(read1_file, num_lines))
    read2_record = list(islice(read2_file, num_lines))
    index1_record = list(islice(index1_file, num_lines))
    index2_record = list(islice(index2_file, num_lines))

    # break if the list is empty. we'll use index1_record to check since all 4 files have the same line count.
    if not index1_record:
        break

    # Extract the sequence lines from each record.
    read1_sequence = read1_record[1].strip()
    read2_sequence = read2_record[1].strip()
    index1_sequence = index1_record[1].strip()
    index2_sequence = index2_record[1].strip()

    # Extract the quality lines from each index.
    index1_qread = index1_record[3].strip()
    index2_qread = index2_record[3].strip()

    # this variable holds the current index-pair header
    index_pair = index1_sequence + '-' + index2_sequence

    # index2_sequence rev comp.
    index2_rev = reverse_complement(index2_sequence)

    # check for an unknown index read (does not match the 24) or below the quality score cutoff.
    if ( (is_known_index(index1_sequence) == False) or (is_known_index(index2_rev) == False) or (is_low_quality(index1_qread, cutoff_score) == True) or (is_low_quality(index2_qread, cutoff_score) == True) ):

        num_unknown_index += 1

        with open("unknown_read1.fq", "a") as f1, open("unknown_read2.fq", "a") as f2:

            for i in range(RECORD_LENGTH):

                # append the index-pair to the header of BOTH reads.
                if (i == 0):
                    r1_header = read1_record[0].strip() + " " + index_pair + "\n"
                    r2_header = read2_record[0].strip() + " " + index_pair + "\n"
        
                    f1.write(r1_header)
                    f2.write(r2_header)
                else:
                    f1.write(read1_record[i])
                    f2.write(read2_record[i])
        continue

    # check for a matched index-pair
    if (is_index_pair(index1_sequence, index2_sequence)):

        num_matched_indexes += 1

        # track the occurrences of each index.
        index_count_dictionary[index1_sequence] += 1

        forward_file_name = forward_files[index1_sequence]
        reverse_file_name = reverse_files[index2_rev]

        with open(forward_file_name, "a") as f1, open(reverse_file_name, "a") as f2:

            for i in range(RECORD_LENGTH):

                # append the index-pair to the header of BOTH reads.
                if (i == 0):
                    r1_header = read1_record[0].strip() + " " + index_pair + "\n"
                    r2_header = read2_record[0].strip() + " " + index_pair + "\n"
        
                    f1.write(r1_header)
                    f2.write(r2_header)
                else:
                    f1.write(read1_record[i])
                    f2.write(read2_record[i])
        continue

    # Indexes that are not low-quality/unknown and are not matched index-pairs are index-hopped.
    with open("ihop_read1.fq", "a") as f1, open("ihop_read2.fq", "a") as f2:

        num_index_hopping += 1

        # Increment permutations dictionary. The rev comp of sequence 2 should be used.
        key = (index1_sequence, index2_rev)
        permutations_dict[key] += 1

        for i in range(RECORD_LENGTH):
                # append the index-pair to the header of BOTH reads.
                if (i == 0):
                    r1_header = read1_record[0].strip() + " " + index_pair + "\n"
                    r2_header = read2_record[0].strip() + " " + index_pair + "\n"
        
                    f1.write(r1_header)
                    f2.write(r2_header)
                else:
                    f1.write(read1_record[i])
                    f2.write(read2_record[i])

# Create outputfile of the occurrences of each index.
create_output_files(index_count_dictionary, num_matched_indexes, num_unknown_index, num_index_hopping, permutations_dict)

read1_file.close()
read2_file.close()
index1_file.close()
index2_file.close()