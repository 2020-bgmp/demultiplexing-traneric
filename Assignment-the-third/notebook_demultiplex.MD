**Bi622 - Demultiplexing**

Goals (from assignment description): Our goal is to look through a lane of sequencing
generated from the 2017 BGMP cohort’s library preps and determine the level of index
swapping and undetermined index-pairs, before and after quality filtering of index reads.
In order to do this, we must first demultiplex the data. In Assignment the first, we will
develop a strategy to de-multiplex samples to create 48 FASTQ files that contain acceptable
index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped
reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs.

**Part 1 Demultiplexing and Index Swapping - Assignment the First**

Objectives: 7/29/20 - 8/01/20

Write quality_distribution.py to generate quality distribution plots and pseudocode
to de-multiplex the samples.

**The python code, sbatch scripts, and output files are stored in the following directory**:

```
$ /home/etran3/bgmp/bioinfo/Bi622/demultiplex/quality_distribution
```

**The four fastq files are stored in the following directory**:

```
/projects/bgmp/shared/2017_sequencing/
```

**The scripts are called**:
```
quality_distribution.py
read1.sh
read2.sh
read3.sh
read4.sh
```

**Command to run sbatch scripts**:
```
$ sbatch read1.sh
$ sbatch read2.sh
$ sbatch read3.sh
$ sbatch read4.sh
```
**Job ID for Each Run**:
```
Note: Forgot to create output files for read2.

sbatch read1.sh = 12606600
sbatch read2.sh = ?
sbatch read3.sh = 12605709
sbatch read4.sh = 12606601
```
**read1.sh usr/bin/time/ summary**:
```
$ nano read1.12606600.err

Command being timed: "python quality_distribution.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -l 101 -r 1"
	User time (seconds): 6768.45
	System time (seconds): 3.59
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:53:21
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 78736
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 53703
	Voluntary context switches: 1455
	Involuntary context switches: 23876
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
**read3.sh usr/bin/time/ summary**:
```
$ nano read3.12605709.err

Command being timed: "./quality_distribution.py"
	User time (seconds): 2193.58
	System time (seconds): 0.64
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 36:39.71
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 75808
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 62486
	Voluntary context switches: 1060
	Involuntary context switches: 8010
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```
**read4.sh usr/bin/time/ summary**:
```
$ nano read4.12606601.err

Command being timed: "python quality_distribution.py -f /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -l 101 -r 4"
	User time (seconds): 6754.57
	System time (seconds): 4.21
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 1:52:56
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 78736
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 53265
	Voluntary context switches: 848
	Involuntary context switches: 58381
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

**Mean Distribution Plot Output files**:
```
read1.png
read2.png
read3.png
read4.png

Uploaded to: https://github.com/2020-bgmp/demultiplexing-traneric/tree/master/Assignment-the-first

Answers to the questions were submitted as "Answers.md" in the above repository.
```

**Writing Pseudocode for Demultiplex**:
```
Pseudocode written and stored in: https://github.com/2020-bgmp/demultiplexing-traneric/tree/master/Assignment-the-first

As "pseudocode.txt"
```

**Unit Tests**:
```
Stored in: https://github.com/2020-bgmp/demultiplexing-traneric

Inside TEST-input_FASTQ directory as Read1.fq, Read2.fq, Read3.fq, and Read4.fq

And results in TEST-output_FASTQ
```

**Part 2 Demultiplexing and Index Swapping - Assignment the Second**

```
Created pseudocode reviews for:

1. Sam Koehler
2. Susan Collins
3. Michal Babinski

These were reported as issues in their Github respositories.
```

**Part 3 Demultiplexing and Index Swapping - Assignment the Third**

Date: 8/05/20 - 8/07/20

**The python script, sbatch scripts, and output files is stored in the following directory**:
```
$ /home/etran3/bgmp/bioinfo/Bi622/demultiplex/demultiplex_code
```

**The name of the Python and sbatch scripts are called**:
```
demultiplex.py
demultiplex.sh
```

**Commands to run scripts and check status**:
```
$ chmod 755 demultiplex.py
$ sbatch demultiplex.sh
$ squeue -u etran3
```

**demultiplex.sh usr/bin/time summary**:
```
$ nano demultiplex.12673600.err

Command being timed: "python demultiplex.py -c 30 -f /projects/bgmp/shared/2017_sequencing/indexes.txt -r1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz -r2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz -i1 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz -i2 /projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
	User time (seconds): 15038.09
	System time (seconds): 12919.95
	Percent of CPU this job got: 90%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 8:33:19
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 21060
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 727125003
	Voluntary context switches: 6703195
	Involuntary context switches: 2159090
	Swaps: 0
	File system inputs: 0
	File system outputs: 0
	Socket messages sent: 0
	Socket messages received: 0
	Signals delivered: 0
	Page size (bytes): 4096
	Exit status: 0
```

**Note: Program Runtime**:

The runtime for my program took much longer than expected. It took 8 hours and 30 mins to complete. Some things to consider to
optimize runtime efficiency is how I am opening my files in the program. I opened/closed the files whenever I appended something
into it. The program may have run much faster if I'd only open and close each file once.

It may have also helped to run the program on more than one node?

**Output files**:
```
A10_read1.fq
A10_read2.fq
A11_read1.fq
A11_read2.fq
A12_read1.fq
A12_read2.fq
A1_read1.fq
A1_read2.fq
A2_read1.fq
A2_read2.fq
A3_read1.fq
A3_read2.fq
A5_read1.fq
A5_read2.fq
A7_read1.fq
A7_read2.fq
A8_read1.fq
A8_read2.fq
B10_read1.fq
B10_read2.fq
B1_read1.fq
B1_read2.fq
B2_read1.fq
B2_read2.fq
B3_read1.fq
B3_read2.fq
B4_read1.fq
B4_read2.fq
B7_read1.fq
B7_read2.fq
B8_read1.fq
B8_read2.fq
B9_read1.fq
B9_read2.fq
C10_read1.fq
C10_read2.fq
C1_read1.fq
C1_read2.fq
C2_read1.fq
C2_read2.fq
C3_read1.fq
C3_read2.fq
C4_read1.fq
C4_read2.fq
C7_read1.fq
C7_read2.fq
C9_read1.fq
C9_read2.fq
demultiplex.12673600.err
demultiplex.12673600.out
ihop_read1.fq
ihop_read2.fq
index_hop_counts.txt
report.txt
unknown_read1.fq
unknown_read2.fq
```

**User report**:
```
$ nano report.text

Matched Indexes: 330766141 (91.06%)
Unknown or Low Quality Indexes: 31801088 (8.75%)
Index Hopping: 679506 (0.19%)

Index Percentage
GTAGCGTA 2.45
CGATCGAT 1.69
GATCAAGG 1.99
AACAGCGA 2.67
TAGCCATG 3.21
CGGTAATC 1.52
CTCTGGAT 10.54
TACCGGAT 23.01
CTAGCTCA 5.22
CACTTCAC 1.26
GCTACTCT 2.23
ACGATCAG 2.4
TATGGCAC 3.37
TGTTCCGT 4.75
GTCCTAAG 2.66
TCGACAAG 1.16
TCTTCGAC 12.7
ATCATGCG 3.04
ATCGTGGT 2.08
TCGAGAGT 3.54
TCGGATTC 1.39
GATCTTGC 1.1
AGAGTCCA 3.41
AGGATAGC 2.61
```

**Index Hop Counts**:
```
('GTAGCGTA', 'CGATCGAT') 201
('GTAGCGTA', 'GATCAAGG') 243
('GTAGCGTA', 'AACAGCGA') 321
('GTAGCGTA', 'TAGCCATG') 313
('GTAGCGTA', 'CGGTAATC') 136
('GTAGCGTA', 'CTCTGGAT') 886
('GTAGCGTA', 'TACCGGAT') 1803
('GTAGCGTA', 'CTAGCTCA') 960
('GTAGCGTA', 'CACTTCAC') 98
('GTAGCGTA', 'GCTACTCT') 1507
('GTAGCGTA', 'ACGATCAG') 768
('GTAGCGTA', 'TATGGCAC') 351
('GTAGCGTA', 'TGTTCCGT') 550
('GTAGCGTA', 'GTCCTAAG') 400
('GTAGCGTA', 'TCGACAAG') 108
('GTAGCGTA', 'TCTTCGAC') 1241
('GTAGCGTA', 'ATCATGCG') 320
('GTAGCGTA', 'ATCGTGGT') 288
('GTAGCGTA', 'TCGAGAGT') 738
('GTAGCGTA', 'TCGGATTC') 151
('GTAGCGTA', 'GATCTTGC') 146
('GTAGCGTA', 'AGAGTCCA') 506
('GTAGCGTA', 'AGGATAGC') 281
('CGATCGAT', 'GTAGCGTA') 200
('CGATCGAT', 'GATCAAGG') 167
('CGATCGAT', 'AACAGCGA') 160
('CGATCGAT', 'TAGCCATG') 226
('CGATCGAT', 'CGGTAATC') 209
('CGATCGAT', 'CTCTGGAT') 1270
('CGATCGAT', 'TACCGGAT') 2008
('CGATCGAT', 'CTAGCTCA') 442
('CGATCGAT', 'CACTTCAC') 120
('CGATCGAT', 'GCTACTCT') 186
('CGATCGAT', 'ACGATCAG') 268
('CGATCGAT', 'TATGGCAC') 206
('CGATCGAT', 'TGTTCCGT') 494
('CGATCGAT', 'GTCCTAAG') 215
('CGATCGAT', 'TCGACAAG') 102
('CGATCGAT', 'TCTTCGAC') 1135
('CGATCGAT', 'ATCATGCG') 190
('CGATCGAT', 'ATCGTGGT') 160
('CGATCGAT', 'TCGAGAGT') 309
('CGATCGAT', 'TCGGATTC') 114
('CGATCGAT', 'GATCTTGC') 121
('CGATCGAT', 'AGAGTCCA') 269
('CGATCGAT', 'AGGATAGC') 297
('GATCAAGG', 'GTAGCGTA') 256
('GATCAAGG', 'CGATCGAT') 195
('GATCAAGG', 'AACAGCGA') 203
('GATCAAGG', 'TAGCCATG') 375
('GATCAAGG', 'CGGTAATC') 114
('GATCAAGG', 'CTCTGGAT') 719
('GATCAAGG', 'TACCGGAT') 1544
('GATCAAGG', 'CTAGCTCA') 364
('GATCAAGG', 'CACTTCAC') 92
('GATCAAGG', 'GCTACTCT') 233
('GATCAAGG', 'ACGATCAG') 543
('GATCAAGG', 'TATGGCAC') 272
('GATCAAGG', 'TGTTCCGT') 377
('GATCAAGG', 'GTCCTAAG') 409
('GATCAAGG', 'TCGACAAG') 143
('GATCAAGG', 'TCTTCGAC') 17990
('GATCAAGG', 'ATCATGCG') 287
('GATCAAGG', 'ATCGTGGT') 426
('GATCAAGG', 'TCGAGAGT') 296
('GATCAAGG', 'TCGGATTC') 111
('GATCAAGG', 'GATCTTGC') 250
('GATCAAGG', 'AGAGTCCA') 272
('GATCAAGG', 'AGGATAGC') 248
('AACAGCGA', 'GTAGCGTA') 455
('AACAGCGA', 'CGATCGAT') 309
('AACAGCGA', 'GATCAAGG') 323
('AACAGCGA', 'TAGCCATG') 342
('AACAGCGA', 'CGGTAATC') 153
('AACAGCGA', 'CTCTGGAT') 1190
('AACAGCGA', 'TACCGGAT') 2460
('AACAGCGA', 'CTAGCTCA') 636
('AACAGCGA', 'CACTTCAC') 174
('AACAGCGA', 'GCTACTCT') 211
('AACAGCGA', 'ACGATCAG') 539
('AACAGCGA', 'TATGGCAC') 458
('AACAGCGA', 'TGTTCCGT') 581
('AACAGCGA', 'GTCCTAAG') 334
('AACAGCGA', 'TCGACAAG') 128
('AACAGCGA', 'TCTTCGAC') 2063
('AACAGCGA', 'ATCATGCG') 670
('AACAGCGA', 'ATCGTGGT') 426
('AACAGCGA', 'TCGAGAGT') 478
('AACAGCGA', 'TCGGATTC') 358
('AACAGCGA', 'GATCTTGC') 607
('AACAGCGA', 'AGAGTCCA') 577
('AACAGCGA', 'AGGATAGC') 511
('TAGCCATG', 'GTAGCGTA') 311
('TAGCCATG', 'CGATCGAT') 173
('TAGCCATG', 'GATCAAGG') 383
('TAGCCATG', 'AACAGCGA') 250
('TAGCCATG', 'CGGTAATC') 185
('TAGCCATG', 'CTCTGGAT') 1013
('TAGCCATG', 'TACCGGAT') 2710
('TAGCCATG', 'CTAGCTCA') 754
('TAGCCATG', 'CACTTCAC') 242
('TAGCCATG', 'GCTACTCT') 220
('TAGCCATG', 'ACGATCAG') 410
('TAGCCATG', 'TATGGCAC') 489
('TAGCCATG', 'TGTTCCGT') 643
('TAGCCATG', 'GTCCTAAG') 414
('TAGCCATG', 'TCGACAAG') 222
('TAGCCATG', 'TCTTCGAC') 1491
('TAGCCATG', 'ATCATGCG') 407
('TAGCCATG', 'ATCGTGGT') 216
('TAGCCATG', 'TCGAGAGT') 503
('TAGCCATG', 'TCGGATTC') 228
('TAGCCATG', 'GATCTTGC') 186
('TAGCCATG', 'AGAGTCCA') 410
('TAGCCATG', 'AGGATAGC') 309
('CGGTAATC', 'GTAGCGTA') 125
('CGGTAATC', 'CGATCGAT') 213
('CGGTAATC', 'GATCAAGG') 164
('CGGTAATC', 'AACAGCGA') 137
('CGGTAATC', 'TAGCCATG') 327
('CGGTAATC', 'CTCTGGAT') 669
('CGGTAATC', 'TACCGGAT') 7089
('CGGTAATC', 'CTAGCTCA') 393
('CGGTAATC', 'CACTTCAC') 196
('CGGTAATC', 'GCTACTCT') 125
('CGGTAATC', 'ACGATCAG') 181
('CGGTAATC', 'TATGGCAC') 268
('CGGTAATC', 'TGTTCCGT') 321
('CGGTAATC', 'GTCCTAAG') 133
('CGGTAATC', 'TCGACAAG') 76
('CGGTAATC', 'TCTTCGAC') 1190
('CGGTAATC', 'ATCATGCG') 169
('CGGTAATC', 'ATCGTGGT') 355
('CGGTAATC', 'TCGAGAGT') 224
('CGGTAATC', 'TCGGATTC') 153
('CGGTAATC', 'GATCTTGC') 132
('CGGTAATC', 'AGAGTCCA') 194
('CGGTAATC', 'AGGATAGC') 322
('CTCTGGAT', 'GTAGCGTA') 1059
('CTCTGGAT', 'CGATCGAT') 2123
('CTCTGGAT', 'GATCAAGG') 750
('CTCTGGAT', 'AACAGCGA') 985
('CTCTGGAT', 'TAGCCATG') 1174
('CTCTGGAT', 'CGGTAATC') 912
('CTCTGGAT', 'TACCGGAT') 13858
('CTCTGGAT', 'CTAGCTCA') 3196
('CTCTGGAT', 'CACTTCAC') 679
('CTCTGGAT', 'GCTACTCT') 2420
('CTCTGGAT', 'ACGATCAG') 1270
('CTCTGGAT', 'TATGGCAC') 1669
('CTCTGGAT', 'TGTTCCGT') 2455
('CTCTGGAT', 'GTCCTAAG') 1311
('CTCTGGAT', 'TCGACAAG') 532
('CTCTGGAT', 'TCTTCGAC') 5943
('CTCTGGAT', 'ATCATGCG') 1244
('CTCTGGAT', 'ATCGTGGT') 1212
('CTCTGGAT', 'TCGAGAGT') 2234
('CTCTGGAT', 'TCGGATTC') 708
('CTCTGGAT', 'GATCTTGC') 594
('CTCTGGAT', 'AGAGTCCA') 1893
('CTCTGGAT', 'AGGATAGC') 1052
('TACCGGAT', 'GTAGCGTA') 1968
('TACCGGAT', 'CGATCGAT') 2080
('TACCGGAT', 'GATCAAGG') 2042
('TACCGGAT', 'AACAGCGA') 2191
('TACCGGAT', 'TAGCCATG') 4518
('TACCGGAT', 'CGGTAATC') 2715
('TACCGGAT', 'CTCTGGAT') 19845
('TACCGGAT', 'CTAGCTCA') 4640
('TACCGGAT', 'CACTTCAC') 1400
('TACCGGAT', 'GCTACTCT') 1708
('TACCGGAT', 'ACGATCAG') 2604
('TACCGGAT', 'TATGGCAC') 4494
('TACCGGAT', 'TGTTCCGT') 5442
('TACCGGAT', 'GTCCTAAG') 3080
('TACCGGAT', 'TCGACAAG') 1147
('TACCGGAT', 'TCTTCGAC') 13123
('TACCGGAT', 'ATCATGCG') 2441
('TACCGGAT', 'ATCGTGGT') 3138
('TACCGGAT', 'TCGAGAGT') 5340
('TACCGGAT', 'TCGGATTC') 1682
('TACCGGAT', 'GATCTTGC') 1314
('TACCGGAT', 'AGAGTCCA') 2649
('TACCGGAT', 'AGGATAGC') 2395
('CTAGCTCA', 'GTAGCGTA') 1188
('CTAGCTCA', 'CGATCGAT') 569
('CTAGCTCA', 'GATCAAGG') 436
('CTAGCTCA', 'AACAGCGA') 569
('CTAGCTCA', 'TAGCCATG') 831
('CTAGCTCA', 'CGGTAATC') 386
('CTAGCTCA', 'CTCTGGAT') 3794
('CTAGCTCA', 'TACCGGAT') 4263
('CTAGCTCA', 'CACTTCAC') 318
('CTAGCTCA', 'GCTACTCT') 1781
('CTAGCTCA', 'ACGATCAG') 762
('CTAGCTCA', 'TATGGCAC') 886
('CTAGCTCA', 'TGTTCCGT') 1051
('CTAGCTCA', 'GTCCTAAG') 641
('CTAGCTCA', 'TCGACAAG') 14739
('CTAGCTCA', 'TCTTCGAC') 2537
('CTAGCTCA', 'ATCATGCG') 1767
('CTAGCTCA', 'ATCGTGGT') 481
('CTAGCTCA', 'TCGAGAGT') 1109
('CTAGCTCA', 'TCGGATTC') 398
('CTAGCTCA', 'GATCTTGC') 363
('CTAGCTCA', 'AGAGTCCA') 1166
('CTAGCTCA', 'AGGATAGC') 567
('CACTTCAC', 'GTAGCGTA') 124
('CACTTCAC', 'CGATCGAT') 170
('CACTTCAC', 'GATCAAGG') 111
('CACTTCAC', 'AACAGCGA') 118
('CACTTCAC', 'TAGCCATG') 3187
('CACTTCAC', 'CGGTAATC') 219
('CACTTCAC', 'CTCTGGAT') 646
('CACTTCAC', 'TACCGGAT') 1013
('CACTTCAC', 'CTAGCTCA') 347
('CACTTCAC', 'GCTACTCT') 98
('CACTTCAC', 'ACGATCAG') 367
('CACTTCAC', 'TATGGCAC') 349
('CACTTCAC', 'TGTTCCGT') 272
('CACTTCAC', 'GTCCTAAG') 127
('CACTTCAC', 'TCGACAAG') 64
('CACTTCAC', 'TCTTCGAC') 909
('CACTTCAC', 'ATCATGCG') 146
('CACTTCAC', 'ATCGTGGT') 99
('CACTTCAC', 'TCGAGAGT') 169
('CACTTCAC', 'TCGGATTC') 135
('CACTTCAC', 'GATCTTGC') 103
('CACTTCAC', 'AGAGTCCA') 182
('CACTTCAC', 'AGGATAGC') 212
('GCTACTCT', 'GTAGCGTA') 624
('GCTACTCT', 'CGATCGAT') 204
('GCTACTCT', 'GATCAAGG') 217
('GCTACTCT', 'AACAGCGA') 245
('GCTACTCT', 'TAGCCATG') 267
('GCTACTCT', 'CGGTAATC') 132
('GCTACTCT', 'CTCTGGAT') 1110
('GCTACTCT', 'TACCGGAT') 1751
('GCTACTCT', 'CTAGCTCA') 737
('GCTACTCT', 'CACTTCAC') 93
('GCTACTCT', 'ACGATCAG') 345
('GCTACTCT', 'TATGGCAC') 308
('GCTACTCT', 'TGTTCCGT') 629
('GCTACTCT', 'GTCCTAAG') 310
('GCTACTCT', 'TCGACAAG') 95
('GCTACTCT', 'TCTTCGAC') 1057
('GCTACTCT', 'ATCATGCG') 291
('GCTACTCT', 'ATCGTGGT') 210
('GCTACTCT', 'TCGAGAGT') 420
('GCTACTCT', 'TCGGATTC') 135
('GCTACTCT', 'GATCTTGC') 165
('GCTACTCT', 'AGAGTCCA') 445
('GCTACTCT', 'AGGATAGC') 238
('ACGATCAG', 'GTAGCGTA') 259
('ACGATCAG', 'CGATCGAT') 351
('ACGATCAG', 'GATCAAGG') 335
('ACGATCAG', 'AACAGCGA') 269
('ACGATCAG', 'TAGCCATG') 481
('ACGATCAG', 'CGGTAATC') 187
('ACGATCAG', 'CTCTGGAT') 984
('ACGATCAG', 'TACCGGAT') 2157
('ACGATCAG', 'CTAGCTCA') 592
('ACGATCAG', 'CACTTCAC') 184
('ACGATCAG', 'GCTACTCT') 332
('ACGATCAG', 'TATGGCAC') 415
('ACGATCAG', 'TGTTCCGT') 539
('ACGATCAG', 'GTCCTAAG') 492
('ACGATCAG', 'TCGACAAG') 326
('ACGATCAG', 'TCTTCGAC') 1240
('ACGATCAG', 'ATCATGCG') 423
('ACGATCAG', 'ATCGTGGT') 226
('ACGATCAG', 'TCGAGAGT') 502
('ACGATCAG', 'TCGGATTC') 168
('ACGATCAG', 'GATCTTGC') 178
('ACGATCAG', 'AGAGTCCA') 489
('ACGATCAG', 'AGGATAGC') 404
('TATGGCAC', 'GTAGCGTA') 340
('TATGGCAC', 'CGATCGAT') 263
('TATGGCAC', 'GATCAAGG') 327
('TATGGCAC', 'AACAGCGA') 656
('TATGGCAC', 'TAGCCATG') 590
('TATGGCAC', 'CGGTAATC') 340
('TATGGCAC', 'CTCTGGAT') 1532
('TATGGCAC', 'TACCGGAT') 4555
('TATGGCAC', 'CTAGCTCA') 860
('TATGGCAC', 'CACTTCAC') 469
('TATGGCAC', 'GCTACTCT') 313
('TATGGCAC', 'ACGATCAG') 535
('TATGGCAC', 'TGTTCCGT') 88393
('TATGGCAC', 'GTCCTAAG') 732
('TATGGCAC', 'TCGACAAG') 265
('TATGGCAC', 'TCTTCGAC') 4221
('TATGGCAC', 'ATCATGCG') 437
('TATGGCAC', 'ATCGTGGT') 306
('TATGGCAC', 'TCGAGAGT') 920
('TATGGCAC', 'TCGGATTC') 978
('TATGGCAC', 'GATCTTGC') 390
('TATGGCAC', 'AGAGTCCA') 498
('TATGGCAC', 'AGGATAGC') 615
('TGTTCCGT', 'GTAGCGTA') 529
('TGTTCCGT', 'CGATCGAT') 516
('TGTTCCGT', 'GATCAAGG') 382
('TGTTCCGT', 'AACAGCGA') 768
('TGTTCCGT', 'TAGCCATG') 589
('TGTTCCGT', 'CGGTAATC') 266
('TGTTCCGT', 'CTCTGGAT') 1865
('TGTTCCGT', 'TACCGGAT') 4482
('TGTTCCGT', 'CTAGCTCA') 998
('TGTTCCGT', 'CACTTCAC') 262
('TGTTCCGT', 'GCTACTCT') 467
('TGTTCCGT', 'ACGATCAG') 539
('TGTTCCGT', 'TATGGCAC') 85318
('TGTTCCGT', 'GTCCTAAG') 427
('TGTTCCGT', 'TCGACAAG') 231
('TGTTCCGT', 'TCTTCGAC') 2619
('TGTTCCGT', 'ATCATGCG') 652
('TGTTCCGT', 'ATCGTGGT') 473
('TGTTCCGT', 'TCGAGAGT') 1003
('TGTTCCGT', 'TCGGATTC') 863
('TGTTCCGT', 'GATCTTGC') 332
('TGTTCCGT', 'AGAGTCCA') 794
('TGTTCCGT', 'AGGATAGC') 539
('GTCCTAAG', 'GTAGCGTA') 361
('GTCCTAAG', 'CGATCGAT') 238
('GTCCTAAG', 'GATCAAGG') 432
('GTCCTAAG', 'AACAGCGA') 434
('GTCCTAAG', 'TAGCCATG') 433
('GTCCTAAG', 'CGGTAATC') 136
('GTCCTAAG', 'CTCTGGAT') 1093
('GTCCTAAG', 'TACCGGAT') 2255
('GTCCTAAG', 'CTAGCTCA') 567
('GTCCTAAG', 'CACTTCAC') 127
('GTCCTAAG', 'GCTACTCT') 271
('GTCCTAAG', 'ACGATCAG') 426
('GTCCTAAG', 'TATGGCAC') 7014
('GTCCTAAG', 'TGTTCCGT') 935
('GTCCTAAG', 'TCGACAAG') 178
('GTCCTAAG', 'TCTTCGAC') 1251
('GTCCTAAG', 'ATCATGCG') 404
('GTCCTAAG', 'ATCGTGGT') 241
('GTCCTAAG', 'TCGAGAGT') 352
('GTCCTAAG', 'TCGGATTC') 263
('GTCCTAAG', 'GATCTTGC') 216
('GTCCTAAG', 'AGAGTCCA') 368
('GTCCTAAG', 'AGGATAGC') 288
('TCGACAAG', 'GTAGCGTA') 96
('TCGACAAG', 'CGATCGAT') 101
('TCGACAAG', 'GATCAAGG') 141
('TCGACAAG', 'AACAGCGA') 92
('TCGACAAG', 'TAGCCATG') 320
('TCGACAAG', 'CGGTAATC') 52
('TCGACAAG', 'CTCTGGAT') 427
('TCGACAAG', 'TACCGGAT') 909
('TCGACAAG', 'CTAGCTCA') 229
('TCGACAAG', 'CACTTCAC') 30
('TCGACAAG', 'GCTACTCT') 103
('TCGACAAG', 'ACGATCAG') 277
('TCGACAAG', 'TATGGCAC') 179
('TCGACAAG', 'TGTTCCGT') 229
('TCGACAAG', 'GTCCTAAG') 229
('TCGACAAG', 'TCTTCGAC') 756
('TCGACAAG', 'ATCATGCG') 7800
('TCGACAAG', 'ATCGTGGT') 93
('TCGACAAG', 'TCGAGAGT') 347
('TCGACAAG', 'TCGGATTC') 126
('TCGACAAG', 'GATCTTGC') 65
('TCGACAAG', 'AGAGTCCA') 134
('TCGACAAG', 'AGGATAGC') 171
('TCTTCGAC', 'GTAGCGTA') 1260
('TCTTCGAC', 'CGATCGAT') 1174
('TCTTCGAC', 'GATCAAGG') 1876
('TCTTCGAC', 'AACAGCGA') 1032
('TCTTCGAC', 'TAGCCATG') 1642
('TCTTCGAC', 'CGGTAATC') 1148
('TCTTCGAC', 'CTCTGGAT') 5829
('TCTTCGAC', 'TACCGGAT') 10758
('TCTTCGAC', 'CTAGCTCA') 2789
('TCTTCGAC', 'CACTTCAC') 945
('TCTTCGAC', 'GCTACTCT') 1080
('TCTTCGAC', 'ACGATCAG') 1682
('TCTTCGAC', 'TATGGCAC') 2815
('TCTTCGAC', 'TGTTCCGT') 3140
('TCTTCGAC', 'GTCCTAAG') 1309
('TCTTCGAC', 'TCGACAAG') 1026
('TCTTCGAC', 'ATCATGCG') 1408
('TCTTCGAC', 'ATCGTGGT') 4398
('TCTTCGAC', 'TCGAGAGT') 2631
('TCTTCGAC', 'TCGGATTC') 1591
('TCTTCGAC', 'GATCTTGC') 1018
('TCTTCGAC', 'AGAGTCCA') 1521
('TCTTCGAC', 'AGGATAGC') 1930
('ATCATGCG', 'GTAGCGTA') 322
('ATCATGCG', 'CGATCGAT') 253
('ATCATGCG', 'GATCAAGG') 497
('ATCATGCG', 'AACAGCGA') 524
('ATCATGCG', 'TAGCCATG') 864
('ATCATGCG', 'CGGTAATC') 184
('ATCATGCG', 'CTCTGGAT') 1453
('ATCATGCG', 'TACCGGAT') 2374
('ATCATGCG', 'CTAGCTCA') 774
('ATCATGCG', 'CACTTCAC') 166
('ATCATGCG', 'GCTACTCT') 306
('ATCATGCG', 'ACGATCAG') 771
('ATCATGCG', 'TATGGCAC') 765
('ATCATGCG', 'TGTTCCGT') 662
('ATCATGCG', 'GTCCTAAG') 630
('ATCATGCG', 'TCGACAAG') 295
('ATCATGCG', 'TCTTCGAC') 1734
('ATCATGCG', 'ATCGTGGT') 568
('ATCATGCG', 'TCGAGAGT') 459
('ATCATGCG', 'TCGGATTC') 196
('ATCATGCG', 'GATCTTGC') 262
('ATCATGCG', 'AGAGTCCA') 579
('ATCATGCG', 'AGGATAGC') 567
('ATCGTGGT', 'GTAGCGTA') 248
('ATCGTGGT', 'CGATCGAT') 215
('ATCGTGGT', 'GATCAAGG') 598
('ATCGTGGT', 'AACAGCGA') 194
('ATCGTGGT', 'TAGCCATG') 212
('ATCGTGGT', 'CGGTAATC') 165
('ATCGTGGT', 'CTCTGGAT') 1322
('ATCGTGGT', 'TACCGGAT') 2036
('ATCGTGGT', 'CTAGCTCA') 428
('ATCGTGGT', 'CACTTCAC') 84
('ATCGTGGT', 'GCTACTCT') 244
('ATCGTGGT', 'ACGATCAG') 300
('ATCGTGGT', 'TATGGCAC') 277
('ATCGTGGT', 'TGTTCCGT') 434
('ATCGTGGT', 'GTCCTAAG') 174
('ATCGTGGT', 'TCGACAAG') 107
('ATCGTGGT', 'TCTTCGAC') 1346
('ATCGTGGT', 'ATCATGCG') 307
('ATCGTGGT', 'TCGAGAGT') 525
('ATCGTGGT', 'TCGGATTC') 143
('ATCGTGGT', 'GATCTTGC') 127
('ATCGTGGT', 'AGAGTCCA') 352
('ATCGTGGT', 'AGGATAGC') 259
('TCGAGAGT', 'GTAGCGTA') 334
('TCGAGAGT', 'CGATCGAT') 337
('TCGAGAGT', 'GATCAAGG') 306
('TCGAGAGT', 'AACAGCGA') 323
('TCGAGAGT', 'TAGCCATG') 415
('TCGAGAGT', 'CGGTAATC') 246
('TCGAGAGT', 'CTCTGGAT') 1791
('TCGAGAGT', 'TACCGGAT') 3528
('TCGAGAGT', 'CTAGCTCA') 750
('TCGAGAGT', 'CACTTCAC') 155
('TCGAGAGT', 'GCTACTCT') 614
('TCGAGAGT', 'ACGATCAG') 651
('TCGAGAGT', 'TATGGCAC') 466
('TCGAGAGT', 'TGTTCCGT') 1064
('TCGAGAGT', 'GTCCTAAG') 353
('TCGAGAGT', 'TCGACAAG') 276
('TCGAGAGT', 'TCTTCGAC') 1943
('TCGAGAGT', 'ATCATGCG') 354
('TCGAGAGT', 'ATCGTGGT') 489
('TCGAGAGT', 'TCGGATTC') 378
('TCGAGAGT', 'GATCTTGC') 203
('TCGAGAGT', 'AGAGTCCA') 461
('TCGAGAGT', 'AGGATAGC') 484
('TCGGATTC', 'GTAGCGTA') 124
('TCGGATTC', 'CGATCGAT') 243
('TCGGATTC', 'GATCAAGG') 123
('TCGGATTC', 'AACAGCGA') 1025
('TCGGATTC', 'TAGCCATG') 173
('TCGGATTC', 'CGGTAATC') 187
('TCGGATTC', 'CTCTGGAT') 574
('TCGGATTC', 'TACCGGAT') 1114
('TCGGATTC', 'CTAGCTCA') 366
('TCGGATTC', 'CACTTCAC') 108
('TCGGATTC', 'GCTACTCT') 116
('TCGGATTC', 'ACGATCAG') 204
('TCGGATTC', 'TATGGCAC') 407
('TCGGATTC', 'TGTTCCGT') 295
('TCGGATTC', 'GTCCTAAG') 151
('TCGGATTC', 'TCGACAAG') 84
('TCGGATTC', 'TCTTCGAC') 1246
('TCGGATTC', 'ATCATGCG') 159
('TCGGATTC', 'ATCGTGGT') 114
('TCGGATTC', 'TCGAGAGT') 282
('TCGGATTC', 'GATCTTGC') 289
('TCGGATTC', 'AGAGTCCA') 207
('TCGGATTC', 'AGGATAGC') 290
('GATCTTGC', 'GTAGCGTA') 151
('GATCTTGC', 'CGATCGAT') 195
('GATCTTGC', 'GATCAAGG') 175
('GATCTTGC', 'AACAGCGA') 167
('GATCTTGC', 'TAGCCATG') 126
('GATCTTGC', 'CGGTAATC') 123
('GATCTTGC', 'CTCTGGAT') 370
('GATCTTGC', 'TACCGGAT') 754
('GATCTTGC', 'CTAGCTCA') 228
('GATCTTGC', 'CACTTCAC') 102
('GATCTTGC', 'GCTACTCT') 125
('GATCTTGC', 'ACGATCAG') 171
('GATCTTGC', 'TATGGCAC') 234
('GATCTTGC', 'TGTTCCGT') 191
('GATCTTGC', 'GTCCTAAG') 163
('GATCTTGC', 'TCGACAAG') 46
('GATCTTGC', 'TCTTCGAC') 711
('GATCTTGC', 'ATCATGCG') 227
('GATCTTGC', 'ATCGTGGT') 95
('GATCTTGC', 'TCGAGAGT') 160
('GATCTTGC', 'TCGGATTC') 81
('GATCTTGC', 'AGAGTCCA') 183
('GATCTTGC', 'AGGATAGC') 359
('AGAGTCCA', 'GTAGCGTA') 1027
('AGAGTCCA', 'CGATCGAT') 274
('AGAGTCCA', 'GATCAAGG') 235
('AGAGTCCA', 'AACAGCGA') 462
('AGAGTCCA', 'TAGCCATG') 376
('AGAGTCCA', 'CGGTAATC') 214
('AGAGTCCA', 'CTCTGGAT') 1113
('AGAGTCCA', 'TACCGGAT') 2480
('AGAGTCCA', 'CTAGCTCA') 1215
('AGAGTCCA', 'CACTTCAC') 188
('AGAGTCCA', 'GCTACTCT') 325
('AGAGTCCA', 'ACGATCAG') 627
('AGAGTCCA', 'TATGGCAC') 495
('AGAGTCCA', 'TGTTCCGT') 762
('AGAGTCCA', 'GTCCTAAG') 363
('AGAGTCCA', 'TCGACAAG') 135
('AGAGTCCA', 'TCTTCGAC') 1429
('AGAGTCCA', 'ATCATGCG') 440
('AGAGTCCA', 'ATCGTGGT') 316
('AGAGTCCA', 'TCGAGAGT') 1014
('AGAGTCCA', 'TCGGATTC') 206
('AGAGTCCA', 'GATCTTGC') 174
('AGAGTCCA', 'AGGATAGC') 556
('AGGATAGC', 'GTAGCGTA') 953
('AGGATAGC', 'CGATCGAT') 265
('AGGATAGC', 'GATCAAGG') 208
('AGGATAGC', 'AACAGCGA') 337
('AGGATAGC', 'TAGCCATG') 350
('AGGATAGC', 'CGGTAATC') 308
('AGGATAGC', 'CTCTGGAT') 941
('AGGATAGC', 'TACCGGAT') 1992
('AGGATAGC', 'CTAGCTCA') 605
('AGGATAGC', 'CACTTCAC') 356
('AGGATAGC', 'GCTACTCT') 374
('AGGATAGC', 'ACGATCAG') 512
('AGGATAGC', 'TATGGCAC') 612
('AGGATAGC', 'TGTTCCGT') 574
('AGGATAGC', 'GTCCTAAG') 283
('AGGATAGC', 'TCGACAAG') 127
('AGGATAGC', 'TCTTCGAC') 1915
('AGGATAGC', 'ATCATGCG') 488
('AGGATAGC', 'ATCGTGGT') 281
('AGGATAGC', 'TCGAGAGT') 519
('AGGATAGC', 'TCGGATTC') 328
('AGGATAGC', 'GATCTTGC') 522
('AGGATAGC', 'AGAGTCCA') 579
```
