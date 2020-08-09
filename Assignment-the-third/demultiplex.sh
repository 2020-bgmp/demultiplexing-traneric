#!/bin/bash
#SBATCH --partition=bgmp                ### Partition (like a queue in PBS)
#SBATCH --job-name=demultiplex          ### Job Name
#SBATCH --output=demultiplex.%j.out     ### File in which to store job output
#SBATCH --error=demultiplex.%j.err      ### File in which to store job error messages
#SBATCH --time=0-20:00:00               ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1                       ### Number of nodes needed for the job
#SBATCH --ntasks-per-node=1             ### Number of tasks to be launched per Node
#SBATCH --cpus-per-task=1               ### Number of tasks to be launched
#SBATCH --account=bgmp                  ### Account used for job submission

indexes=/projects/bgmp/shared/2017_sequencing/indexes.txt
read1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz
index1=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz
index2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz
read2=/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz

conda activate bgmp_py37

/usr/bin/time -v python demultiplex.py -c 30 -f $indexes -r1 $read1 -r2 $read2 -i1 $index1 -i2 $index2