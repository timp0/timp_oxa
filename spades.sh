#!/bin/bash

#first download the correct version of spades
wget http://spades.bioinf.spbau.ru/release3.9.0/SPAdes-3.9.0-Linux.tar.gz
tar -xzf SPAdes-3.9.0-Linux.tar.gz

#if you are on an AWS instance, you may need to download python as well


#to run (for example, here is sample 79A):
~/SPAdes-3.9.0-Linux/bin/spades.py -1 160108_SSSSS/79A_S1_L001_R1_001.fastq.gz -2 160108_SSSSS/79A_S1_L001_R2_001.fastq.gz -o 79A_Spades

#-1 read 1
#-2 read 2
#-o output directory
