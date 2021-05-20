#!/usr/bin/python

# Import all the needed libraries
from Bio import SeqIO
import io
import os, sys
import shutil
import glob

# Define dictionary for sequence ids and sequence data
seqs = {}

path = "/Users/abayomi/Desktop/heather/"
f = "/Users/abayomi/Desktop/heather/ref_genomes.afa"

fa_file = SeqIO.parse(f, "fasta")

items = f.split("/")
seq_filename = items[len(items)-1].rstrip('.afa')

# Store sequence and ids in a dictionary
for seq_data in fa_file:
		seqs[seq_data.id] = seq_data.seq
		ids = seq_data.id


# Information for window analysis
slide_size = 500
step_size = 100

# Make path to new directory for sliding window outputs
sw_output_path = path + "/sliding_window/"

# Make directory if it does not exist
if not os.path.isdir(sw_output_path):
		os.makedirs(sw_output_path)


# make path to new directory for MAFFT re-alignment outputs
mafft_output_path = path + "/mafft/"


# make directory if it does not exist
if not os.path.isdir(mafft_output_path):
		os.makedirs(mafft_output_path)


# make path to new directory for TN93 matrices
tn93_output_path = path + "/tn93/"


# make directory if it does not exist
if not os.path.isdir(tn93_output_path):
		os.makedirs(tn93_output_path)



## try this instead - Art
seqlen = len(list(seqs.values())[0])

for i in range(0, seqlen - slide_size, step_size):

		# TODO: add check for window length

		mafft_infile = "{}{}-sw{}.fa".format(sw_output_path, seq_filename, i)

		# prepare output file
		outfile = open(mafft_infile, 'w')

		# For each batch of sliding window sequences, create a new output file
		for header, seq in iter(seqs.items()):  # if Python 2, you might substitute .iter_items()
				subseq = seq[i:(i+slide_size)]
				outfile.write('>{}\n{}\n'.format(header, subseq))

		outfile.close()

		# MAFFT re-alignment step
		mafft_outfile = "{}{}-sw{}.mafft".format(mafft_output_path, seq_filename, i)
		cmd = 'mafft' + " --thread 10 --op 3 --ep 0.123 --maxiterate 1000 --retree 2 " + mafft_infile + " >" + mafft_outfile
		os.system(cmd)


