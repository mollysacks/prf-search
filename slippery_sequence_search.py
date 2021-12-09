#!/bin/env python

#SBATCH -c 1                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -t 00-00:10         # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p eddy             # Partition to submit to
#SBATCH --mem=10000         # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o ../logs/script_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e ../logs/script_%j.err  # File to which STDERR will be written, %j inserts jobid


import argparse
import re
import os
import shutil

def slippery(string):
	# Returns true if sequence is slippery
	if string[0] != string[1]:
		return False
	if string[1] != string[2]:
		return False
	if string[3] == 'C':
		return False
	if string[3] == 'G':
		return False
	if string[3] != string[4]:
		return False
	if string[4] != string[5]:
		return False
	if string[6] == 'G':
		return False
	return True

def search(reference, length, min_struct_len, upstream_len, space):
	# searches through sequence for slippery sequences
	# outputs list of tuples with (loca, sequence + L bases)
	
	segments = list() #initialize list of segments
	i = 0 # start at first base
	while i <= len(reference) - 7: # go until last 7 bases
		if slippery(reference[i:i+7]): # check if slippery sequence
			if len(reference) - i > length: # if there are > length bases left
				new_seq = reference[i + space :i + space + length] # create new sequence
				if i > upstream_len:
					upstream = reference[i - upstream_len : i]
				else:
					upstream = reference[:i]
				segments.append((i, new_seq, upstream)) # add to segments
			elif len(reference) - (i + space) > min_struct_len: # if there are < length bases left
				if i > upstream_len:
					upstream = reference[i - upstream_len : i]
				else:
					upstream = reference[:i]
				segments.append((i, reference[i + space:len(reference)], upstream)) # add remaining bases in sequence
		i = i + 3 # go to next codon
	return segments

def realign_newlines(seq):
	# write new sequence with properly spaced \n
	new_seq = ''
	i = 0
	while i < len(seq):
		j = i + 79
		if j > len(seq):
			new_seq = new_seq + seq[i:]
		else:
			new_seq = new_seq + seq[i:j] + '\n'
		i += 79
	return new_seq

def find_slippery_sequences(rootdir, length, min_struct_len, upstream_len, space):
	# rootdir contains a directory for each cDNA sequence in input file
	# Each of these folders contains a fasta w/ just that cDNA sequence
	for directory in os.listdir(rootdir):
		for fasta in os.listdir(os.path.join(rootdir, directory)):
			with open(os.path.join(rootdir, directory, fasta)) as f:
				reference = f.read()
			f.close()
			lines = reference.split("\n")
			header = lines[0]
			sequence = ''.join(lines[1:]) 
			segments = search(sequence, length, min_struct_len, upstream_len, space)
			if segments == []:
				shutil.rmtree(os.path.join(rootdir, directory))
			for loc, seq, upstream in segments:
				# realign newlines to improve readability
				seq = realign_newlines(seq)
				os.mkdir(os.path.join(rootdir, directory, '{}-{}'.format(directory, loc)))
				with open(os.path.join(rootdir, directory, '{}-{}'.format(directory, loc),'{}-{}.fa'.format(directory, loc)), "w+") as f:
					f.write(header + "\n")
					f.write(seq)
					f.close()
				with open(os.path.join(rootdir, directory, '{}-{}'.format(directory, loc),'{}-{}.upstream'.format(directory, loc)), "w+") as f:
					f.write(upstream)
					f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-D', default=None, type=str, help='directory w cDNA folders')
	parser.add_argument('-L', default=None, type=int, help='length of sequence after slippery to include in output')
	parser.add_argument('-M', default=18, type=int, help='min length of downstream sequence to include in output')
	parser.add_argument('-U', default=90, type=int, help='max length of upstream segment to search for SD and nascent polypeptide')
	parser.add_argument('-B', default=17, type=int, help='distance between slippery start and potential downstream structure')
	args = parser.parse_args()
	find_slippery_sequences(args.D, args.L, args.M, args.U, args.B)


