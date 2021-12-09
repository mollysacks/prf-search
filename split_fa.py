import argparse
import os
import datetime

def split_fasta(input_fasta, output_folder):
	# Splits fasta with cDNAs into folder for each cDNA
	# Input: fasta with cDNA
	# Output: directory with folders for each cDNA sequence
	
	with open(input_fasta) as f:
		combined_fasta = f.read()
	f.close()
	#initialize file_name variable
	file_name = "init"
	#initialize dictionary with name/ content of each new fasta
	name_content_dict = {}
	cDNAs_list = combined_fasta.split('>')
	for cDNA in cDNAs_list:
		# extract name of cDNA from header
		file_name = cDNA.split(' ')[0]
		if len(file_name) == 0:
			continue

		# put cDNA name and sequence (header included) into dictionary
		name_content_dict[file_name] = '>' + cDNA
	#create new files/ directories
	os.mkdir(output_folder)
	for name, content in name_content_dict.items():
		os.mkdir('{}/{}'.format(output_folder, name))
		with open('{}/{}/{}.fa'.format(output_folder, name, name), 'w+') as f:
			f.write(content)
			f.close()

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-F', default=None, type=str, help='Path to fasta')
	parser.add_argument('-O', default=None, type=str, help='output folder')
	args = parser.parse_args()
	split_fasta(args.F, args.O)