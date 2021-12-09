import argparse
import pandas as pd
import json

def reformat_output(output):
	# Convert report into TSV
	# Input: file with list of dictionaries
	# Output: TSV file with data about each potential PRF site

	dict_list = []
	with open(output, 'r') as f:
		for line in f:
			line = line.replace("'", '"')
			line_dict = json.loads(line)
			dict_list.append(line_dict)
	f.close()
	df = pd.DataFrame(dict_list)
	first_column = df.pop('cDNA')
	second_column = df.pop('Slippery Sequence Start')
	df.insert(0, 'Slippery Sequence Start', second_column)
	df.insert(0, 'cDNA', first_column)
	df.to_csv(f'{output}.tsv', sep="\t", index=False)
	return

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-O', default=None, type=str, help='Path to report file')
	args = parser.parse_args()
	reformat_output(args.O)
