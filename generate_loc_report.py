import argparse
import RNA
import localcider
from Bio.Seq import Seq
import shutil

def nascent_polypeptide(fa, proline_threshold):
	# Checks upstream bases for nascent polypeptide arrest sequence
	# Input: upstream bases, % proline threshold
	# Output: dictionary with info about nascent chain

	with open(fa) as f:
		seq = f.read()
	dna = Seq(seq)
	amino_acids = dna.translate()
	SeqOb = localcider.sequenceParameters.SequenceParameters(str(amino_acids))
	kappa = SeqOb.get_kappa()
	NCPR = SeqOb.get_NCPR()
	prolines = 0
	for char in amino_acids:
		if char == 'P':
			prolines = prolines + 1
	pct_proline = 100 * prolines/len(amino_acids) 
	positive_residues = False
	if (kappa > .5) ^ (NCPR > 0):
		positive_residues = True
	nascent_info = {}
	if positive_residues or (pct_proline > proline_threshold):
		nascent_info["Arrest Sequence"] = str(amino_acids)
		nascent_info["Kappa"] = kappa
		nascent_info["NCPR"] = NCPR
		nascent_info["Percent Prolines"] = pct_proline
		return nascent_info
	else:
		nascent_info["Arrest Sequence"] = "-"
		nascent_info["Kappa"] = "-"
		nascent_info["NCPR"] = "-"
		nascent_info["Percent Prolines"] = "-"
		return nascent_info

def calc_free_energy(segment, aSD):
	# Finds SD interactions in a given sequence 
	# Input: segment sequence, anti-SD sequence
	# Output: list of dictionaries with SD interactions

	fc = RNA.fold_compound(segment[0] + '&' + aSD)
	# compute MFE and MFE structure
	(mfe_struct, mfe) = fc.mfe()
	SD_interactions = []
	if mfe <= -4.4:
		SD_interactions.append({"SD-like Sequence": segment[0], "SD distance from shift": segment[1],"SD MFE": mfe})
	return SD_interactions

def SD(sequence, aSD, candidate_SD_region_len):
	# Finds the minimum MFE SD interaction in candidate SD region
	# Input: upstream sequeence, anti-SD sequence, candidate region length
	# Output: dictionary with most likely SD interaction

	with open(sequence) as f:
		reference = f.read()
	f.close()
	candidate_SD_region = reference[-candidate_SD_region_len:]
	for i in range(len(candidate_SD_region) - 4):
		loc = len(candidate_SD_region) - i #define starting loc for possible SD motif, this is also distance until slippery site
		SD_interactions = calc_free_energy((candidate_SD_region[i:i+4], loc), aSD)
	if SD_interactions == []:
		return {"SD-likesequence": "-", "SD distance from shift": "-","SD MFE": "-"}
	min_mfe = 0
	min_mfe_dict = {}
	# need to account for ties later
	for SD_dict in SD_interactions:
		if SD_dict["SD MFE"] < min_mfe:
			min_mfe_dict = SD_dict
	return min_mfe_dict

def eval_downstream_structure(sto, ref, bases):
	# Looks for downstream structure
	# Input: CaCoFold fold output, original sequence, number of unpaired bases
	# Output: Whether there is a downstream structure that will stall translation, and CaCoFold structure

	# Read in reference
	with open(ref) as f:
		reference = f.read()
	f.close()

	# Read in CaCoFold output
	with open(sto) as f:
		CaCoFold_out = f.read()
	f.close()
	#evaluate structure
	for line in CaCoFold_out.split('#'):
		if line[0:11] == "=GC SS_cons":
			for string in line.split(" "):
				if len(string) > 7:
					structure = string
	check = structure[0:bases]
	struct = False
	for char in check:
		if char != ":":
			struct = True
			break
	return struct, structure

def generate_output(sto, ref, num_bases, upstream_bases, len_SD_region, aSD, proline_threshold, report_file, loc):
	# Generates dictionary with information about this slippery sequence
	# Only writes to output if there is a downstream structure that will stall translation
	
	struct, structure_string = eval_downstream_structure(sto, ref, num_bases)
	if struct == True:
		# Look for SD-like sequence
		min_mfe_dict = SD(upstream_bases, aSD, len_SD_region)
		# Look for nascent polypeptide arrest sequence
		nascent = nascent_polypeptide(upstream_bases, proline_threshold)
		loc_dict = nascent.copy()
		for key, value in min_mfe_dict.items():
			loc_dict[key] = value
		loc_dict["Slippery Sequence Start"]	= loc.split('-')[2]
		loc_dict["cDNA"] = loc.split('/')[1]
		loc_dict["Structure"] = structure_string[:-1]
		with open(report_file, 'a') as f:
			f.write(str(loc_dict))
			f.write("\n")
		f.close()
		return
	else:
		return


if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument('-S', default=None, type=str, help='Path to CaCoFold output')
	parser.add_argument('-F', default=None, type=str, help='Path to reference sequence')
	parser.add_argument('-B', default=20, type=int, help='# downstream bases to evaluate for structure')
	parser.add_argument('-U', default=None, type=str, help='Path to upstream bases')
	parser.add_argument('-R', default=20, type=int, help='Length of candidate SD region')
	parser.add_argument('-A', default='AUCACCUCCUUU', type=str, help='Anti SD sequence')
	parser.add_argument('-P', default=None, type=int, help='Proline threshold')
	parser.add_argument('-O', default=None, type=str, help='Path to report file')
	parser.add_argument('-L', default=None, type=str, help='Slippery sequence loc')
	args = parser.parse_args()
	generate_output(args.S, args.F, args.B, args.U, args.R, args.A, args.P, args.O, args.L)

