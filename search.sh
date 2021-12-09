# inputs:
# 	genome to search for frameshifts
# 	db file
# 	Anti SD sequence
# 	number of upstream bases to consider for stemloop/ pseudoknot

#declare default values
upstream=90 #number of upstream bases to evaluate for arrest sequence
downstream=100 #number of downstream bases to evaluate for conserved structure
sd_frame=20 #number of upstream bases to consider for SD-like sequence
structure_start=12 #number of bases after beginning of slippery sequence to start looking for structure

while getopts o:q:d:l:u:s:r:p:i: flag
		do
		        case "${flag}" in
		                o) out=${OPTARG};;
		                q) fa=${OPTARG};;
		                d) db=${OPTARG};;
						l) downstream=${OPTARG};;
						u) upstream=${OPTARG};;
						s) sd_frame=${OPTARG};;
						r) rscape_path=${OPTARG};;
						p) path_to_prf_search=${OPTARG};;
						i) structure_start=${OPTARG};;
		                *) echo "Invalid option: -$flag" ;;
		        esac
		done

 
up=$(($upstream % 3))

if [[ $up -ne 0 ]] ; then 
	echo "Error: Upstream bases (-u) must be divisible by 3."
	exit 1
fi

time_stamp=$(date "+%Y_%m_%d-%H_%M_%S")
root=$out$time_stamp

python3 split_fa.py -F ${fa} -O ${root}
wait
python3 slippery_sequence_search.py -D ${root} -L ${downstream} -U ${upstream} -B ${structure_start}
wait


find_features () {
	for loc in ${cDNA}/*
	# for each slippery sequence
	do
		if [[ -d $loc ]]
        then
        	# find fasta file
        	for seq in $loc/*.fa
        	do
        		# build multiple alignment
        		nhmmer -A ${seq}.sto -o ${seq}.txt ${seq} ${db}
        		
        		# make sure there were hits before proceeding
        		if [ -s ${seq}.sto ]
        		then
        			# continue building alignment
        			hmmbuild -o ${seq}.iter1.hmm.txt ${seq}.iter1.hmm ${seq}.sto 
					nhmmer -A ${seq}.iter1.sto -o ${seq}.iter1.txt ${seq}.iter1.hmm ${db}
					hmmbuild -o ${seq}.iter2.hmm.txt ${seq}.iter2.hmm ${seq}.iter1.sto
					nhmmer -A ${seq}.iter2.sto -o ${seq}.iter2.txt ${seq}.iter2.hmm ${db}
					hmmbuild -o ${seq}.iter3.hmm.txt ${seq}.iter3.hmm ${seq}.iter2.sto
					nhmmer -A ${seq}.iter3.sto -o ${seq}.iter3.txt ${seq}.iter3.hmm ${db}
					
					# fold
					mkdir ${loc}/CaCoFold
					cd ${loc}/CaCoFold
					${rscape_path} --fold ${path_to_prf_search}/${seq}.iter3.sto &> stdout
					pwd
					cp *.fold.sto ${path_to_prf_search}/${loc}
					cd ..
					
					# check for structure and look for other (optional) features
					python3 ${path_to_prf_search}/generate_loc_report.py \
						-S *.fold.sto \
						-F *.fa \
						-B 20 \
						-U *.upstream\
						-R 20 \
						-A ACCUCCU\
						-P 5 \
						-O ${path_to_prf_search}/${root}.report \
						-L ${loc}
					cd ${path_to_prf_search}
				else
        			echo "No hits satisfy inclusion thresholds; no alignment saved for ${loc}"
        			echo "Removing ${loc} from candidate PRF sites"
        			rm -R $loc
        			continue
     			fi
				done
        fi
	done
}

for cDNA in ${root}/*
do 
  	# call find_features for each cDNA
  	find_features ${cDNA} &
done
wait

# reformat output into .tsv
python3 reformat_output.py -O ${root}.report
rm ${root}.report
open ${root}.report.tsv


