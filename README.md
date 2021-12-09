#### PRF Search
### By Molly Sacks

This codebase allows users to search through a bacterial genome to find potential -1 progammed ribosomal frameshifting events.
Currently, this can be run locally with small queries (<50 cDNAs). Cloud support is being added soon!

To run this code, you will need working installations of:
1. hmmer
2. RScape
3. python3 and packages:
  a. argparse
  b. os
  c. datetime
  d. pandas
  e. re
  f. shutil
  g. json
  h. RNA (python wrapper for ViennaRNA)
  i. localcider
  j. Bio.Seq (I used Bioconda for this)

Any recent/ stable release should work.

`E_coli_small_2021_12_08-21_53_38.report.tsv` contains output from searching through `small_ecoli.fa`, aligning to database `bacteria.1236.1.genomic.fna`. Intermediate files (including hmmer and RScape/CaCoFold output) are in E_coli_small_2021_12_08-21_53_38.

Here is the command I used to generate this output. I used all the default parameters
`bash search.sh -o E_coli_small_ -q /path/to/small_ecoli.fa -d /path/to/bacteria.1236.1.genomic.fna -r /path/to/rscape_v1.6.1/bin/R-scape -p /path/to/prf-search`
