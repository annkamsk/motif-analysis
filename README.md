Task for Introduction to Computational Biology at MIM UW (http://regulomics.mimuw.edu.pl/wp/2020/05/wbo-assignment-2-zadanie-2/)

`data/protein_fragments.fa` contains protein fragments from 2 tested groups: A and B. We speculate (based on the empirical evidence), that these groups of genes should have different regulatory mechanisms.

The script performs a series of 3 tasks:

1.  Identifies the closest related protein-coding sequences in the E. coli genome for each of the input protein sequences from `data/protein_fragments.fa`. 
The result is a csv file `task1.csv` containing: input sequence id, best matching E. coli gene id and the associated e-value.

E. coli protein-coding sequences are read from: `data/genes_e_coli.fa`.

2. For each of the identified E. coli genes, finds the associated promoter DNA sequence. Motif analysis is done with MEME-suite.

The result is two sets of motif position-specific matrices in a .pfm format.

E. coli promoter sequences are read from: `data/proms_e_coli.fa`.

3. Given the two sets of motifs, selects only the motifs specific to group A or group B by obtaining a number of positions in these promoter sequences that have a log-odds score higher than 0 and using the binomial test. 

The result is a csv file `task3.csv` containing motifs, their associated number of hits in the promoters from group A, and promoters from group B and the associated p-values for enrichment in group A and group B.

# Usage

First, download [NCBI BLAST+ suite](http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec125). The code will run the following command for creating a local BLAST db with this command:
```
makeblastdb -in data/protein_e_coli.fa -parse_seqids -blastdb_version 5 -taxid 1 -title "ecoli" -dbtype prot
```
where `data/protein_e_coli.fa` will contain proteins obtained by translation of `data/genes_e_coli.fa`.

Install dependencies:
```
pip install -r requirements.txt
```

Then run the script with:
```
python main.py PATH_TO_BLAST
```
where PATH_TO_BLAST is a path to NCBI BLAST+ suite.