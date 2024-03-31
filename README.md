# Human_Genes
This repository contains the code and resulting output file for the approximate current list of uncharacterized human protein coding genes. The code gets the list of human protein coding genes from a regularly updated data file on the HGNC website, 'protein-coding_gene.txt'. For each gene symbol, if it is not in the title of any papers on PubMed, that gene is put into the output file 'uncharacterized_human_genes.txt'. 

## uncharacterized_genes.py 
Code file - puts gene symbols with no PubMed papers into uncharacterized_human_genes.txt. 

## uncharacterized_human_genes.txt 
Output text file - contains uncharacterized human protein coding gene symbols and any additional data retrieved. 

## protein-coding_gene.txt 
Input text file - protein-coding_gene.txt is retrieved from HGNC and contains the current list of approved human protein coding genes. This file is created by running uncharacterized_genes.py, it does not need to be in the working directory before executing the code. 

## characterized_symbols_2024.txt 
Input text file - contains gene symbols that were found to be in the title of at least one PubMed paper during the first run of uncharacterized_genes.py in 2024. 
If 'characterized_symbols_YYYY.txt' for the current year is not in the current working directory, it is created by running uncharacterized_genes.py. For every gene symbol, the code retrieves the list of PubMed papers with the gene symbol in its title. PubMed API is rate limited to one request every 0.33 seconds, so for >19,250 genes, this takes a theoretical ideal execution time of 2 hours. In practice, the current code takes about 4 hours to do this full retrieval. 
