# -*- coding: utf-8 -*-
'''
Created Dec 13 2023 
Python 3.10 
Jasper Bellefeuille - belle172@umn.edu 
Repository: Human_Genes/uncharacterized_genes.py 

This script outputs a file, uncharacterized_human_genes.txt, listing uncharacterized human 
protein-coding genes calculated based on data retrieved from HGNC, Wikipedia, and PubMed. 
    The code retrieves from HGNC the regularly updated list of all human protein-coding genes, 
retrieves from Wikipedia every english protein and gene wikidata item, and retrieves from PubMed 
the titles of articles that contain a gene symbol. 

Features of input data expected: 
    1) Wikidata items whose HGNC ID is not in the data obtained from HGNC correspond to withdrawn 
       HGNC records. 
    2) Genes with a wikipedia page have a linked wikidata item. 
    3) Wikidata items with a Wikipedia page means the HGNC ID of that wikidata item is a gene with 
       a wikipedia page. 
    4) A gene symbol being in the title of any pubmed papers means there is research on that gene. 
    5) Once added to the list of gene symbols with pubmed papers, that gene will not go back to 
       being uncharacterized. 

Code for wikidata SQL query generated from https://w.wiki/8VAx, and some code taken from the 
script for writing the wikipedia list of human genes, which can be found at https://w.wiki/8aYR 

To customize the wikipedia user-agent: https://www.mediawiki.org/wiki/Manual:Pywikibot/User-agent 

If you do not have the sparqlwrapper installed, run: pip install sparqlwrapper 
https://rdflib.github.io/sparqlwrapper/ 
''' 

from SPARQLWrapper import SPARQLWrapper, JSON 
from datetime import datetime 
import ftplib 
import io 
import os 
import requests 

# =============================================================================
# Retrieving HGNC file of all human protein-coding genes 
# from https://www.genenames.org/download/statistics-and-files/
# ============================================================================= 
def downloadGeneFile(readFile = 'protein-coding_gene.txt'): # Save 'protein-coding_gene.txt' 
	ftp = ftplib.FTP('ftp.ebi.ac.uk') 
	ftp.login() 
	ftp.cwd('/pub/databases/genenames/new/tsv/locus_groups') 
	with io.open(readFile, 'wb') as data: 
		ftp.retrbinary('RETR protein-coding_gene.txt', data.write) 

downloadGeneFile() 

# =============================================================================
# SQL query retrieving all wikidata gene and protein items matching a HGNC 
# protein-coding gene 
# ============================================================================= 
query = '''SELECT DISTINCT ?gene ?geneLabel ?HGNC_ID ?HGNCsymbol ?protein ?proteinLabel ?wd_gene_item_article_link ?wd_protein_item_article_link
 {
   ?gene wdt:P31 wd:Q7187 .
   ?gene wdt:P703 wd:Q15978631 .
   ?gene wdt:P279 wd:Q20747295 .
   ?gene wdt:P354 ?HGNC_ID .
   ?gene wdt:P353 ?HGNCsymbol .
   ?gene wdt:P688 ?protein .
    OPTIONAL { 
    ?article 	schema:about ?gene ;
                schema:name ?wd_gene_item_article_link ;
 			    schema:isPartOf <https://en.wikipedia.org/> .
    }
    OPTIONAL { 
    ?article 	schema:about ?protein ;
                schema:name ?wd_protein_item_article_link ;
 			    schema:isPartOf <https://en.wikipedia.org/> .
    }
   SERVICE wikibase:label { bd:serviceParam wikibase:language "en" } .
 }''' 

# Function for calling wikipedia API with query 
#   string endpoint_url: the API endpoint with data to retrieve 
def get_results(endpoint_url, query): 

    # example: user_agent = 'CoolBot/0.0 (https://example.org/coolbot/; coolbot@example.org)' 
    user_agent = 'DysoticBot/0.0 (https://en.wikipedia.org/wiki/User:Dysotic)' 

    sparql = SPARQLWrapper(endpoint_url, agent=user_agent)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)
    return sparql.query().convert() 

# returned lists of headers and values for all wikidata human protein and gene items 
results = get_results('https://query.wikidata.org/sparql', query) # includes duplicates 


# =============================================================================
# Create dictionary of genes with their data from HGNC and wikipedia 
# =============================================================================
genes_dict = {} 

HGNC_file = open('protein-coding_gene.txt', encoding='utf-8') 
headers = HGNC_file.readline().strip().split('\t') 

# for each gene symbol from HGNC, add it to the data structures 
for line in HGNC_file: 
    values = line.strip().split('\t') 

    # make dictionary genes_dict where the HGNC ID of each gene is the dictionary's keys 
    ID = values[0].lstrip('HGNC:') 
    genes_dict[ID] = {} 
    genes_dict[ID]['wiki_bool'] = False 
    genes_dict[ID]['wikidatas'] = [] 
    for header, value in enumerate(values): 
        if value != '': 
            genes_dict[ID][headers[header]] = value 

HGNC_file.close() 

# number of wikidata items with no wikipedia page 
wikidata_no_wiki, failed_wikidatas = 0, 0 

# iterate through all wikidata results and append each to genes_dict 
for result in results['results']['bindings']: 
    current_gene = {} 
    for key in result: 
        current_gene[key] = result[key]['value'] 

    try: # check for the wikidata item corresponding to an approved HGNC gene 
        genes_dict[current_gene['HGNC_ID']]['wikidatas'].append(current_gene) 

        try: # if this doesn't throw an error, the wikidata item has a wikipedia page 
            current_gene['wd_gene_item_article_link'] 
            genes_dict[current_gene['HGNC_ID']]['wiki_bool'] = True 

        except KeyError: # the wikidata item doesn't have a gene wikipedia 
            try: # now check that the gene's protein also doesn't have a wikipedia page 
                current_gene['wd_protein_item_article_link'] 
                genes_dict[current_gene['HGNC_ID']]['wiki_bool'] = True 
            except KeyError: 
                wikidata_no_wiki += 1 

    except KeyError:          # if a wikidata's HGNC ID is not in the HGNC list of gene IDs, the 
        failed_wikidatas += 1 # wikidata item is often using a withdrawn gene ID, likely outdated  

no_wikis = {} # make dictionary of genes with no wikipedia page 
for gene in genes_dict: 
    if genes_dict[gene]['wiki_bool'] == False: 
        no_wikis[gene] = genes_dict[gene] 

# =============================================================================
# Retrieve file of characterized genes 
# ============================================================================= 
# Create file name with current year, then check for current year file 
filename = 'characterized_symbols_' + str(datetime.now().year) + '.txt' 

# Read in the list of gene symbols that appeared in at least one pubmed title after previous runs 
# of this code. These genes will be skipped during the actual API retrieval, because the pubmed 
# rate limit means it takes at least 2 hours to retrieve queries for all genes. 
characterized = [] 
newYear = False
try: 
    characterized_file = open(filename) 
    for gene in characterized_file: 
        characterized.append(gene.strip()) 
    characterized_file.close() 
except FileNotFoundError: 
    print('File listing characterized genes not found. The working directory is', os.getcwd(), 
          '\nIf you run the code, it will take several hours to retrieve the pubmed results for', 
          'all genes. Continue? [y/n]') 

    # Get confirmation input before proceeding with full pubmed retrieval 
    if input() == 'y': 
        newYear = True 

# ============================================================================= 
# Retrieve titles of pubmed papers for each gene symbol 
# =============================================================================
no_pubmeds = [] # list of gene HGNC IDs 
no_pubmeds_dict = {} 
partial_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term=' 
prev = datetime.now() 

# For each gene, if its not in the list of characterized symbols from last time, retrieve titles 
for gene in genes_dict: 
    if genes_dict[gene]['symbol'] not in characterized: 
        url = partial_url + genes_dict[gene]['symbol'] + '[title]' 

        # wait for 0.34 seconds to pass, the rate limit of pubmed 
        delta = datetime.now() - prev 
        while delta.seconds < 0.34: 
            delta = datetime.now() - prev 

        # get response of the list of pubmed titles with the gene symbol 
        response = requests.get(url) 
        prev = datetime.now() 

        # Add to list of genes with no pubmed papers, since the symbol is not in any titles 
        if '<IdList>' not in response.text: 
            no_pubmeds.append(genes_dict[gene]['symbol']) 
            no_pubmeds_dict[gene] = genes_dict[gene] 

# TODO: make a timestamp of how long full retrieval takes 
#       try to speed up code by splitting for loop into 
#       saving dictionary of responses[gene] = response.text every 0.34 seconds 
#       for gene in responses: if idlist not in response add to no_pubmeds 


# If running full PubMed retrieval of all genes, create the list of characterized gene symbols 
if newYear == True: 
    characterized_file = open(filename, 'w', encoding='utf-8') 
    for gene in genes_dict: 
        if genes_dict[gene]['symbol'] not in no_pubmeds: 
            print(genes_dict[gene]['symbol'], file=characterized_file) 
    characterized_file.close() 


# =============================================================================
# text processing the data of uncharacterized genes for output file 
# =============================================================================
output_data = {} 

process = ['locus_group', 'locus_type', 'status',                    # same for every gene 
           'agr', 'date_modified', 'gencc', 'location', 'wiki_bool', # duplicate of other data 
           'alias_symbol', 'hgnc_id', 'name', 'symbol', 'wikidatas', # cleanup text and reorder 

           # fields omitted 
           'ena', 'rgd_id', 'vega_id', 'mane_select', 'mgd_id', 'entrez_id', 'iuphar', 'ucsc_id', 
           'date_approved_reserved', 'ccds_id', 'lncipedia'] 

for gene in no_pubmeds_dict: 
    output_data[gene] = {} # initialize gene with garuanteed data items 
    output_data[gene]['symbol'] = genes_dict[gene]['symbol'] 
    output_data[gene]['name'] = genes_dict[gene]['name'] 

    for key in genes_dict[gene]: 
        if key == 'alias_symbol': # change the header name from alias_symbol to aliases 
            output_data[gene]['aliases'] = genes_dict[gene]['alias_symbol'] 

        if key == 'hgnc_id': # process hgnc_id field, which all start with 'HGNC:' 
            output_data[gene]['HGNC_id'] = genes_dict[gene][key].lstrip('HGNC:') 

        # Get the wikidata item IDs instead of full urls 
        if (key == 'wikidatas') and (genes_dict[gene][key] != []): 
            items = '' 
            for item in genes_dict[gene][key]: 
                gene_item = item['gene'].split('entity/')[1] 
                protein_item = item['protein'].split('entity/')[1] 
                items += gene_item + '|' + protein_item 
            output_data[gene][key] = items 

        if key not in process: # copy all remaining keys and corresponding data 
            output_data[gene][key] = genes_dict[gene][key] 

# write output file 
uncharacterized_file = open('uncharacterized_human_genes.txt', 'w', encoding='utf-8') 

# file information 
print('# There are', str(len(genes_dict)), 'human protein-coding genes retrieved from HGNC. This', 
      'file lists the', str(len(no_pubmeds)), 'genes that do not appear in the title of any', 
      'papers on pubmed. This is to serve as an updating estimate of the current number of', 
      'uncharacterized human protein-coding genes.\n# This file was written', 
      str(datetime.now().date()), '(year-month-day).\n', file = uncharacterized_file) 

for gene in output_data: 
    data = '' 
    for key in output_data[gene]: 
        data += key + ': ' + str(output_data[gene][key]) + '\t' 
    print(data + '\n', file = uncharacterized_file) 

uncharacterized_file.close() 


# =============================================================================
# FILE COMPARATOR 
# ============================================================================= 
# TODO: make a function so I can just comment out one line to skip this 
shorter = open(os.getcwd() + '\\uncharacterized_human_genes.txt', encoding='utf-8') 
shorter_matrix = '' 
for line in shorter: 
    shorter_matrix += line 
shorter.close() 
shorter_matrix = shorter_matrix.split('\n') 

shorter_symbols = [] 
for gene in shorter_matrix: 
    if gene != '': 
        shorter_symbols.append(gene.split('\t')[0]) 

longer = open(os.getcwd() + '\\uncharacterized_human_genes_04102024.txt', encoding='utf-8') 
longer_matrix = '' 
for line in longer: 
    longer_matrix += line
longer.close() 
longer_matrix = longer_matrix.split('\n') 

longer_symbols = [] 
for gene in longer_matrix: 
    if gene != '': 
        longer_symbols.append(gene.split('\t')[0])

index = 0 
for i in range(len(longer_symbols)): 
    gene = longer_symbols[index] 
    if gene in shorter_symbols: 
        longer_symbols.remove(gene) 
        shorter_symbols.remove(gene) 
    else: 
        index += 1 

