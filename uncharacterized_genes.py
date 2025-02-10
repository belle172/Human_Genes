# -*- coding: utf-8 -*-
'''
Jasper Bellefeuille - jasperbellefeuille@gmail.com 
Repository: Human_Genes/uncharacterized_genes.py 
Created Dec 13 2023 
Python 3.10 

This script outputs a file, uncharacterized_human_genes.txt, listing uncharacterized human 
protein-coding genes calculated based on data retrieved from HGNC, Wikipedia, and PubMed. 
    The code retrieves from HGNC the regularly updated list of all human protein-coding genes, 
retrieves from Wikipedia every english protein and gene wikidata item, and retrieves from PubMed 
the titles of articles that contain a gene symbol. 

Features of input data expected: 
    1) Wikidata items whose HGNC ID is not in the data obtained from HGNC correspond to withdrawn 
           HGNC records. 
    2) Genes with a Wikipedia page have a linked wikidata item. 
    3) Wikidata items with a Wikipedia page means the HGNC ID of that wikidata item is a gene with 
           a Wikipedia page. 
    4) A gene's gene symbol or previous gene symbol being in the title of any PubMed papers means 
           there is research characterizing that gene. 
    5) Once added to the list of gene symbols with PubMed papers, that gene will not go back to 
           being uncharacterized. 

Code for wikidata SQL query generated from https://w.wiki/8VAx, and some code taken from the 
script for writing the Wikipedia list of human genes, which can be found at https://w.wiki/8aYR 

To customize the Wikipedia user-agent: https://www.mediawiki.org/wiki/Manual:Pywikibot/User-agent 

If you do not have the SPARQLWrapper installed, run: pip install sparqlwrapper 
https://rdflib.github.io/sparqlwrapper/ 
''' 

from SPARQLWrapper import SPARQLWrapper, JSON 
from datetime import datetime 
import ftplib 
import io 
import os 
import requests # Allows sending HTTP requests to a web page 

os.chdir('C:\\Users\\jaspe\\GitHub\\Human_Genes') 

# ==================================================================================== 
# Retrieving HGNC file of all human protein-coding genes 
# from https://www.genenames.org/download/statistics-and-files/ 
# Documentation on database: https://ftp.ebi.ac.uk/pub/databases/genenames/README.txt 
# ==================================================================================== 
def downloadGeneFile(readFile = 'protein-coding_gene.txt'): # Save 'protein-coding_gene.txt' 
    ftp = ftplib.FTP('ftp.ebi.ac.uk') 
    ftp.login() 

    # ftp.cwd('/pub/databases/genenames/new/tsv/locus_groups') # original path, deprecated 
    # ftp.cwd('/pub/databases/genenames/hgnc/tsv/locus_groups') # path stated in README, not retrievable 

    ftp.cwd('/pub/databases/genenames/out_of_date_hgnc/tsv/locus_groups') # navigatable path 

    # TODO: retrieve file from Google Storage Bucket of HGNC download files 
    # storage.googleapis.com path? 
    # ftp.cwd('/public-download-files/hgnc/tsv/tsv/locus_groups') 

    with io.open(readFile, 'wb') as data: 
        ftp.retrbinary('RETR protein-coding_gene.txt', data.write) 

downloadGeneFile() 

# =============================================================================
# Call Wikipedia API with SQL query retrieving all wikidata gene and protein 
# items matching a HGNC protein-coding gene 
# ============================================================================= 
def get_results(endpoint_url, query): 
    #   string endpoint_url: the API endpoint with data to retrieve 

    # example: user_agent = 'CoolBot/0.0 (https://example.org/coolbot/; coolbot@example.org)' 
    user_agent = 'DysoticBot/0.1 (https://en.wikipedia.org/wiki/User:Dysotic; https://github.com/belle172/Human_Genes)' 

    sparql = SPARQLWrapper(endpoint_url, agent=user_agent)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON) 
    return sparql.query().convert() 

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

# Returned lists of headers and values for all wikidata human protein and gene items 
results = get_results('https://query.wikidata.org/sparql', query) # includes duplicates 

genes_dict = {} # Create dictionary of genes with their data from HGNC and Wikipedia 

HGNC_file = open('protein-coding_gene.txt', encoding='utf-8') 
headers = HGNC_file.readline().strip().split('\t') 

for line in HGNC_file: # for each gene symbol from HGNC, add it to the data structures 
    values = line.strip().split('\t') 

    # Make dictionary genes_dict where the HGNC ID of each gene is the dictionary's keys 
    ID = values[0].lstrip('HGNC:') 
    genes_dict[ID] = {} 
    genes_dict[ID]['wiki_bool'] = False 
    genes_dict[ID]['wikidatas'] = [] 
    for header, value in enumerate(values): 
        if value != '': genes_dict[ID][headers[header]] = value 

HGNC_file.close() 

# Number of wikidata items with no Wikipedia page 
wikidata_no_wiki, failed_wikidatas = 0, 0 

# Iterate through all wikidata results and append each to genes_dict 
for result in results['results']['bindings']: 
    current_gene = {} 
    for key in result: current_gene[key] = result[key]['value'] 

    try: # Check if the gene has a wikidata item corresponding to an approved HGNC gene 
        genes_dict[current_gene['HGNC_ID']]['wikidatas'].append(current_gene) 

        try: # Check if the gene's wikidata item has a Wikipedia page 
            current_gene['wd_gene_item_article_link'] 
            genes_dict[current_gene['HGNC_ID']]['wiki_bool'] = True 

        except KeyError: # The wikidata item doesn't have a gene Wikipedia, so 
            try: # check that the gene's protein also doesn't have a Wikipedia page 
                current_gene['wd_protein_item_article_link'] 
                genes_dict[current_gene['HGNC_ID']]['wiki_bool'] = True 
            except KeyError: wikidata_no_wiki += 1 

    except KeyError:          # if a wikidata's HGNC ID is not in the HGNC list of gene IDs, the 
        failed_wikidatas += 1 # wikidata item is often using a withdrawn gene ID, likely outdated  

# no_wikis = {} # make dictionary of genes with no Wikipedia page 
# for gene in genes_dict: 
#     if genes_dict[gene]['wiki_bool'] == False: 
#         no_wikis[gene] = genes_dict[gene] 


# =============================================================================
# Retrieve file of characterized genes 
# ============================================================================= 
# Create file name with current year, then check for current year file 
filename = 'characterized_symbols_' + str(datetime.now().year) + '.txt' 

# Read in the list of gene symbols that appeared in at least one PubMed title after previous runs 
# of this code. These genes will be skipped during the actual API retrieval, because the PubMed 
# rate limit means it takes at least 2 hours to retrieve queries for all genes. 
characterized = [] 
newYear = False 
try: 
    characterized_file = open(filename) 
    for gene in characterized_file: characterized.append(gene.strip()) 
    characterized_file.close() 
except FileNotFoundError: 
    print('File listing characterized genes not found. The working directory is', os.getcwd(), 
          '\nIf you run the code, it will take several hours to retrieve the PubMed results for', 
          'all genes. Continue? [y/n]\n') 

    # Set NewYear True if it is January, else get confirmation input before proceeding 
    if datetime.now().month == 1: newYear = True    # with full PubMed retrieval 
    elif input() == 'y': newYear = True 

# =============================================================================
# Retrieve genes with a previous gene symbol 
# =============================================================================
changed = {} # make dictionary 
for gene in genes_dict: 
    if genes_dict[gene]['symbol'] not in characterized: 
        try: 
            genes_dict[gene]['prev_symbol'] 
            changed[gene] = genes_dict[gene] 
        except KeyError: 1 


# ============================================================================= 
# Retrieve titles of PubMed papers for each gene symbol 
# =============================================================================
start = datetime.now() 
prev = datetime.now() 
no_pubmeds = [] # list of gene HGNC IDs 
no_pubmeds_dict = {} 
partial_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?term=' 
headers = {'user_agent': 'DysoticBot/0.1 (https://en.wikipedia.org/wiki/User:Dysotic; https://github.com/belle172/Human_Genes)'} 

# For each gene not in the list of characterized symbols, retrieve titles from PubMed 
for gene in genes_dict: 
    if genes_dict[gene]['symbol'] not in characterized: 
        url = partial_url + genes_dict[gene]['symbol'] + '[title]' 

        # Wait for the rate limit of PubMed API requests 
        delta = datetime.now() - prev 
        while delta.seconds < 0.334: delta = datetime.now() - prev 

        # Get response of the list of PubMed titles with the gene symbol 
        response = requests.get(url, headers=headers).text 
        prev = datetime.now() 

        # Add gene symbols not in the title of any PubMed papers to no_pubmeds 
        if '<IdList>' not in response: 
            try: # Check if the gene has a previous gene symbol 
                url = partial_url + genes_dict[gene]['prev_symbol'] + '[title]' 

                # Check for PubMed papers with the gene's previous gene symbol 
                delta = datetime.now() - prev # Wait for the rate limit of PubMed API requests 
                while delta.seconds < 0.334: delta = datetime.now() - prev 
                response = requests.get(url, headers=headers).text 
                prev = datetime.now() 

                # If the gene's previous symbol is also not in any papers, add to the list of 
                if '<IdList>' not in response: # uncharacterized genes 
                    no_pubmeds.append(genes_dict[gene]['symbol']) 
                    no_pubmeds_dict[gene] = genes_dict[gene] 

            # Add the gene to the list of uncharacterized genes, because it doesn't have a 
            except KeyError: # previous gene symbol 
                no_pubmeds.append(genes_dict[gene]['symbol']) 
                no_pubmeds_dict[gene] = genes_dict[gene] 

time = datetime.now() - start # Timestamp of how long full retrieval takes 

# If running full PubMed retrieval of all genes, create the list of characterized gene symbols 
if newYear == True: 
    characterized_file = open(filename, 'w', encoding='utf-8') 
    for gene in genes_dict: 
        if genes_dict[gene]['symbol'] not in no_pubmeds: 
            print(genes_dict[gene]['symbol'], file=characterized_file) 
    characterized_file.close() 


# =============================================================================
# Text processing the data of uncharacterized genes for output file 
# =============================================================================
process = ['locus_group', 'locus_type', 'status',                    # same for every gene 
           'agr', 'date_modified', 'gencc', 'location', 'wiki_bool', # duplicate of other data 
           'alias_symbol', 'hgnc_id', 'name', 'symbol', 'wikidatas', # cleanup text and reorder 

           # fields omitted 
           'ccds_id', 'date_approved_reserved', 'ena', 'iuphar', 'lncipedia', 'mane_select', 
           'mgd_id', 'rgd_id', 'vega_id'] 

output_data = {} 
for gene in no_pubmeds_dict: 
    output_data[gene] = {} # initialize gene with garuanteed data items 
    output_data[gene]['Symbol'] = genes_dict[gene]['symbol'] 
    output_data[gene]['Name'] = genes_dict[gene]['name'] 

    for key in genes_dict[gene]: 
        if key == 'alias_symbol': # change the header name from alias_symbol to Aliases 
            output_data[gene]['Aliases'] = genes_dict[gene]['alias_symbol'] 

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
      'papers on PubMed. This is to serve as an updating estimate of the current number of', 
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
def get_list(filename): 
    file = open(os.getcwd() + filename, encoding='utf-8') 
    matrix = '' 
    for line in file: matrix += line 
    file.close() 
    matrix = matrix.split('\n') 
    matrix.pop(0) 
    matrix.pop(0) 

    symbols = [] 
    for gene in matrix: 
        if gene != '': symbols.append(gene.split('\t')[0]) 

    return symbols 

shorter_symbols = get_list('\\uncharacterized_human_genes.txt') 
longer_symbols = get_list('\\uncharacterized_human_genes_07102024.txt') 

index = 0 
for i in range(len(longer_symbols)): 
    gene = longer_symbols[index] 
    if gene in shorter_symbols: 
        longer_symbols.remove(gene) 
        shorter_symbols.remove(gene) 
    else: index += 1 


# =============================================================================
# Write file of genes removed from uncharacterized list due to version update 
# =============================================================================
update = False 

if update == True: 
    update_file = open('update_0.1.txt', 'w', encoding='utf-8') 
    
    # file information 
    print('# This file lists the', str(len(longer_symbols)), 'genes excluded from the list of', 
          'uncharacterized human protein-coding genes between code versions 0.0 -> 0.1. This', 
          'version update removes genes whose previous gene symbol is in the title of any papers', 
          'on PubMed. The uncharacterized_human_genes files from February 2024 to July 2024 are', 
          'created from code version 0.0.\n# This file was written', str(datetime.now().date()), 
          '(year-month-day).\n', file = update_file) 

    for gene in longer_symbols: print(gene, file = update_file) 
    update_file.close() 

