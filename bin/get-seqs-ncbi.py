#!/usr/bin/env python3
import argparse
from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import pickle
import concurrent.futures

def arg_parser():
    parser = argparse.ArgumentParser(description="Get sequences from NCBI")
    parser.add_argument("-a", "--api", help="Gene name", default = os.getenv("NCBI_API_KEY"))
    parser.add_argument("-e", "--email", help="Gene name", default = os.getenv("NCBI_EMAIL"))
    parser.add_argument("-o", "--output", help="Output directory", default = "/mnt/bigdata/linuxhome/millican/ref_db/trait-files")
    return parser.parse_args()

def get_features(record):
    org = ''.join(record.features[0].qualifiers["organism"])
    taxid = ''.join(record.features[0].qualifiers["db_xref"])
    gene = "".join(record.features[2].qualifiers['gene'])
    product = "".join(record.features[1].qualifiers['product'])
    ec = "".join(record.features[1].qualifiers['EC_number'])
    return org, taxid, gene, product, ec

def gene_fetch(gene_key, gene_dict, dir_path):
    ID_list = gene_dict[gene_key]
    seq_dict = {}
    with Entrez.efetch(db="protein", id=ID_list, rettype="gb", retmode="text") as h:
        for record in SeqIO.parse(h, "genbank"):
            seq = record.seq
            name = record.name
            rcid = record.id
            desc = record.description
            org, taxid, gene, product, ec = get_features(record)
            seq_dict[rcid] = {"header": rcid, "name": name, "description": desc, "sequence": seq, "organism": org, "taxid": taxid, "gene": gene, "product": product, "ec": ec}
            with open(f"{dir_path}/{gene_key}/{rcid}.faa", 'a') as f:
                record = SeqRecord(Seq(seq), id=rcid, description='')
                SeqIO.write(record, f, "fasta")
    return seq_dict

def gene_search(gene_name, API):
    gene_dict = {}
    search_term = f"bacteria[Organism] OR archaea[Organism] AND {gene_name}[Gene Name]"
    handle = Entrez.esearch(db="protein", term=search_term, api_key=API, retmax='10000000')
    record = Entrez.read(handle)
    gene_dict[gene] = record["IdList"]
    return gene_dict

def main():
    args = arg_parser()
    Entrez.email = args.email
    with open(f"{args.output}/ncyc_list.pkl", 'rb') as f:
        ncyc_list = pickle.load(f)
    
    gene_dict = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(gene_search, gene, args.api) for gene in ncyc_list]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            gene_dict.update(result)
    
    with open(f"{args.output}/gene_dict.pkl", 'wb') as f:
        pickle.dump(gene_dict, f, pickle.HIGHEST_PROTOCOL)
    
    seq_dict = {}
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(gene_fetch, key, gene_dict, args.output) for key in gene_dict.keys()]
        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            seq_dict.update(result)
    
    with open(f"{args.output}/seq_dict.pkl", 'wb') as f:
        pickle.dump(seq_dict, f, pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':
    main()