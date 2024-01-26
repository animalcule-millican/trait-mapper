#!/usr/bin/env python3
# Date: 2023-12-04
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pickle
import os
import concurrent.futures
## Reconstructing trait-database to integrate new onotology and taxonomic information.
uid_file = "/home/glbrc.org/millican/ref_db/trait-files/uid_pgpt_database.csv"
ont_file = "/home/glbrc.org/millican/repos/trait-mapper/bin/new_pgpt_ontology.txt"
tax_file = "/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax-hits.txt"
ref_dir = "/home/glbrc.org/millican/repos/trait-mapper/reference_database"

def sort_keys(key, uid_dict, tax_dict, ont_dict):
    res_dict = {}
    res_dict = uid_dict
    tax_set = set(tax_dict.keys())
    lookup_val = res_dict[key]['pgptid']
    if lookup_val in tax_set:
        res_dict[key]['taxonomy'] = tax_dict[lookup_val]
        print(f"Lookup val: {lookup_val}")
        print(f"tax_dict[lookup_val]: {tax_dict[lookup_val]}")
        print(f"res_dict[key]['taxonomy']: {res_dict[key]['taxonomy']}")
    #for key2 in tax_dict.keys():
    #        if key2 in set(res_dict[key]['pgptid'].values()):
    #            res_dict[key]['taxonomy'] = tax_dict[key2]
    res_dict[key]['ontology'] = ont_dict[res_dict[key]['ontology']]
    return res_dict

def process_key(key, uid_dict):
    file_lvl1, file_lvl2, file_lvl3, seq, seqid = parse_seqs_ontology(uid_dict, key)
    seq_rec = SeqRecord(Seq(seq), id=seqid, description='')
    return key, [uid_dict[key]['gene'], uid_dict[key]['ko'], uid_dict[key]['ontology'], uid_dict[key]['taxonomy']], [(file_lvl1, seq_rec), (file_lvl2, seq_rec), (file_lvl3, seq_rec)]

def parse_taxonomy(input):
    tax = input.strip().split(";")[1:]
    taxonomy = [item.split("_")[1] for item in tax]
    lineage = ";".join(taxonomy)
    return lineage

def parse_seqs_ontology(input, key):
    ref_dir = "/home/glbrc.org/millican/repos/trait-mapper/reference_database"
    ont = input[key]['ontology'].split("|")
    print(ont)
    file_lvl1 = f"{ref_dir}/lvl1/{ont[0]}.faa"
    file_lvl2 = f"{ref_dir}/lvl2/{ont[0]}-{ont[1]}.faa"
    file_lvl3 = f"{ref_dir}/lvl3/{ont[0]}-{ont[1]}-{ont[2]}.faa"
    print(f"Creating file paths: {file_lvl1}, {file_lvl2}, {file_lvl3}")
    seq = input[key]['seq']
    seqid = key
    return file_lvl1, file_lvl2, file_lvl3, seq, seqid

def build_uid(uid_file = "/home/glbrc.org/millican/ref_db/trait-files/uid_pgpt_database.csv"):
    uid_dict = {}
    with open(uid_file, 'r') as f:
        for i, line in enumerate(f):
            if i == 0:
                continue
            row = line.strip().split(',')
            uid_dict[row[0]] = {"seq": row[1], "desc": row[2], "pgptid": row[3], "gene": row[4], "ko": row[5], "ontology": row[6]}
    return uid_dict

def build_tax(tax_file = "/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax-hits.txt"):
    tax_dict = {}
    with open(tax_file, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            pgptid = row[0].split("-")[0]
            lin = parse_taxonomy(row[7])
            tax_dict[pgptid] = lin
    return tax_dict

def build_ont(ont_file = "/home/glbrc.org/millican/repos/trait-mapper/bin/new_pgpt_ontology.txt"):
    ont_dict = {}
    with open(ont_file, 'r') as f:
        for line in f:
            row = line.strip().split('\t')
            ont_dict[row[0]] = row[1]
    return ont_dict

def main():
    ontology_dict = {}
    database_dict = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=3) as executor:
        f1 = executor.submit(build_uid)
        f2 = executor.submit(build_tax)
        f3 = executor.submit(build_ont)

    uid_dict = f1.result()
    tax_dict = f2.result()
    ont_dict = f3.result()

    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        futures = []
        for key in uid_dict.keys():
            future = {executor.submit(sort_keys, key, uid_dict, tax_dict, ont_dict)}
            futures.append(future)
    print("Sort keys complete")
    # If you need to wait for all functions to complete and gather their results:
    dict_results = [future.result() for future in futures]
    for result in dict_results:
        ontology_dict.update(result)
    
    with open("/home/glbrc.org/millican/repos/trait-mapper/reference_database/ontology-map.pkl", 'wb') as f:
        pickle.dump(ontology_dict, f, protocol = pickle.HIGHEST_PROTOCOL)

    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        futures = {executor.submit(process_key, key, ontology_dict): (key, ontology_dict) for key in ontology_dict.keys()}

    for future in futures:
        results.append(future.result())      
    print("Processing keys complete")
    results = [future.result() for future in new_futures]
    
    with open("/home/glbrc.org/millican/repos/trait-mapper/reference_database/process_key_results.pkl", 'wb') as f:
        pickle.dump(results, f, protocol = pickle.HIGHEST_PROTOCOL)

    for key, database_entry, file_entries in results:
        database_dict[key] = database_entry
        for file_name, seq_rec in file_entries:
            print(f"Writing {file_name}")
            print(f"Number of sequences: {len(seq_rec)}")
            with open(file_name, 'a') as out:
                SeqIO.write(seq_rec, out, 'fasta')

    #for file_name, seq_recs in file_data.items():
    #    print(f"Writing {file_name}")
    #    print(f"Number of sequences: {len(seq_recs)}")
    #    with open(file_name, 'w') as out:
    #        SeqIO.write(seq_recs, out, 'fasta')

    with open("/home/glbrc.org/millican/repos/trait-mapper/reference_database/database-map.pkl", 'wb') as f:
        pickle.dump(database_dict, f, protocol = pickle.HIGHEST_PROTOCOL)
    
if __name__ == "__main__":
    main()