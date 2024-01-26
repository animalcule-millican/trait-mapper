#!/usr/bin/env python3
import random
import pandas as pd
import string
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pickle
from concurrent.futures import ProcessPoolExecutor

def create_uid(strings):
    uid = ''.join(random.choice(string.ascii_letters) + ''.join(random.choices(string.ascii_letters + string.digits, k=7))).upper()
    while uid in strings:
        uid = ''.join(random.choice(string.ascii_letters) + ''.join(random.choices(string.ascii_letters + string.digits, k=7))).upper()
    strings.add(uid)
    return uid

def get_items(record, pgp_check):
    split_id = record.id.split("-")
    pgptid = split_id[0]
    # Check if the current pgptid is already been processed.
    if pgptid in pgp_check:
        return None # If the current pgptid is already been processed, return None. This will lead to skipping to next record
    # If there are less than 2 items after split, skip to next record.
    if len(split_id) < 2:
        print(split_id)
        return "LESS THAN TWO"
    # Otherwise, check if the ko is part of the record id. If so, there will be 3 items after split.
    elif len(split_id) > 2:
        pgptid, gene, ko = record.id.split("-")
    # If no KO assigned to record, there should only be two items after split.
    elif len(split_id) == 2:
        pgptid, gene = record.id.split("-")
        ko = 'NA'
    # parse functional ontology from record description
    try:
        ontology = record.description.split(" ")[1]
    except IndexError:
        print("IndexError")
        print(record.description)
        return None
    pgp_check.add(pgptid)
    return pgptid, gene, ko, ontology

def create_sequence_dict(in_file):
    with open(in_file, "r") as handle:
        seq_dict = {}
        pgp_check = set()
        strings = set()
        #seq_dict = {create_uid(strings): {"id": uid, "seq": record.seq, "desc": record.description, 'pgptid': pgptid, 'gene': gene, 'KO': ko, 'ontology': ontology} for record in SeqIO.parse(handle, "fasta") if (items := get_items(record, pgp_check)) is not None and items != "ONLY TWO" and (pgptid, gene, ko, ontology := items) is not None}
        for record in SeqIO.parse(handle, "fasta"):
            pgptid = None
            items = get_items(record, pgp_check)
            if items is None:
                continue
            elif items == "LESS THAN TWO":
                continue
            else:
                pgptid, gene, ko, ontology = items
            if pgptid is None:
                print("Error: pgptid is None")
                print("pgptid should have been reassigned.")
                continue
            uid = create_uid(strings)
            seq_dict[uid] = {"id": uid, "seq": record.seq, "desc": record.description, 'pgptid': pgptid, 'gene': gene, 'KO': ko, 'ontology': ontology}
        return seq_dict

def save_dict(seq_dict, pickle_file):        
    with open(pickle_file, 'wb') as handle:
        pickle.dump(seq_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

def write_dict_2_csv(seq_dict, csv_file):
    df = pd.DataFrame.from_dict(seq_dict, orient='index')
    df.to_csv(csv_file, index=False)

def write_dict_fasta(seq_dict, out_file):
    seq_records = [SeqRecord(Seq(seq_dict[key]["seq"]), id=seq_dict[key]["id"], description = '') for key in seq_dict.keys()]
    with open(out_file, "w") as output_handle:
        SeqIO.write(seq_records, output_handle, "fasta")

def main():
    in_file="/home/glbrc.org/millican/ref_db/pgp/pgpt_trait_seqs.faa"
    out_file="/home/glbrc.org/millican/ref_db/trait-files/pgpt_database.fasta"
    pickle_file="/home/glbrc.org/millican/ref_db/trait-files/pgpt_database_seq_dict.pkl"
    csv_file="/home/glbrc.org/millican/ref_db/trait-files/uid_pgpt_database.csv"
    print("starting create sequence dict now....")
    seq_dict = create_sequence_dict(in_file)
    print("finished create sequence dict now....")
    if not seq_dict:
        print("Error: seq_dict is empty.")
        print("Killing script.")
        return
    else:
        print("seq_dict is not empty.")
        print("Starting executor now...")
        save_dict(seq_dict, pickle_file)
        write_dict_2_csv(seq_dict, csv_file)
        write_dict_fasta(seq_dict, out_file)
        #with ProcessPoolExecutor() as executor:
        #    future1 = executor.submit(save_dict, seq_dict, pickle_file)
        #    future2 = executor.submit(write_dict_2_csv, seq_dict, csv_file)
        #    future3 = executor.submit(write_dict_fasta, seq_dict, out_file)
        #    try:
        #        result1 = future1.result()
        #        result2 = future2.result()
        #        result3 = future3.result()
        #    except Exception as e:
        #        print("An exception occurred: ", e)
        #with ProcessPoolExecutor() as executor:
        #    executor.submit(save_dict, seq_dict, pickle_file)
        #    executor.submit(write_dict_2_csv, seq_dict, csv_file)
        #    executor.submit(write_dict_fasta, seq_dict, out_file)
        print("Finished executor now...")

if __name__ == "__main__":
    main()