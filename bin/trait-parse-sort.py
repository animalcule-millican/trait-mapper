#!/usr/bin/env python3
import os 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import uuid
import pandas as pd
import pickle

directory = "/home/glbrc.org/millican/ref_db/trait-files/major-functions"

def process_file(filepath, id_list):
    seq_dict = {}
    with open(filepath, 'r') as f:
        for rec in SeqIO.parse(f, 'fasta'):
            desc = rec.description.split()[1]
            ids = rec.id.split("-")
            nid = str(uuid.uuid4().hex)[:8]
            while nid in id_list:
                nid = str(uuid.uuid4().hex)[:8]
            seq_dict[nid] = {"trait_id": nid, "pgpt_id": ids[0], "gene_name": ids[1], "KO": ids[2], "description": desc, "sequence": rec.seq}
            print(f"for {ids[0]}:")
            print(ids[1])
            print(desc)
            print(" ")
            id_list.append(nid)
            nid = ""
    return seq_dict, id_list

def process_key(key, seq_dict):
    desc = seq_dict[key]["description"].split("/")
    for name in desc:
        with lock:
            with open(f"/home/glbrc.org/millican/ref_db/trait-files/renamed/{name}", 'a') as o:
                record = SeqRecord(Seq(seq_dict[key]["sequence"]), id=seq_dict[key]["trait_id"], description=f"{seq_dict[key]['gene_name']} | {name}")
                SeqIO.write(record, o, 'fasta')

def save_map_key(map_dict):
    with open("/home/glbrc.org/millican/ref_db/trait-files/tax-2-id-map-key.csv", 'a') as mp:
        df = pd.DataFrame.from_dict(map_dict, orient='index')
        df.to_csv(mp, header=True, index=False)

    with open("/home/glbrc.org/millican/ref_db/trait-files/tax-2-id-map-key.pkl", 'wb') as f:
        pickle.dump(map_dict, f)
    
if __name__ == '__main__':
    map_dict = {}
    id_list = []
    out_dict = {}
    out_list = []
    filepaths = [os.path.join(directory, filename) for filename in os.listdir(directory) if filename.endswith(".faa")]
    for files in filepaths:
        out_dict, out_list = process_file(files, id_list)
        map_dict.update(out_dict)
        id_list.extend(out_list)
    
    keys = list(map_dict.keys())
    for key in keys:
        process_key(key, map_dict)
    
    save_map_key(map_dict)