#!/usr/bin/env python3
import os 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import uuid
import pickle

def seq_dict(directory):
    seq_dict = {}
    nid_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".faa"):
            filepath = os.path.join(directory, filename) 
            with open(filepath, 'r') as f:
                for rec in SeqIO.parse(f, 'fasta'):
                    desc = rec.description.split()[1]
                    ids = rec.id.split("-")
                    nid = str(uuid.uuid4().hex)[:8]
                    while nid in nid_list:
                        nid = str(uuid.uuid4().hex)[:8]
                    seq_dict[nid] = {"trait_id": nid, "pgpt_id": ids[0], "gene_name": ids[1], "KO": ids[2], "description": desc, "sequence": rec.seq}
                    nid_list.append(nid)
                    nid = ""
    for key in seq_dict.keys():
        desc = seq_dict[key]["description"].split("/")
        for name in desc:
            with open(f"/home/glbrc.org/millican/ref_db/trait-files/renamed/{name}", 'a') as o:
                record = SeqRecord(Seq(seq_dict[key]["sequence"]), id=seq_dict[key]["trait_id"], description=f"{seq_dict[key]['gene_name']} | {name}")
                SeqIO.write(record, o, 'fasta')
    return seq_dict

def gen_id_seq_dict(seq_dict, id_set):
    new_dict = {}
    for key in seq_dict.keys():
        for ids in id_set:
            new_dict[ids] = {"trait_id": ids, "pgpt_id": seq_dict[key]["pgpt_id"], "gene_name": seq_dict[key]["gene_name"], "KO": seq_dict[key]["KO"], "description": seq_dict[key]["description"], "sequence": seq_dict[key]["sequence"]}
    return new_dict


def rename_seqs(new_dict):
    for key in new_dict.keys():
        desc = new_dict[key]["description"].split("/")
        for name in desc:
            with open(f"/home/glbrc.org/millican/ref_db/trait-files/renamed/{name}", 'a') as o:
                record = SeqRecord(Seq(new_dict[key]["sequence"]), id=new_dict[key]["trait_id"], description=f"{new_dict[key]['gene_name']} | {name}")
                SeqIO.write(record, o, 'fasta')


def generate_id(number_of_ids):
    num_ids = int(number_of_ids)
    id_set = set()
    for i in range(num_ids):
        id_str = str(uuid.uuid4().hex)[:8]
        while id_str in id_set:
            id_str = str(uuid.uuid4().hex)[:8]
        id_set.add(id_str)
    return id_set

def count_lines(input_file):
    with open(input_file, 'r') as f:
        count = 0
        for line in f:
            if line.startswith(">"):
                count += 1
        return count

def count_all_seqs(directory):
    total_seqs = 0
    for filename in os.listdir(directory):
        if filename.endswith(".faa"):
            filepath = os.path.join(directory, filename)
            num_lines = count_lines(filepath)
            total_seqs += num_lines
    return total_seqs

#def rename_seqs(directory, new_dict):
#    for filename in os.listdir(directory):
#        if filename.endswith(".faa"):
#            name = os.path.basename(filename)
#            filepath = os.path.join(directory, filename) 
#            output_file = f"/home/glbrc.org/millican/ref_db/trait-files/renamed/{name}"
#            with open(filepath, 'r') as f:
#                for rec in SeqIO.parse(f, 'fasta'):
#                    ids = rec.id.split("-")[0]
#                    if ids in new_dict["pgpt_id"]:
#                        rec.id = new_dict["trait_id"]
#                        with open(output_file, 'a') as o:
#                            SeqIO.write(rec, o, 'fasta')

if __name__ == "__main__":

    directory = "/home/glbrc.org/millican/ref_db/trait-files/major-functions"
    seqs_count = count_all_seqs(directory)
    seqs_dict = seq_dict(directory)

    id_set = generate_id(seqs_count)
    names_dict = gen_id_seq_dict(seqs_dict, id_set)

    rename_seqs(directory, names_dict)

    with open("/home/glbrc.org/millican/ref_db/trait-files/map-key.pkl", 'wb') as f:
        pickle.dump(names_dict, f)

