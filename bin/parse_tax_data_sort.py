#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import uuid
import pickle
import pandas as pd

fasta_file = "/mnt/bigdata/linuxhome/millican/ref_db/trait-files/tax/trait.faa"
dir_path = "/mnt/bigdata/linuxhome/millican/ref_db/trait-files"
tax_file = "/mnt/bigdata/linuxhome/millican/ref_db/trait-files/tax/trait-tax-hits.txt"

"/home/glbrc.org/millican/ref_db/trait-files/major-functions"




tax_dict = {}
with open(tax_file, 'r') as f:
    for line in f:
        new_ID = str(uuid.uuid4().hex)[:8]
        query,qheader,qseq,qset,qsetid,taxid,taxname,taxlineage,pident,evalue,tcov,qlen,tlen,target,theader,tseq,tset,tsetid = line.strip().split('\t')
        seq = qseq
        header = qheader.split(' ')[1]
        print(header)
        for group in header.split("/"):
            if "|" in group:
                groups = group.split("|")
                for g in groups:
                    with open(f"{dir_path}/sorted/{g}.faa", 'a') as f:
                        record = SeqRecord(Seq(seq), id=new_ID, description="")
                        SeqIO.write(record, f, "fasta")
            else:
                with open(f"{dir_path}/sorted/{group}.faa", 'a') as f:
                    record = SeqRecord(Seq(seq), id=new_ID, description="")
                    SeqIO.write(record, f, "fasta")
        else:
            print(qheader.split('\t'))
            try:
                name = qheader.split('\t')[1]
            except IndexError:
                try:
                    name = qheader.split(' ')[1]
                except IndexError:
                    name = qheader
            with open(f"{dir_path}/sorted/{name}.faa", 'a') as f:
                record = SeqRecord(Seq(seq), id=new_ID, description="")
                SeqIO.write(record, f, "fasta")
        tax_dict[new_ID] = [query,qheader,qseq,taxid,taxname,taxlineage]


with open(f"{dir_path}/tax_dict.csv", 'w') as mp:
        df = pd.DataFrame.from_dict(tax_dict, orient='index')
        df.to_csv(mp, header=True, index=False)

with open(f"{dir_path}/tax_dict.pkl", 'wb') as f:
    pickle.dump(tax_dict, f, pickle.HIGHEST_PROTOCOL)
