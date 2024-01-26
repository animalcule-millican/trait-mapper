#!/usr/bin/env python3
import random
import pandas as pd
import string
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import gzip

in_file="/home/glbrc.org/millican/ref_db/trait-files/mgPGPT-db_Feb2022_ul_dwnld.fasta"
out_file="/home/glbrc.org/millican/ref_db/trait-files/uid_mgPGPT-db.fasta"

with open(in_file, "r") as handle:
    strings = set()
    seq_dict = {}
    for record in SeqIO.parse(handle, "fasta"):
        uid = ''.join(random.choice(string.ascii_letters) + ''.join(random.choices(string.ascii_letters + string.digits, k=7))).upper()
        while uid not in strings:
            uid = ''.join(random.choice(string.ascii_letters) + ''.join(random.choices(string.ascii_letters + string.digits, k=7))).upper()
        strings.add(uid)
        record.id = uid
        with open(out_file, "a") as out_handle:
            SeqIO.write(record, out_handle, "fasta")
        try:
            ko = record.description.split('-')[2]
        except IndexError:
            ko = 'NA'
        seq_dict[uid] = {"id": uid, "seq": record.seq, "desc": record.description, 'pgptid': record.description.split('-')[0], 'gene': record.description.split('-')[1], 'KO': ko, 'ontology': record.description.split(' ')[1]}
    df = pd.DataFrame.from_dict(seq_dict, orient='index')
    df.to_csv("/home/glbrc.org/millican/ref_db/trait-files/uid_mgPGPT-db.csv")
