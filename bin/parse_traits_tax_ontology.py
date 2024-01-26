#!/usr/bin/env python3
import pickle
import uuid
import concurrent.futures

tax_hits = "/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax-hits.txt"
ont_file = "/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology.pkl"

with open(ont_file, 'rb') as f:
    ont_dict = pickle.load(f)

trait_tax = {}
uid_list = []

def process_line(line, uid_set):
    query,qheader,qseq,qset,qsetid,taxid,taxname,taxlineage,pident,evalue,tcov,qlen,tlen,target,theader,tseq,tset,tsetid = line.strip().split('\t')
    qid = query.split('_')[0]
    uid = uuid.uuid4()
    while uid in uid_set:
        uid = uuid.uuid4()
    uid_set.add(uid)
    trait_tax[query] = {'uid': uid, 'pgpt': qid, 'pgptid': query, 'descrip': qheader, 'seq': qseq, 'taxid': taxid, 'taxname': taxname, 'taxlineage': taxlineage}

with open(tax_hits, 'r') as f:
    lines = f.readlines()

with concurrent.futures.ThreadPoolExecutor() as executor:
    for i in range(0, len(lines), 1000):
        futures = []
        for line in lines[i:i+1000]:
            futures.append(executor.submit(process_line, line))
        concurrent.futures.wait(futures)

with open("/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax-hits.pkl", 'wb') as f:
    pickle.dump(trait_tax, f, protocol=pickle.HIGHEST_PROTOCOL)

for key in trait_tax:
    if trait_tax[key]['pgpt'] in ont_dict:
        trait_tax[key].update(ont_dict[key])
        print("MATCH!")

with open("/home/glbrc.org/millican/ref_db/trait-files/trait_tax_ontology.pkl", 'wb') as f:
    pickle.dump(trait_tax, f, protocol=pickle.HIGHEST_PROTOCOL)