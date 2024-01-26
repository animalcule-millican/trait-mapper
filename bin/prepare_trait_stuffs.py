import pickle
import pandas as pd
import uuid

seqdb = "/home/glbrc.org/millican/ref_db/trait-files/pickles/seq_database.pkl"
tont = "/home/glbrc.org/millican/ref_db/trait-files/pickles/trait-ontology.pkl"


headers = "/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_headers.txt"

head_dict = {}

with open(headers, 'r') as f:
    for line in f:
        if line.startswith(">PGPT"):
            head,ont = line.strip().split(" ")
            head = head.replace(">", "")
        else:
            header = line.strip().split(">")[-1]
            head,ont = header.strip().split(" ")
        traitid,gene,ko = head.split("-")
        pgpt = traitid.split("_")[0]
        head_dict[traitid] = {'trait_id': traitid, 'pgpt': pgpt, 'gene': gene, 'ko': ko, 'ontology': ont}

with open("/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_headers.pkl", 'wb') as f:
    pickle.dump(head_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

df = pd.DataFrame.from_dict(head_dict, orient='index')
df.to_csv("/home/glbrc.org/millican/ref_db/trait-files/pgpt-ontology_headers.txt", sep="\t")




with open(seqdb, 'rb') as f:
    seqdb = pickle.load(f)

with open(tont, 'rb') as f:
    tont = pickle.load(f)

nid_list = []
new_dict = {}
for key in seqdb.keys():
    if seqdb[key]['pgpt'] in tont.keys():
        nid = str(str(uuid.uuid4().hex)[:12])
        while nid in nid_list:
            nid = str(str(uuid.uuid4().hex)[:12])
        nid_list.append(nid)
        new_dict[nid] = {'id': nid, 'pgptID': key, 'traid_id': tont[seqdb[key]["pgpt"]]["trait_id"], 'pgpt': seqdb[key]['pgpt'], 'seq': seqdb[key]['seq'],'gene': tont[seqdb[key]["gene"]]["trait_id"], 'ko': tont[seqdb[key]["pgpt"]]["ko"], 'class': tont[seqdb[key]["pgpt"]]["class"], 'category': tont[seqdb[key]["pgpt"]]["category"], 'subcategory': tont[seqdb[key]["pgpt"]]["subcategory"], 'group': tont[seqdb[key]["pgpt"]]["group"], 'function': tont[seqdb[key]["pgpt"]]["function"]}
        nid_list.append(nid)