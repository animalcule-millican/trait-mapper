#!/usr/bin/env python3
import pickle

def parse_lineage(taxonomy):
    lineage = taxonomy.split(';')
    tax_list = []
    for item in lineage:
        if item.startswith('d_'):
            tax_list.append(item)
        # If lineage item starts with "p_" and the domain has been found then add item to list, else return list
        if item.startswith('p_') and len(taxonomy) == 1:
            tax_list.append(item.replace('p_', ''))  # Remove "p_" from phylum
        else:
            return ";".join(taxonomy)
        # If lineage item starts with "c_" and the domain,phylum have been found then add item to list, else return list
        if item.startswith('c_') and len(taxonomy) == 2:
            tax_list.append(item.replace('c_', ''))  # Remove "c_" from class
        else:
            return ";".join(taxonomy)
        # If lineage item starts with "o_" and the domain,phylum,class have been found then add item to list, else return list
        if item.startswith('o_') and len(taxonomy) == 3:
            tax_list.append(item.replace('o_', ''))  # Remove "o_" from order
        else:
            return ";".join(taxonomy)
        # If lineage item starts with "f_" and the domain,phylum,class,order have been found then add item to list, else return list
        if item.startswith('f_') and len(taxonomy) == 4:
            tax_list.append(item.replace('f_', ''))  # Remove "f_" from family
        else:
            return ";".join(taxonomy)
        # If lineage item starts with "g_" and the domain,phylum,class,order,family have been found then add item to list, else return list
        if item.startswith('g_') and len(taxonomy) == 5:
            tax_list.append(item.replace('g_', ''))  # Remove "g_" from genus
        else:
            return ";".join(taxonomy)
        # If lineage item starts with "s_" and the domain,phylum,class,order,family,genus have been found then add item to list, else return list
        if item.startswith('s_') and len(taxonomy) == 6:
            tax_list.append(item.replace('s_', ''))  # Remove "s_" from species
        else:
            return ";".join(taxonomy)
    return ";".join(taxonomy)

input_file = "/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax-hits.txt"
output_file = "/home/glbrc.org/millican/ref_db/trait-files/tax/trait-tax.pkl"
hits_dict = {}

with open(input_file, 'r') as f:
    for line in f:
        hits = line.strip().split('\t')
        lineage = parse_lineage(hits[7])
        if len(lineage.split(';')) == 7:
            species = hits[6]
        ncbi = hits[13]
        taxid = hits[5]
        pgpt_seqid = hits[0].split('-')[0]
        seq = hits[2]
        hits_dict[pgpt_seqid] = {'taxid': taxid, 'ncbi': ncbi, 'lineage': lineage, 'organism': lineage.split(";")[-1], 'seq': seq}

with open(output_file, 'wb') as f:
    pickle.dump(hits_dict, f, protocol=pickle.HIGHEST_PROTOCOL)