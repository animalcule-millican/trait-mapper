#!/usr/bin/env python3
import pandas as pd
import pickle
import sys

def main():
    input_list = sys.argv[1]
    pickle_file = sys.argv[2]
    output_file = sys.argv[3]

    with open(pickle_file, 'rb') as f:
        pgp_map = pickle.load(f)

    #query,qheader,qseq,target,theader,tseq,pident,evalue
    for file in input_list:
        with open(file, 'r') as f:
            for line in f:
                row = line.strip().split('\t')
                data_dict[row[1]] = {'query_id': row[0], 'query_header': row[1], 'annotation_id': row[3], 'query_seq': row[2], 'eval': row[7], 'gene': pgp_map[row[3]]['gene'], 'ko': pgp_map[row[3]]['ko'], 'ontology': pgp_map[row[3]]['ontology']}

    df = pd.DataFrame.from_dict(data_dict, orient='index')
    df.to_csv(output_file, index=False)

if __name__ == '__main__':
    main()