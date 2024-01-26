#!/usr/bin/env python3
from tqdm import tqdm
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import multiprocessing
import os
import uuid
import pickle
import tempfile
import shutil

def arg_parser():
    # parse the command line arguments
    parser = argparse.ArgumentParser(description='Parse a large fasta file')
    parser.add_argument('-i', '--input_fasta', help='Input fasta file', default = "/home/glbrc.org/millican/ref_db/pgp/mgPGPT-db.fasta")
    parser.add_argument('-t', '--trait_file', help='Input fasta file', default="/home/glbrc.org/millican/repos/trait-mapper/pgpt-ontology.txt")
    parser.add_argument('-p', '--output_directory', help='Output directory', default="/home/glbrc.org/millican/ref_db/trait-files")
    args = parser.parse_args()
    return args


def get_trait_pickle(trait_file, output_file):
    if os.path.isfile(output_file):
        with open(output_file, 'rb') as f:
            trait_dict = pickle.load(f)
    else:
        trait_dict = {}
        with open(trait_file, 'r') as f, open(output_file, 'wb') as out:
            for line in f:
                row = line.strip().split('\t')
                trait_dict[row[0]] = {'traid_id': row[0], 'gene': row[1], 'ko': row[2], 'class': row[3].lower(), 'category': row[4].lower(), 'subcategory': row[5].lower(), 'group': row[6].lower(), 'function': row[7].lower()}
            pickle.dump(trait_dict, out, protocol=pickle.HIGHEST_PROTOCOL)
    return trait_dict

# Create the chunks of fasta records for multiprocessing
def get_seq_chunks(fasta_file, chunk_size):
    # create a list of fasta records and split it into chunks
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if chunk_size == 0:
        chunk_size = 1
    chunks = [records[i:i+chunk_size] for i in range(0, len(records), chunk_size)]
    return chunks

# Process a chunk of fasta records
def process_chunk(chunk,trait_dict,uid_list):
    # process a chunk of fasta records and write them to a separate output file
    seq_dict = {}
    for rec in tqdm(chunk, desc='Processing chunk'):
        seqid = rec.id.split('-')[0]
        seq_dict[seqid] = {'id': rec.id, 'seq': rec.seq, 'descrip': rec.description}
        for key in seq_dict.keys():
            # get the traits for each sequence
            seq_dict.update(get_traits_seqs(key, seq_dict, trait_dict,uid_list))
    return seq_dict

# Get the traits for each sequence
def get_traits_seqs(seq_key, seq_dict, trait_dict,uid_list):
    new_key = seq_key.split('_')[0]
    for key in trait_dict.keys():
        if key in new_key:
            print("********************")
            print("MATCH!!!!!!!!!!!!!!!")
            print("********************")
            uid = uuid.uuid4()
            while uid in uid_list:
                uid = uuid.uuid4()
            uid_list.append(uid)
            seq_dict[seq_key].update(trait_dict[key])
            seq_dict[seq_key]['uid'] = uid
    return seq_dict

def count_fasta_sequences(fasta_file):
    """
    Count the number of sequences in a fasta file.
    """
    num_sequences = 0
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            num_sequences += 1
    return num_sequences


def worker_write_fasta(seq_key, seq_dict):
    tempfile_list = []
    temp_dir = tempfile.mkdtemp()
    seq = seq_dict[seq_key]['seq']
    ids = seq_dict[seq_key]['uid']
    descrip = f"{seq_dict[seq_key]['gene']};{seq_dict[seq_key]['ko']};{seq_dict[seq_key]['class']};{seq_dict[seq_key]['category']};{seq_dict[seq_key]['subcategory']};{seq_dict[seq_key]['group']};{seq_dict[seq_key]['function']}"
    seq_rec = SeqRecord(Seq(seq), id=ids, description=descrip)
    class_file = f"{temp_dir}/{seq_dict[seq_key]['class']}.faa"
    category_file = f"{temp_dir}/{seq_dict[seq_key]['category']}.faa"
    subcategory_file = f"{temp_dir}/{seq_dict[seq_key]['subcategory']}.faa"
    group_file = f"{temp_dir}/{seq_dict[seq_key]['group']}.faa"
    function_file = f"{temp_dir}/{seq_dict[seq_key]['function']}.faa"
    tempfile_list = [class_file, category_file, subcategory_file, group_file, function_file]
    with open(class_file, 'a') as cf, open(category_file, 'a') as catf, open(subcategory_file, 'a') as subcatf, open(group_file, 'a') as gf, open(function_file, 'a') as ff:
        SeqIO.write(seq_rec, cf, "fasta")
        SeqIO.write(seq_rec, catf, "fasta")
        SeqIO.write(seq_rec, subcatf, "fasta")
        SeqIO.write(seq_rec, gf, "fasta")
        SeqIO.write(seq_rec, ff, "fasta")
    return tempfile_list



# Main function
def main():
    # parse the command line arguments
    args = arg_parser()
    temp_dir = tempfile.mkdtemp()
    seq_name = os.path.splitext(os.path.basename(args.input_fasta))[0]
    trait_name = os.path.splitext(os.path.basename(args.trait_file))[0]
    # count available CPUs
    num_workers = os.cpu_count()
    # set the chunk size based on the number of available CPUs and total number of seqquences
    seq_num = count_fasta_sequences(args.input_fasta)
    chunky_size = round(int(seq_num) / int(num_workers))
    # create pickle file names
    trait_pickle = f"{args.output_directory}/{trait_name}.pkl"
    seq_pickle = f"{args.output_directory}/{seq_name}.pkl"
    # load or create trait pickle, then return trait_dict
    trait_dict = get_trait_pickle(args.trait_file, trait_pickle)
    # if the sequence pickle exists, load it, otherwise create it
    if os.path.isfile(seq_pickle): # loading the pickle file
        with open(seq_pickle, 'rb') as f:
            seq_dict = pickle.load(f)
    else: # creating the pickle file
        # create a shared list using multiprocessing.Manager()
        manager = multiprocessing.Manager()
        uid_list = manager.list()
        chunks = get_seq_chunks(args.input_fasta, chunk_size = chunky_size)
        seq_dict = {}
        with multiprocessing.Pool(num_workers) as pool:
            results = pool.starmap(process_chunk, ((chunk, trait_dict, uid_list) for chunk in chunks))
            for chunk_seq_dict in results:
                seq_dict.update(chunk_seq_dict)

    print(seq_dict)

    #with open(seq_pickle, 'wb') as f:
    #    pickle.dump(seq_dict, f, pickle.HIGHEST_PROTOCOL)

    #with multiprocessing.Pool(num_workers) as pool:
    #    file_list_results = pool.starmap(worker_write_fasta, ((key, seq_dict) for key in seq_dict.keys()))
    #    for file_list in file_list_results:
    #        for file in file_list:
    #            bname = os.path.splitext(os.path.basename(file))[0]
    #            with open(f"{args.output_directory}/{bname}.faa", 'a') as out:
    #                with open(file, 'r') as f:
    #                    shutil.copyfileobj(f, out)


if __name__ == '__main__':
    main()