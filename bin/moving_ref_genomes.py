#!/usr/bin/env python3
import os
import glob
import shutil
from multiprocessing import Pool
import pandas as pd
import sys

def create_acc_from_tax(file):
    acc_list = []
    with open(file, 'r') as f:
        for line in f:
            acc_list.append(line.strip().split(",")[0])
    return set(acc_list)

def create_acc_list(file_list):
    dirpath = "/home/glbrc.org/millican/repos/trait-charmer/workflow/reference/taxonomy"
    acc_list = []
    for file in file_list:
        with open(f"{dirpath}/{file}", 'r') as f:
            for line in f:
                acc_list.append(line.strip().split(",")[0])
    return set(acc_list)

def move_file(acc, taxa):
    try:
        file_path = glob.glob(f"/home/glbrc.org/millican/repos/trait-charmer/workflow/reference/genome/{acc}*fna.gz")[0]
        if ".genomic.fna.gz" in file_path:
            file_path_f = file_path.replace(".genomic.fna.gz", "_genomic.fna.gz")
        else:
            file_path_f = file_path
        shutil.move(file_path, f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/{taxa}")
    except (IndexError, IOError) as e:
        print(f"Error moving file {acc}: {e}")

def move_file_taxa(acc, taxa):
    try:
        file_path = glob.glob(f"/home/glbrc.org/millican/repos/trait-charmer/workflow/reference/genome/{taxa}/{acc}*fna.gz")[0]
        if ".genomic.fna.gz" in file_path:
            file_path_f = file_path.replace(".genomic.fna.gz", "_genomic.fna.gz")
        else:
            file_path_f = file_path
        shutil.move(file_path, f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/{taxa}")
        return
    except (IndexError, IOError) as e:
        print(f"Error moving file {acc}: {e}")
        return

def get_genome_info(list_obj, taxa):
    data_dict = {}
    for item in list_obj:
        try:
            gfile=glob.glob(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/{taxa}/{item}*fna.gz")[0]
            data_dict[item] = {'name': item, 'genome_filename': gfile, 'protein_filename': ""}
        except IndexError as e:
            print(f"Error {e} with {item}")
            continue

def save_genome_info(dict_obj, taxa):
    df = pd.DataFrame.from_dict(dict_obj, orient='index')
    df.to_csv(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_signatures/{taxa}/genome_info_{taxa}.csv")

def main():
    taxa = sys.argv[1]
    if taxa == "archaea":
        arch_file_list = ["archaea_gtdb_taxonomy.csv", "archaea-genbank_ncbi_taxonomy.csv", "archaea-refseq_ncbi_taxonomy.csv"]
        arch_acc_list = create_acc_list(arch_file_list)
        with Pool(32) as p:
            p.starmap(move_file, [(acc, "archaea") for acc in arch_acc_list])
        arch_dict = get_genome_info(arch_acc_list, "archaea")
        save_genome_info(arch_dict, "archaea")
    elif taxa == "bacteria":
        bact_file_list = ["bacteria_gtdb_taxonomy.csv", "bacteria-genbank_ncbi_taxonomy.csv" , "bacteria-refseq_ncbi_taxonomy.csv"]
        bact_acc_list = create_acc_list(bact_file_list)
        with Pool(32) as p:
            p.starmap(move_file, [(acc, "bacteria") for acc in bact_acc_list])
        bact_dict = get_genome_info(bact_acc_list, "bacteria")
        save_genome_info(bact_dict, "bacteria")
    elif taxa == "fungi":
        fung_tax = "/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/fungi_lineage.csv"
        fung_acc_list = create_acc_from_tax(fung_tax)
        with Pool(32) as p:
            p.starmap(move_file_taxa, [(acc, "fungi") for acc in fung_acc_list])
        fung_dict = get_genome_info(fung_acc_list, "fungi")
        save_genome_info(fung_dict, "fungi")
    elif taxa == "invertebrate":
        invert_tax = "/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/invertebrate_lineage.csv"
        inver_acc_list = create_acc_from_tax(invert_tax)
        with Pool(32) as p:
            p.starmap(move_file_taxa, [(acc, "invertebrate") for acc in inver_acc_list])
        inver_dict = get_genome_info(inver_acc_list, "invertebrate")
        save_genome_info(inver_dict, "invertebrate")
    elif taxa == "plant":
        plant_tax = "/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/plant_lineage.csv"
        plant_acc_list = create_acc_from_tax(plant_tax)
        with Pool(32) as p:
            p.starmap(move_file_taxa, [(acc, "plant") for acc in plant_acc_list])
        plant_dict = get_genome_info(plant_acc_list, "plant")
        save_genome_info(plant_dict, "plant")
    elif taxa == "protozoa":
        proto_tax = "/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/protozoa_lineage.csv"
        proto_acc_list = create_acc_from_tax(proto_tax)
        with Pool(32) as p:
            p.starmap(move_file_taxa, [(acc, "protozoa") for acc in proto_acc_list])
        proto_dict = get_genome_info(proto_acc_list, "protozoa")
        save_genome_info(proto_dict, "protozoa")
    elif taxa == "viral":
        viral_tax = "/home/glbrc.org/millican/ref_db/reference_genomes/taxonomy/viral_lineage.csv"
        viral_acc_list = create_acc_from_tax(viral_tax)
        with Pool(32) as p:
            p.starmap(move_file_taxa, [(acc, "viral") for acc in viral_acc_list])
        viral_dict = get_genome_info(viral_acc_list, "viral")
        save_genome_info(viral_dict, "viral")

if __name__ == "__main__":
    for taxa in ["archaea", "bacteria", "fungi", "invertebrate", "plant", "protozoa", "viral"]:
        if not os.path.exists(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/{taxa}"):
            os.makedirs(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_genomes/{taxa}")
        if not os.path.exists(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_signatures/{taxa}"):
            os.makedirs(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_signatures/{taxa}")
            os.makedirs(f"/home/glbrc.org/millican/projects/Inter_BRC/reference_files/reference_signatures/{taxa}/signatures")
    main()


