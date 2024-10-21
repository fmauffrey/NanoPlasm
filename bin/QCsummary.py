#!/usr/bin/python3

import pandas as pd
import re
from statistics import mean

""" Compile statistics on reads and assembly """

def nanoq_parsing(files):
    " Reads and extracts information from nanoq summary files "

    # Creation of the dataframe
    nanoq_df = pd.DataFrame({"Sample":[], 
                           "Number of reads":[], 
                           "Number of bases":[],
                           "N50 reads length":[], 
                           "Longest read":[],
                           "Mean read length":[],
                           "Mean read quality":[]})

    # Parsing of all nanoq files
    for f in files:
        sample = f.split("/")[1].split("_")[0] # extracted from the path
        text = open(f, "r").read()
        reads_number = re.search(r'(?<=Number of reads:      )[0-9]*', text).group()
        base_number = re.search(r'(?<=Number of bases:      )[0-9]*', text).group()
        N50 = re.search(r'(?<=N50 read length:      )[0-9]*', text).group()
        longest_read = re.search(r'(?<=Longest read:         )[0-9]*', text).group()
        mean_reads_length = re.search(r'(?<=Mean read length:     )[0-9]*', text).group()
        mean_reads_quality = re.search(r'(?<=Mean read quality:    )[0-9]*', text).group()

        # Add info to the dataframe
        nanoq_df.loc[len(nanoq_df)] = [sample, reads_number, base_number, N50, longest_read, mean_reads_length, mean_reads_quality]

    return nanoq_df


def flye_parsing(files, df, min_size, max_size):
    " Reads and extracts information from flye summary files "

    # Creation of the dataframe
    flye_df = pd.DataFrame({"Sample":[], 
                           "Number of plasmids":[], 
                           "Number of circular plasmids":[],
                           "Chromosome coverage":[],
                           "Mean plasmids coverage":[]})

    # Parsing of all flye assembly info files
    for f in files:
        sample = f.split("/")[1].split("_")[0] # extracted from the path
        plasm, circular_plasm, chrom_cov, cov_list = 0, 0, 0, []
        
        text = open(f, "r").readlines()
        for line in text[1:]:
            seq_size = int(line.split("\t")[1])
            cov = int(line.split("\t")[2])
            circ = line.split("\t")[3]
            if min_size <= seq_size <= max_size:
                plasm += 1
                cov_list.append(cov)
                if circ == "Y":
                    circular_plasm += 1
            else:
                chrom_cov = cov
        
        # Check if at least one plasmid
        if not cov_list:
            mean_cov = 0
        else:
            mean_cov = round(mean(cov_list))

        # Add info to the dataframe
        flye_df.loc[len(flye_df)] = [sample, plasm, circular_plasm, chrom_cov, mean_cov]

    # Merge dataframes
    new_df = df.merge(flye_df, how="outer", on=["Sample"])

    return new_df


if __name__ == "__main__":
    # Load files list
    qc_files = snakemake.input["nanoq"]
    flye_files = snakemake.input["flye"]
    plasmid_min_size = snakemake.params["plasmid_min_size"]
    plasmid_max_size = snakemake.params["plasmid_max_size"]

    # If only one sample, must converted into list
    if len(qc_files) == 1:
        qc_files = [str(qc_files)]
        flye_files = [str(flye_files)]
    
    # Parse files and generate dataframe
    df_nanoq = nanoq_parsing(qc_files)
    df = flye_parsing(flye_files, df_nanoq, plasmid_min_size, plasmid_max_size)
    
    # Save file
    df.to_csv(snakemake.output[0], sep="\t")