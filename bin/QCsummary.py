#!/usr/bin/python3

import pandas as pd
import re

""" Compile statistics on reads and assembly """

def nanoq_parsing(files):
    " Reads and extracts information from nanoq summary files "

    # Creation of the dataframe
    qc_df = pd.DataFrame({"Sample":[], 
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
        qc_df.loc[len(qc_df)] = [sample, reads_number, base_number, N50, longest_read, mean_reads_length, mean_reads_quality]

    return qc_df


if __name__ == "__main__":
    # Load files list
    qc_files = snakemake.input["nanoq"]

    # If only one sample, must converted into list
    if len(qc_files) == 1:
        qc_files = [str(qc_files)]
    
    # Parse files and generate dataframe
    df = nanoq_parsing(qc_files)
    
    # Save file
    df.to_csv(snakemake.output[0], sep="\t")