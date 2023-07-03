#!/usr/bin/python3

import pandas as pd

""" Gather results from Resfinder and mobsuite """

def add_resfinder(master_df, res_files):
    " Add resfinder results to the classification file "

    # Creation of dataframe for resfinder results
    res_df = pd.DataFrame({"Sample":[], "Contig ID":[], "ResFinder Resistance genes":[], "ResFinder phenotype":[]})

    # Parsing of all Resfinder files
    for f in res_files:
        sample = f.split("/")[1] # extracted from the path
        results = {} # stores the results per contig
        text = open(f, "r").readlines()[1:] # open and skip header
        for line in text:
            gene = line.split("\t")[0]
            contig = line.split("\t")[5]
            phenotype = line.split("\t")[7].split(", ") # list of phenotypes
            if contig not in results:
                results[contig] = [gene, phenotype] # list of genes in pos 0 and list of phenotype in pos 1
            else:
                results[contig][0] += f", {gene}" # add the new gene to the string in position 0 (list of genes)
                results[contig][1].extend(x for x in phenotype if x not in results[contig][1]) # add phenotype if not already present

        # Looping over results dict to complete dataframe. The label used is the length of the dataframe, creating a simple index.
        for elt in results:
            res_df.loc[len(res_df)] = [sample, elt, results[elt][0], ", ".join(results[elt][1])]

    # Merging the classification dataframe with the resfinder dataframe
    new_df = master_df.merge(res_df, how="outer", on=["Sample", "Contig ID"])

    return new_df

def add_mobsuite(master_df, mob_files):
    " Add mobsuite results to the master file "

    # Creation of dataframe for mobsuite results
    mob_df = pd.DataFrame({"Sample":[], "Contig ID":[], "Mob-suite replicon types":[], "Mob-suite relaxase types":[]})

    # Parsing of all mobsuite files
    for f in mob_files:
        sample = f.split("/")[1].split("_")[0] # extracted from the path
        results = {} # stores the results per contig
        text = open(f, "r").readlines()[1:] # open and skip header
        for line in text:
            contig = line.split("\t")[0]
            replicon = line.split("\t")[5]
            relaxase = line.split("\t")[7]
            results[contig] = [replicon, relaxase]

        # Looping over results dict to complete dataframe. The label used is the length of the dataframe, creating a simple index.
        for elt in results:
            mob_df.loc[len(mob_df)] = [sample, elt, results[elt][0], results[elt][1]]

    # Merging the classification dataframe with the mobsuite dataframe
    new_df = master_df.merge(mob_df, how="outer", on=["Sample", "Contig ID"])

    return new_df


if __name__ == "__main__":

    # Load dataframe and list of files from the Snakefile
    class_file = pd.read_csv(snakemake.input["classification"], sep='\t', header=0)
    class_file = class_file.astype({"Sample": str}) # Convert first column to string in case samples ID are just numbers
    res_files = snakemake.input["resfinder"]
    mob_files = snakemake.input["mobsuite"]

    # If only one sample, must converted into list
    if len(res_files) > 1:
        res_list = res_files
        mob_list = mob_files
    else:
        res_list = [str(res_files)]
        mob_list = [str(mob_files)]
        

    # Add information for each analysis
    infos_resfinder = add_resfinder(class_file, res_list)
    infos_mobsuite = add_mobsuite(infos_resfinder, mob_list)

    # Export the final dataframe
    final_df = infos_mobsuite.sort_values(by=["Sample", "Type", 'Length'])
    final_df.to_csv(snakemake.output[0], sep="\t")