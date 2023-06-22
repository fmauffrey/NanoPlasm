#!/usr/bin/python3

# Isolate putative plasmids contigs from assemblies

from Bio import SeqIO
import os

min_size = snakemake.config["classify_contigs"]["min_size"]
max_size = snakemake.config["classify_contigs"]["max_size"]
output = snakemake.output[0]
assemblies = snakemake.input[0]

# If only one sample, must converted into list
if type(assemblies) != "list":
    assemblies = [assemblies]

# Create folders
if not os.path.exists("Sequences/plasmids"):
    os.mkdir("Sequences/plasmids")
if not os.path.exists("Sequences/chromosomes"):
    os.mkdir("Sequences/chromosomes")

with open(output, "w") as report:
    # Write report first line
    report.write("Sample\tContig ID\tType\tLength\n")
    
    for file in assemblies:
        sample = file.split("/")[1] # sample ID (second term in path)
        assembly = SeqIO.parse(file, "fasta")
        for seq in assembly:
            # If plasmid
            if min_size <= len(seq.seq) <= max_size:
                with open(f"Sequences/plasmids/{sample}_plasmid_{seq.id}.fasta", "w") as out:
                    SeqIO.write([seq], out, "fasta")
                report.write(f"{sample}\t{seq.id}\tPlasmid\t{len(seq.seq)}\n")
            
            # If chromosome
            elif max_size < len(seq.seq):
                with open(f"Sequences/chromosomes/{sample}_chromosome_{seq.id}.fasta", "w") as out:
                    SeqIO.write([seq], out, "fasta")
                report.write(f"{sample}\t{seq.id}\tChromosome\t{len(seq.seq)}\n")