#!/usr/bin/python3

# Isolate putative plasmids contigs from assemblies

from Bio import SeqIO
import os

min_size = snakemake.config["classify_contigs"]["min_size"]
max_size = snakemake.config["classify_contigs"]["max_size"]
assembly = SeqIO.parse(snakemake.input[0], "fasta")
log = snakemake.output["log"]
hidden_log = snakemake.output["hidden_log"]
sample = snakemake.params["sample"]

# Create folders if not present
if not os.path.exists("06-contigs/plasmids"):
    os.mkdir("06-contigs/plasmids")
if not os.path.exists("06-contigs/chromosomes"):
    os.mkdir("06-contigs/chromosomes")

# Iterate over fasta assembly and classify sequence according to size
chrom = 1 # chromosomes index
plasm = 1 # plasmids index

with open(log, "w") as log, open(hidden_log, "w") as h_log:
    # Write report first line
    log.write("Sample\tPlasmid\tLength\n")
    
    for seq in assembly:
        # If plasmid
        if min_size <= len(seq.seq) <= max_size:
            with open(f"06-contigs/plasmids/{sample}_plasmid#{plasm}.fasta", "w") as out:
                SeqIO.write([seq], out, "fasta")
            log.write(f"{sample}\tPlasmid#{plasm}\t{len(seq.seq)}\n")
            h_log.write(f"{sample}_plasmid#{plasm}.fasta\n")
            plasm += 1
        
        # If chromosome
        elif max_size < len(seq.seq):
            with open(f"06-contigs/chromosomes/{sample}_chrom#{chrom}.fasta", "w") as out:
                SeqIO.write([seq], out, "fasta")
            log.write(f"{sample}\tChromosome#{chrom}\t{len(seq.seq)}\n")
            chrom += 1