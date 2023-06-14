#!/usr/bin/python3

# Isolate putative plasmids contigs from Medaka assembly

from Bio import SeqIO

circularization = snakemake.config["get_plasmid"]["circularization"]
min_size = snakemake.config["get_plasmid"]["min_size"]
max_size = snakemake.config["get_plasmid"]["max_size"]
assembly = snakemake.input["assembly"]
info = snakemake.input["info"]

plasmids = {}

# Parse info and filter out potential plasmids sequences
for line in open(info, "r").readlines()[1:]:
    contig, length, cov, circ, repeat, *_ = line.split("\t")
    if min_size <= int(length) <= max_size and circ in circularization:
        plasmids[contig] = (length, cov, circ)

# Generate a new fasta file with plasmids sequences
pl_contigs = []
fasta = SeqIO.parse(assembly, "fasta")

for seq in fasta:
    if seq.id in plasmids.keys():
        seq.description = f"length_flye_{plasmids[seq.id][0]}_cov_{plasmids[seq.id][1]}_circ_{plasmids[seq.id][2]}"
        pl_contigs.append(seq)

with open(snakemake.output[0], "w") as out:
    SeqIO.write(pl_contigs, out, "fasta")
