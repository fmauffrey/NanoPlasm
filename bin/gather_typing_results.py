#!/usr/bin/python3

# Gather results from AMR and mob typing

samples = snakemake.config["samples"]
output = snakemake.output[0]

typing = {}  # master dict to store typing results per sample
amr = {}  # master dict to store amr results per sample

for sample in samples:
    mob_results = open(f"mob_typer/{sample}_mobtyping.tsv", "r").readlines()
    typing[sample] = {}
    for line in mob_results[1:]:
        features = line.split("\t")
        id = features[0].split(" ")[0]
        size, gc, rep, relax, mpf = features[2], round(float(features[3]), 2), features[5], features[7], features[9]
        typing[sample][id] = f"{size}\t{gc}\t{rep}\t{relax}\t{mpf}"

    amr_results = open(f"amrfinder/{sample}_AMR.txt", "r").readlines()

    amr[sample] = {}
    for line in amr_results[1:]:
        features = line.split("\t")
        if features[8] == "AMR":
            id, resist = features[1], features[5]
            if id in amr[sample]:
                if resist not in amr[sample][id]:
                    amr[sample][id].append(resist)
            else:
                amr[sample][id] = [resist]


with open(output, "w") as out:
    out.write("Sample\tSize\tGC\tReplicon\tRelaxase\tMpf\tAMR\n")
    for sample in typing:
        out.write(f"{sample}\t")
        for pos, plasmid in enumerate(typing[sample]):
            if pos == 0:
                try:
                    out.write(f"{typing[sample][plasmid]}\t{amr[sample][plasmid]}\n")
                except KeyError:
                    out.write(f"{typing[sample][plasmid]}\t-\n")
            else:
                try:
                    out.write(f"\t{typing[sample][plasmid]}\t{amr[sample][plasmid]}\n")
                except KeyError:
                    out.write(f"\t{typing[sample][plasmid]}\t-\n")
        out.write("\n")
