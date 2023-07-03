#!/usr/bin/python3

""" Parse minimap2 alignment and karga results for linking karga profil to contigs """

def parse(karga, alignment, f_out):
    # Open files
    k_file = open(karga, "r").readlines()
    a_file = open(alignment, "r").readlines()

    # Output object as lists of resistances per contigs
    output = {}

    # Parse alignment file and save results per read
    reads_contigs = {}
    for line in a_file:
        read = line.split("\t")[0]
        contig = line.split("\t")[5]
        reads_contigs[read] = contig

    # Parse KARGA output, extract resistance per reads and assign them to contigs
    for line in k_file[1:]:
        read = line.split(",")[0].split(" ")[0]
        read_ID = read.replace("@", "")
        try:
            contig = reads_contigs[read_ID]
        except KeyError:
            pass

        resistance = line.split(",")[2]
        try:
            res_cat = resistance.split("|")[1]
            res_gene = resistance.split("|")[4].replace("\n", "")
        except IndexError:
            pass

        # Create nested dictionary checking if value already present
        output.setdefault(contig, {}).setdefault(res_cat, []).append(res_gene)

    # Save results in tsv format
    with open(f_out, "w") as out:
        for contig in output:
            result = ""
            for cat, genes in output[contig].items():
                result += f'{cat}:{genes} '
            result = result.strip()

            out.write(f"{contig}\t{result}\n")

if __name__ == "__main__":

    # Load files
    karga = snakemake.input["karga"]
    alignment = snakemake.input["alignment"]
    f_out = snakemake.output[0]

    # Parse files
    output = parse(karga, alignment, f_out)