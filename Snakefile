import os

# Config file with run parameters
configfile: "config.yaml"


############### Overall rules ##################

# Trigger run start with all samples
rule run:
    message: "Starting Nanoplasm"
    input: "Typing_results.tsv",
            "08-mge-cluster/mge-cluster_results.csv",
            expand("09-karga/{sample}_nanofilt_KARGA_mappedReads.csv", sample=config["samples"]),
            expand("09-karga/{sample}_contigs_amr_profile.tsv", sample=config["samples"]),
            ".annotations.txt"

# Trigger the quality check rules
rule quality_check:
    message: "Quality check"
    input: expand("QC/{sample}_nanoq_report.txt", sample=config["samples"]),
           "QCsummary.txt"


############### Rules common to long and hybrid mode ###############


# LONG/HYBRID MODE -- Reads filtering with Nanofilt
rule nanofilt:
    message: "NanoFilt: {wildcards.sample}"
    input: config["path_nano"]+"/{sample}.fastq"
    output: "01-NanoFilt/{sample}_nanofilt.fastq"
    log: "01-NanoFilt/{sample}_nanofilt.log"
    container: "docker://mcfonsecalab/nanofilt"        
    threads: 2
    params:
        quality = config["nanofilt"]["quality"],
        length = config["nanofilt"]["length"]
    shell:
        "NanoFilt -q {params.quality} -l {params.length} {input} --logfile {log}> {output}"

# LONG/HYBRID MODE -- Assembly polishing with Medaka
rule medaka:
    message: "Medaka: {wildcards.sample}"
    input:
        reads = "01-NanoFilt/{sample}_nanofilt.fastq",
        assembly = "03-assembly/{sample}_long/assembly.fasta" if config["mode"] == "long"  else "03-assembly/{sample}_hybrid/assembly.fasta"
    output: "04-Medaka/{sample}/consensus.fasta"
    log: "04-Medaka/{sample}/{sample}_medaka_log.txt"
    container: "docker://ontresearch/medaka"
    threads: 2
    params:
        model = config["medaka"]["model"]
    shell:
        "medaka_consensus -i {input.reads} -d {input.assembly} -t {threads} -m {params.model} -o 04-Medaka/{wildcards.sample} > {log} 2>&1"

############### Rules for long mode ###############


# LONG MODE -- Reads quality with Nanoq
rule nanoq:
    message: "Nanoq: {wildcards.sample}"
    input: "01-NanoFilt/{sample}_nanofilt.fastq"
    output: "QC/{sample}_nanoq_report.txt"
    container: "docker://jimmyliu1326/nanoq"        
    shell:
        "nanoq -i {input} -s -vvv 2> {output}"

# LONG MODE -- QC summary
rule qc_summary:
    message: "Qc summary"
    input: 
        nanoq = expand("QC/{sample}_nanoq_report.txt", sample=config["samples"]),
        flye = expand("03-assembly/{sample}_long/assembly_info.txt", sample=config["samples"])
    output: "QCsummary.txt"
    params:
        plasmid_min_size = config["classify_contigs"]["min_size"],
        plasmid_max_size = config["classify_contigs"]["max_size"]
    script:
        "bin/QCsummary.py"

# LONG MODE -- Assembly with Flye
rule flye:
    message: "Flye: {wildcards.sample}"
    input: "01-NanoFilt/{sample}_nanofilt.fastq"
    output: 
        assembly = "03-assembly/{sample}_long/assembly.fasta",
        info = "03-assembly/{sample}_long/assembly_info.txt"
    log: "03-assembly/{sample}_long/{sample}_flye_log.txt"
    container: "docker://staphb/flye"
    threads: 12
    params:
        others = config["flye"]["others"]
    shell:
        "flye --nano-hq {input} -t {threads} -o 03-assembly/{wildcards.sample}_long {params.others} > {log} 2> /dev/null"

# LONG MODE -- Homopolish
rule homopolish:
    message: "Homopolish: {wildcards.sample}"
    input: "04-Medaka/{sample}/consensus.fasta"
    output: "05-Homopolish/{sample}/consensus_homopolished.fasta"
    log: "05-Homopolish/{sample}/{sample}_homopolish_log.txt"
    container: "docker://staphb/homopolish"
    threads: 12
    params:
        model = config["homopolish"]["model"],
        db = config["homopolish"]["db"]
    shell:
        "homopolish polish -a {input} -s {workflow.basedir}/{params.db} -m R10.3.pkl -o 05-Homopolish/{wildcards.sample} > {log} 2>&1"


############### Rules for hybrid mode ###############


# HYBRID MODE -- Reads filtering with filFastp
rule fastp:
    message: "Fastp: {wildcards.sample}"
    input: 
        R1 = config["path_ill"]+"/{sample}_R1.fastq",
        R2 = config["path_ill"]+"/{sample}_R2.fastq"
    output: 
        R1 = "02-Fastp/{sample}_R1_fastp.fastq",
        R2 = "02-Fastp/{sample}_R2_fastp.fastq",
        html = "02-Fastp/{sample}_fastp.html",
        json = "02-Fastp/{sample}_fastp.json"
    log: "02-Fastp/{sample}_fastp_log.txt"
    container: "docker://staphb/fastp"
    threads: 1
    params:
        qual = config["fastp"]["qualified_quality_phred"],
        unqual = config["fastp"]["unqualified_percent_limit"],
        ave_qual = config["fastp"]["average_qual"],
        length = config["fastp"]["length_limit"]
    shell:
        "fastp -i {input.R1} -I {input.R2} -o {output.R1} -O {output.R2} -q {params.qual} "
        "-u {params.unqual} -e {params.ave_qual} -l {params.length} -h {output.html} -j {output.json} 2> {log}"

# HYBRID MODE -- Reads quality check with FastQC
rule fastqc:
    message: "FastQC: {wildcards.sample}"
    input: 
        config["path_ill"]+"/{sample}_R1.fastq",
        config["path_ill"]+"/{sample}_R2.fastq"
    output: "QC/{sample}_R1_fastqc.zip"
    container: "docker://staphb/fastqc"
    threads: 1
    shell:
        "mkdir QC/{wildcards.sample};"
        "fastqc {input} -o QC > /dev/null 2>&1"

# HYBRID MODE -- Assembly with Unicycler
rule unicycler:
    message: "Unicycler: {wildcards.sample}"
    input: 
        short_1 = "02-Fastp/{sample}_R1_fastp.fastq",
        short_2 = "02-Fastp/{sample}_R2_fastp.fastq",
        long = "01-NanoFilt/{sample}_nanofilt.fastq"
    output: "03-assembly/{sample}_hybrid/assembly.fasta"
    container: "docker://staphb/unicycler"
    threads: 12
    params:
        mode = config["unicycler"]["mode"]
    shell:
        "unicycler -1 {input.short_1} -2 {input.short_2} -l {input.long} -o 03-assembly/{wildcards.sample}_hybrid -t {threads} --mode {params.mode} "
        "--verbosity 0 --keep 0"

# HYBRID MODE -- Short reads alignments for Polypolish
rule polypolish_prep:
    message: "Reads alignment for polypolish: {wildcards.sample}"
    input: 
        draft = "04-Medaka/{sample}/consensus.fasta",
        short_1 = "02-Fastp/{sample}_R1_fastp.fastq",
        short_2 = "02-Fastp/{sample}_R2_fastp.fastq"
    output: 
        alignments_1 = "05-polypolish/{sample}/alignments_1.sam",
        alignments_2 = "05-polypolish/{sample}/alignments_2.sam"
    log: "05-polypolish/{sample}/{sample}_BWA.log"
    threads: 12
    shell:
        "bwa index {input.draft} > /dev/null 2>&1;"
        "bwa mem -t {threads} -a {input.draft} {input.short_1} > {output.alignments_1} 2> {log};"
        "bwa mem -t {threads} -a {input.draft} {input.short_1} > {output.alignments_2} 2> {log}"

# HYBRID MODE -- Assembly polyshing with polypolish
rule polypolish:
    message: "Polypolish: {wildcards.sample}"
    input:
        draft = "04-Medaka/{sample}/consensus.fasta",
        alignments_1 = "05-polypolish/{sample}/alignments_1.sam",
        alignments_2 = "05-polypolish/{sample}/alignments_2.sam"
    output: 
        alignments_filt_1 = "05-polypolish/{sample}/filtered_1.sam",
        alignments_filt_2 = "05-polypolish/{sample}/filtered_2.sam",
        assembly = "05-polypolish/{sample}/assembly.fasta"
    log: "05-polypolish/{sample}/{sample}_polypolish.log"
    container: "docker://staphb/polypolish"
    threads: 1
    shell:
        "polypolish_insert_filter.py --in1 {input.alignments_1} --in2 {input.alignments_2} --out1 {output.alignments_filt_1} --out2 {output.alignments_filt_2} > /dev/null 2>&1;"
        "polypolish {input.draft} {output.alignments_filt_1} {output.alignments_filt_2} > {output.assembly} 2> {log}"


############### Rules for plasmids analysis ###############

# ResFinder - Antimicrobial resistance genes
rule resfinder:
    message: "Resfinder: {wildcards.sample}"
    input: "05-Homopolish/{sample}/consensus_homopolished.fasta" if config["mode"] == "long" else "05-polypolish/{sample}/assembly.fasta"
    output: "06-resfinder/{sample}/ResFinder_results_tab.txt"
    container: "docker://genomicepidemiology/resfinder"
    params:
        id_threshold = config["resfinder"]["id_threshold"],
        coverage_threshold = config["resfinder"]["coverage_threshold"],
        db = config["resfinder"]["db"]
    shell:
        "python -m resfinder -ifa {input} --nanopore -acq -o 06-resfinder/{wildcards.sample} -l {params.coverage_threshold} -t {params.id_threshold} > /dev/null 2>&1"

# KARGA - Antimicrobial resistance genes
rule karga:
    message: "KARGA: {wildcards.sample}"
    input: "01-NanoFilt/{sample}_nanofilt.fastq"
    output: "09-karga/{sample}_nanofilt_KARGA_mappedReads.csv"
    threads: 24
    params:
        db = config["karga"]["db"]
    shell:
        """java -cp {workflow.basedir}/KARGA/openjdk-8/KARGA KARGA {input} d:{workflow.basedir}/{params.db} -Xmx16GB;
        mv 01-NanoFilt/{wildcards.sample}_nanofilt_KARGA_mappedReads.csv 09-karga;
        mv 01-NanoFilt/{wildcards.sample}_nanofilt_KARGA_mappedGenes.csv 09-karga"""

# KARGA - Align fastq on assembly with minimap2 to get reads-plasmids relationship
rule minimap2:
    message: "Minimap2: {wildcards.sample}"
    input: 
        fastq = "01-NanoFilt/{sample}_nanofilt.fastq",
        assembly = "03-assembly/{sample}_long/assembly.fasta"
    output: "09-karga/{sample}_contigs_on_assembly.paf"
    threads: 24
    container: "docker://nanozoo/minimap2"
    shell:
        "minimap2 {input.assembly} {input.fastq} -o {output} -t {threads}"

rule reads_to_contigs:
    message: "Link reads and contigs: {wildcards.sample}"
    input: 
        alignment = "09-karga/{sample}_contigs_on_assembly.paf",
        karga = "09-karga/{sample}_nanofilt_KARGA_mappedReads.csv"
    output: "09-karga/{sample}_contigs_amr_profile.tsv"
    script:
        "bin/reads_to_contigs.py"

# Mob-suite - plasmid typing
rule mobsuite:
    message: "Mob-Suite"
    input: "05-Homopolish/{sample}/consensus_homopolished.fasta" if config["mode"] == "long" else "05-polypolish/{sample}/assembly.fasta"
    output: "07-mobsuite/{sample}_typing.txt"
    container: "docker://kumalpha/mob-suite"
    threads: 24
    log: "07-mobsuite/{sample}_typing.log"
    shell:
        "mob_typer -i {input} -o {output} -n {threads} --multi > {log} 2>&1"

# Contigs classification
rule classification:
    message: "Contigs classification"
    input: expand("05-Homopolish/{sample}/consensus_homopolished.fasta", sample=config["samples"]) if config["mode"] == "long" 
           else expand("05-polypolish/{sample}/assembly.fasta", sample=config["samples"])
    output: 
        class_file = "Sequences/contigs_classification.tsv",
        plasmids_list = ".plasmids_list.txt"
    script:
        "bin/classify_contigs.py"

# Mge-cluster - clustering plasmids sequences based on unitigs
rule mgecluster:
    message: "Mge-cluster"
    input: ".plasmids_list.txt"
    output: "08-mge-cluster/mge-cluster_results.csv"
    params:
        perplexity = config["mge-cluster"]["perplexity"],
        min_cluster = config["mge-cluster"]["min_cluster"],
        kmer = config["mge-cluster"]["kmer"]
    threads: 24
    shell:
        "mge_cluster --create --input {input} --outdir 08-mge-cluster --perplexity {params.perplexity} --min_cluster {params.min_cluster} --kmer {params.kmer} --threads {threads}"

# Gathering results from Mobsuite and Resfinder
rule gather_results:
    message: "Gathering results"
    input: mobsuite = expand("07-mobsuite/{sample}_typing.txt", sample=config["samples"]),
           resfinder = expand("06-resfinder/{sample}/ResFinder_results_tab.txt", sample=config["samples"]),
           karga = expand("09-karga/{sample}_contigs_amr_profile.tsv", sample=config["samples"]),
           classification = "Sequences/contigs_classification.tsv"
    output: "Typing_results.tsv"
    script:
        "bin/gather_typing_results.py"

# PROKKA - annotation of plasmids
rule prokka:
    message: "Prokka"
    input: "Sequences/contigs_classification.tsv"
    output: ".annotations.txt"
    container: "docker://staphb/prokka"
    params:
        options = config["prokka"]["options"]
    threads: 48
    shell:
        """
        mkdir -p 10-prokka;
        for i in $(ls Sequences/plasmids); 
            do prokka --force --outdir 10-prokka/$i -cpus {threads} {params.options} Sequences/plasmids/$i;
        done;
        touch {output}
        """