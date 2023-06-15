# Config file with run parameters
configfile: "config.yaml"


############### Overall rules ##################


# Trigger run start with all samples
if config["mode"] == "long":
    final_input = expand("05-Homopolish/{sample}/consensus_homopolished.fasta", sample=config["samples"])
elif config["mode"] == "hybrid":
    final_input = expand("05-polypolish/{sample}/assembly.fasta", sample=config["samples"])

rule run:
    message: "Starting Nanoplasm"
    input: final_input

# Trigger the quality check rules
rule quality_check:
    message: "Quality check"
    input: expand("QC/{sample}_nanoq_report.txt", sample=config["samples"]),
           expand("QC/{sample}_R1_fastqc.zip", sample=config["samples"]) if config["mode"] == "hybrid" else ""


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
    threads: 12
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

# LONG MODE -- Assembly with Flye
rule flye:
    message: "Flye: {wildcards.sample}"
    input: "01-NanoFilt/{sample}_nanofilt.fastq"
    output: "03-assembly/{sample}_long/assembly.fasta"
    log: "03-assembly/{sample}/{sample}_flye_log.txt"
    container: "docker://staphb/flye"
    threads: 12
    params:
        others = config["flye"]["others"]
    shell:
        "flye --nano-hq {input} -t {threads} -o 03-assembly/{wildcards.sample} {params.others} > {log}"

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
        "homopolish polish -a {input} -s {workflow.basedir}/{params.db} -m R10.3.pkl -o 05-Homopolish/{wildcards.sample} > {log}"


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
        "unicycler -1 {input.short_1} -2 {input.short_2} -l {input.long} -o 03-assembly/{wildcards.sample} -t {threads} --mode {params.mode} "
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
    threads: 12
    shell:
        "bwa index {input.draft} > /dev/null 2>&1;"
        "bwa mem -t {threads} -a {input.draft} {input.short_1} > {output.alignments_1} 2> /dev/null;"
        "bwa mem -t {threads} -a {input.draft} {input.short_1} > {output.alignments_2} 2> /dev/null"

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
    container: "docker://staphb/polypolish"
    threads: 1
    shell:
        "polypolish_insert_filter.py --in1 {input.alignments_1} --in2 {input.alignments_2} --out1 {output.alignments_filt_1} --out2 {output.alignments_filt_2} > /dev/null 2>&1;"
        "polypolish {input.draft} {output.alignments_filt_1} {output.alignments_filt_2} > {output.assembly} 2> /dev/null"