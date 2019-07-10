# analysis process of Bismark software
# Prerequisite:
#           samtools, bowtie2
# Download and install:
#           tar -zxvf bismark_version.tar.gz

rule all:
    input:
        ""

rule genome_preparation:
# under install folder
    shell:
        "./bismark_genome_preparation --path_to_bowtie2 {BOWTIE2} --verbose {GENOME}"

rule alignment:
    input:
        "raw.fastq/read1.fq.gz",
        "raw.fastq/read2.fq.gz"
    output:
        "alignment/"
    log:
        "alignment/align.log"
    shell:
        "./bismark --genome {GENOME} -1 {input[0]} -2 {input[1]} -p {CORES} \
        -o {output} 2>{log}"

rule methylation_extract:
    input:
        "alignment"
