# files:
#     main_dir
#         |--raw.fastq
#               |--ATAC_seq_rep1_R1.fq.gz
#               |--ATAC_seq_rep1_R2.fq.gz
#               |--ATAC_seq_rep2_R1.fq.gz
#               |--ATAC_seq_rep2_R2.fq.gz
 
# draw topology graph:
#   snakemake -s file.py --dag | dot -Tpdf > test_dag.pdf
# run snakemake:
#   snakemake -s file.py -p -j 2 &
## -j 2: run 2 missions simultaneously

REP_INDEX = {"rep1", "rep2"}
INDEX_BT2 = "directory of bowtie2 genome index"
PICARD = "directory of picard.jar"
ADAPTER_1 = "AGCT..."
ADAPTER_2 = "ACGT..."
NUM_CORES = "4"
NUM_MEMORY = "5"

rule all:
    input:
        expand{"bam/ATAC_seq_{rep}_bt2_hg38_sort.bam", rep = REP_INDEX}
        expand{"bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.bam", rep = REP_INDEX}
        expand{"macs_results/ATAC_seq_{rep}_peaks.narrowPeak", rep = REP_INDEX}

rule cutadapt:
    input: 
        "raw.fastq/ATAC_seq_{rep}_R1.fq.gz",
        "raw.fastq/ATAC_seq_{rep}_R2.fq.gz",
    output: 
        "fix.fasta/ATAC_seq_{rep}_R1.fq.gz",
        "fix.fasta/ATAC_seq_{rep}_R2.fq.gz"
    log:
        "fix.fastq/ATAC_seq_{rep}_cutadapt.log"
    shell:
        "cutadapt -j {NUM_CORES} --times 1 -e 0.1 -O 3 --quality-cutoff 25 \
        -m 50 -a {ADAPTER_1} -A {ADAPTER_2} -o {output[0]} -p {output[1]} \
        {input[0]} {input[1]} > {log} 2>&1"

rule bt2_mapping:
    input:
        "fix.fasta/ATAC_seq_{rep}_R1.fq.gz",
        "fix.fasta/ATAC_seq_{rep}_R2.fq.gz"
    output:
        "bam/ATAC_seq_{rep}_bt2_hg38.sam"
    log:
        "bam/ATAC_seq_{rep}_bt2_hg38.log"
    shell:
        "bowtie2 -x {INDEX_BT2} -p {NUM_CORES} -1 {input[0]} -2 {input[1]} \
        -S {output} > {log} 2>&1"

rule bam_file_sort:
    input:
        "bam/ATAC_seq_{rep}_bt2_hg38.sam"
    output:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort.bam"
    log:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort.log"
    shell:
        "samtools sort -O BAM -o {output} -T {output}.temp \
        -@ {NUM_CORES} -m {NUM_MEMORY}G {input} 2>{log}"

rule remove_duplication:
    input:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort.bam"
    output:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.bam",
        "bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.matrix"
    log:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.log"
    shell:
        "java -Xms{NUM_MEMORY}g -Xmx{NUM_MEMORY}g -XX:ParallelGCThreads={NUM_CORES} \
        -jar {PICARD} MarkDuplicates \
        I={input} O={output[0]} M={output[1]} \
        ASO=coordinate REMOVE_DUPLICATEs=true 2>{log}"

rule call_peak:
    input:
        "bam/ATAC_seq_{rep}_bt2_hg38_sort_rmdup.bam"
    output:
        "macs_results/ATAC_seq_{rep}_peaks.narrowPeak"
    params:
        "ATAC_seq_{rep}",
        "macs_results"
    log:
        "macs_results/ATAC_seq_{rep}_peaks.log"
    shell:
        "macs2 callpeak -t {input} -F BAM -g hs --outdir {params[1]} \
        -n {params[0]} -m 2 100 > {log} 2>&1"