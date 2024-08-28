import os
import glob 

SAMPLE,=glob_wildcards("RawReads/{sample}.fastq.gz")

rule all: 
    input: 
        expand("RawQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html", "zip"]), 
#        expand("FilteredReads/{sample}.fastq.gz", sample=SAMPLE),
#        expand("FilteredQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html","zip"]),
        expand("Kraken_report/{sample}.txt", sample=SAMPLE), 
        expand("Bracken_report/{sample}.txt", sample=SAMPLE)


rule fastqc_rawreads: 
    input: 
        rawread="RawReads/{sample}.fastq.gz"
    output: 
        zip="RawQC/{sample}_fastqc.zip",
        html="RawQC/{sample}_fastqc.html"
    log: "logs/fastqc_{sample}.log"
    threads:
        1
    params:
        path="RawQC"
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path} 2>{log}
        """ 

#rule chopper_run: 
#    input:
#        rawread="RawReads/{sample}.fastq.gz"
#    output:
#        FilteredRead="FilteredReads/{sample}.fastq.gz"
#    params:
#        threads=5
#    shell: 
#        """
#        gunzip -c {input.rawread} | \
#        chopper -q 10 -l 100 --tailcrop 10 --threads {params.threads} | \
#        gzip > {output.FilteredRead}
#        """ 

#rule fastqc_filtered: 
#    input:
#        FilteredRead="FilteredReads/{sample}.fastq.gz"
#    output:
#        zip="FilteredQC/{sample}_fastqc.zip",
#        html="FilteredQC/{sample}_fastqc.html"
#    
#    params:
#        path="FilteredQC"
#    threads:
#        5
#        
#    shell: 
#        """
#        fastqc {input.FilteredRead} --threads {threads} -o {params.path} 
 #       """

rule kraken2_run: 
    input: 
        rawread="RawReads/{sample}.fastq.gz"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        threads=5
    output:
        output="Kraken_output/{sample}.txt",
        report="Kraken_report/{sample}.txt" 
    log: "logs/kraken2_{sample}.log"
    shell:
        """
        kraken2 --db {params.database} --gzip-compressed {input.rawread} \
        --output {output.output} --report {output.report} \
        --threads {params.threads} \
        --confidence 0.01 \
        --use-names  2>{log} 
        """
rule bracken_run:
    input:
        report="Kraken_report/{sample}.txt"
    output:
        report="Bracken_report/{sample}.txt"
    params:
        database="../../Databases/k2_standard_08gb_20240605",
        length=100
    log: "logs/bracken_{sample}.log"

    shell:
        """
        bracken -d {params.database} -i {input.report} -o {output.report} -r {params.length}  2>{log}
        """


