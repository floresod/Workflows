import os
import glob 

SAMPLE,=glob_wildcards("RawReads/{sample}.fastq.gz")

rule all: 
    input: 
        expand("RawQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html", "zip"]), 
        expand("FilteredReads/{sample}.fastq.gz", sample=SAMPLE)


rule fastqc_rawreads: 
    input: 
        rawread="RawReads/{sample}.fastq.gz"
    output: 
        zip="RawQC/{sample}_fastqc.zip",
        html="RawQC/{sample}_fastqc.html"
    threads:
        1
    params:
        path="RawQC"
    shell:
        """
        fastqc {input.rawread} --threads {threads} -o {params.path}
        """ 

rule chopper_run: 
    input:
        rawread="RawReads/{sample}.fastq.gz"
    output:
        FilteredRead="FilteredReads/{sample}.fastq.gz"
    params:
        threads=5
    shell: 
        """
        gunzip -c {input.rawread} | \
        chopper -q 10 -l 100 --tailcrop 10 --threads {params.threads} | \
        gzip > {output.FilteredRead}
        """
