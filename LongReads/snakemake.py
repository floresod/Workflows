import os
import glob 

SAMPLE,=glob_wildcards("RawReads/{sample}.fastq.gz")

rule all: 
    input: 
        expand("RawQC/{sample}_fastqc.{extension}", sample=SAMPLE, extension=["html", "zip"])

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
