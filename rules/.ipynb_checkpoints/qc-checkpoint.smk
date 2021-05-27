# Snakemake rules imported in the main Snakefile to analyze read quality with fastqc and multiqc

rule fastqc:
    ''' Older versions has qc before cutadapt / this change will allow to get rid of N and low quality reads in the right and left '''
    input: 
        config["path"] + "02_quality/{sample}_{i}.fastq.gz"
    output: 
        html= temp(config["path"] + "QC/{sample}_{i}_fastqc.html"),
        zip = temp(config["path"] + "QC/{sample}_{i}_fastqc.zip")
    wrapper: "0.72.0/bio/fastqc"

rule multiqc:
    input: 
        expand([config["path"] + "QC/{sample}_{i}_fastqc.html", 
                config["path"] + "QC/{sample}_{i}_fastqc.zip",
                config["path"] + "QC/{sample}.qc.txt"],sample=sample_id,i=i)
    output: 
        config["path"] + "QC/multiqc_report.html"
    wrapper:  "0.72.0/bio/multiqc"
