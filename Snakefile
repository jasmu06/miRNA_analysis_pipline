SAMPLES = ["tumor1", "tumor2", "tumor3", "control1", "control2", "control3"]







rule fastqc:
    input:
        "tumor_mirna/data/raw_fastq/{sample}.fastq"
    output:
        html = "tumor_mirna/data/fastqc_reports/{sample}_fastqc.html",
        zip = "tumor_mirna/data/fastqc_reports/{sample}_fastqc.zip"
    shell:
        "fastqc {input} --outdir=tumor_mirna/data/fastqc_reports"




rule trim_reads:
    input:
        "tumor_mirna/data/raw_fastq/{sample}.fastq"
    output:
        "tumor_mirna/data/trimmed_fastq/{sample}_trimmed.fastq"
    params:
        adapter = "TGGAATTCTCGGGTGCCAAGG"  # example: universal miRNA adapter, replace if needed
    shell:
        "cutadapt -a {params.adapter} -o {output} {input}"


 
