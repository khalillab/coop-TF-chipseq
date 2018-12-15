#!/usr/bin/env python

localrules:
    bowtie2_build,

basename = os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0]

rule bowtie2_build:
    input:
        os.path.abspath(config["genome"]["fasta"]),
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    conda: "../envs/bowtie2.yaml"
    log: "logs/bowtie2_build-{basename}.log"
    shell: """
        (bowtie2-build {input} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand("{directory}/{bn}.{num}.bt2", directory = config["bowtie2"]["index-path"],
                                             bn = basename,
                                             num = [1,2,3,4]),
        expand("{directory}/{bn}.rev.{num}.bt2", directory = config["bowtie2"]["index-path"],
                                             bn = basename,
                                             num = [1,2]),
        r1 = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.r1.fastq.gz",
        r2 = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.r2.fastq.gz",
    output:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
        unaligned_fastq = f"fastq/unaligned/{{sample}}_{FACTOR}-chipseq-unaligned.fastq.gz",
        aligned_fastq = f"fastq/aligned/{{sample}}_{FACTOR}-chipseq-aligned.fastq.gz",
        log = "logs/align/align_{sample}.log"
    params:
        idx_path = config["bowtie2"]["index-path"],
        minmapq = config["bowtie2"]["minmapq"],
        min_fraglength = config["bowtie2"]["min_fraglength"],
        max_fraglength = config["bowtie2"]["max_fraglength"],
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config["threads"]
    shell: """
        (bowtie2 --minins {params.min_fraglength} --maxins {params.max_fraglength} --fr --no-mixed --no-discordant --al-conc-gz {output.aligned_fastq} --un-conc-gz {output.unaligned_fastq} -p {threads} -x {params.idx_path}/{basename} -1 {input.r1} -2 {input.r2}  | samtools view -buh -q {params.minmapq} - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        """

##indexing is required for separating species by samtools view
#rule index_bam:
#    input:
#        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam"
#    output:
#        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam.bai"
#    log : "logs/index_bam/index_bam-{sample}.log"
#    shell: """
#        (samtools index {input}) &> {log}
#        """

