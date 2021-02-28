#!/usr/bin/env python

localrules:
    build_combined_genome,
    bowtie2_build,

basename = "{exp_name}_{exp_fasta}_{si_name}_{si_fasta}".format(exp_name = config["genome"]["name"],
                                                                exp_fasta = os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0],
                                                                si_name = config["spike_in"]["name"],
                                                                si_fasta = os.path.splitext(os.path.basename(config["spike_in"]["fasta"]))[0]) if SISAMPLES else os.path.splitext(os.path.basename(config["genome"]["fasta"]))[0]

rule build_combined_genome:
    input:
        experimental = os.path.abspath(config["genome"]["fasta"]),
        spikein = config["spike_in"]["fasta"] if SISAMPLES else []
    output:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(config["genome"]["fasta"]))[0], bn=basename),
    params:
        exp_name = config["genome"]["name"],
        si_name = config["spike_in"]["name"] if SISAMPLES else []
    log: "logs/build_combined_genome.log"
    shell: """
        (sed 's/>/>{params.exp_name}_/g' {input.experimental} | \
        cat - <(sed 's/>/>{params.si_name}_/g' {input.spikein}) > {output}) &> {log}
        """

rule bowtie2_build:
    input:
        "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(config["genome"]["fasta"]))[0], bn=basename) if SISAMPLES else os.path.abspath(config["genome"]["fasta"]),
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    conda:
        "../envs/bowtie2.yaml"
    log:
        "logs/bowtie2_build-{basename}.log"
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
        unaligned_r1 = f"fastq/unaligned/{{sample}}_{FACTOR}-chipseq-unaligned.r1.fastq.gz",
        unaligned_r2 = f"fastq/unaligned/{{sample}}_{FACTOR}-chipseq-unaligned.r2.fastq.gz",
        aligned_r1 = f"fastq/aligned/{{sample}}_{FACTOR}-chipseq-aligned.r1.fastq.gz",
        aligned_r2 = f"fastq/aligned/{{sample}}_{FACTOR}-chipseq-aligned.r2.fastq.gz",
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
        (bowtie2 --minins {params.min_fraglength} --maxins {params.max_fraglength} --fr --no-mixed --no-discordant --al-conc-gz fastq/aligned/{wildcards.sample}_{FACTOR}-chipseq-aligned.fastq.gz --un-conc-gz fastq/unaligned/{wildcards.sample}_{FACTOR}-chipseq-unaligned.fastq.gz -p {threads} -x {params.idx_path}/{basename} -1 {input.r1} -2 {input.r2}  | \
         samtools view -buh -q {params.minmapq} - | \
         samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        mv fastq/aligned/{wildcards.sample}_{FACTOR}-chipseq-aligned.fastq.1.gz {output.aligned_r1}
        mv fastq/aligned/{wildcards.sample}_{FACTOR}-chipseq-aligned.fastq.2.gz {output.aligned_r2}
        mv fastq/unaligned/{wildcards.sample}_{FACTOR}-chipseq-unaligned.fastq.1.gz {output.unaligned_r1}
        mv fastq/unaligned/{wildcards.sample}_{FACTOR}-chipseq-unaligned.fastq.2.gz {output.unaligned_r2}
        """

#indexing is required for separating species by samtools view
rule remove_duplicates:
    input:
        f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates.bam",
        index = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates.bam.csi",
        markdup_log = "logs/remove_duplicates/remove_duplicates_duplicate_stats_{sample}.log"
    log:
        "logs/remove_duplicates/remove_duplicates-{sample}.log"
    threads:
        config["threads"]
    shell: """
        (samtools collate -O -u --threads {threads} {input} | \
                samtools fixmate -m -u --threads {threads} - - | \
                samtools sort -u -@ {threads} | \
                samtools markdup -r -f {output.markdup_log} -d 100 -m t --threads {threads} --write-index - {output.bam}) &> {log}
        """

rule bam_separate_species:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates.bam",
        index = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates.bam.csi",
        fasta = "{directory}/{bn}.fa".format(directory = os.path.split(os.path.abspath(config["genome"]["fasta"]))[0], bn=basename) if SISAMPLES else [],
    output:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-{{species}}.bam",
        index = f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-{{species}}.bam.csi",
    params:
        filterprefix = lambda wc: config["spike_in"]["name"] if wc.species=="experimental" else config["genome"]["name"],
        prefix = lambda wc: config["genome"]["name"] if wc.species=="experimental" else config["spike_in"]["name"]
    threads:
        config["threads"]
    log:
        "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h -@ {threads} {input.bam} $(faidx {input.fasta} -i chromsizes | \
                                                     grep {params.prefix}_ | \
                                                     awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | \
         grep -v -e 'SN:{params.filterprefix}_' | \
         sed 's/{params.prefix}_//g' | \
         samtools view -bh -@ {threads} -o {output.bam} -) &> {log}
        """

