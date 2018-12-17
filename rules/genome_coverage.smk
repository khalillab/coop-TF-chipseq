#!/usr/bin/env python

# bam must be sorted by name for bedpe.
# Not doing this in the bowtie step is an artifact of samtools indexing required position-sorted bam for splitting species in workflows with spike-ins.
rule get_fragments:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-fragments.bedpe"
    log:
        "logs/get_fragments/get_fragments_{sample}.log"
    threads:
        config["threads"]
    shell: """
        (samtools sort -n -T .{wildcards.sample} -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule midpoint_coverage:
    input:
        bedpe = f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-fragments.bedpe",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        f"coverage/counts/{{sample}}_{FACTOR}-chipseq-counts-midpoints.bedgraph"
    log:
        "logs/midpoint_coverage/midpoint_coverage_{sample}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g <(faidx {input.fasta} -i chromsizes) -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule protection_coverage:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        f"coverage/counts/{{sample}}_{FACTOR}-chipseq-counts-protection.bedgraph"
    log:
        "logs/protection_coverage/protection_coverage_{sample}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -bga -pc | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = f"coverage/counts/{{sample}}_{FACTOR}-chipseq-counts-{{readtype}}.bedgraph",
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        normalized = f"coverage/libsizenorm/{{sample}}_{FACTOR}-chipseq-libsizenorm-{{readtype}}.bedgraph"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage_{sample}-{readtype}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam}) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor*1e6; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bedgraph",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw"
    log:
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig_{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-midpoints.bw"
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-midpoints_smoothed.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    conda:
        "../envs/smooth_coverage.yaml"
    log:
        "logs/smoothed_midpoint_coverage/smoothed_midpoint_coverage_{sample}-{norm}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """
