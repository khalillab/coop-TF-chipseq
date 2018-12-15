#!/usr/bin/env python

#bam must be sorted by name for bedpe. We don't do this in the bowtie step since samtools indexing required position-sorted bam.
rule get_fragments:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-fragments.bedpe"
    log:
        "logs/get_fragments/get_fragments_{sample}-{species}.log"
    threads:
        config["threads"]
    shell: """
        (samtools sort -n -T .{wildcards.sample}_{wildcards.species} -@ {threads} {input.bam} | bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule midpoint_coverage:
    input:
        bedpe = f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-fragments.bedpe",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        f"coverage/counts/{{sample}}_{FACTOR}-chipseq-midpoint-counts.bedgraph"
    log:
        "logs/midpoint_coverage/midpoint_coverage_{sample}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | sort -k1,1 -k2,2n | bedtools genomecov -i stdin -g <(faidx {input.fasta} -i chromsizes) -bga | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule whole_fragment_coverage:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        f"coverage/counts/{{sample}}_{FACTOR}-chipseq-wholefrag-counts.bedgraph"
    log:
        "logs/whole_fragment_coverage/whole_fragment_coverage_{sample}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -bga -pc | sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = f"coverage/counts/{{sample}}_{FACTOR}-chipseq-{readtype}-counts.bedgraph",
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam",
    output:
        normalized = f"coverage/libsizenorm/{{sample}}_{FACTOR}-{{readtype}}-libsizenorm.bedgraph"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage_{sample}-{readtype}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam}) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{readtype}}-{{norm}}.bedgraph",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-{{readtype}}-{{norm}}.bw"
    log:
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig_{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-midpoint-{{norm}}.bw"
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-midpoint_smoothed-{{norm}}.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    conda:
        "../envs/smooth_coverage.yaml"
    log:
        "logs/smoothed_midpoint_coverage/smoothed_midpoint_coverage_{sample}-{norm}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """
