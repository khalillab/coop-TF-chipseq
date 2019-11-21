#!/usr/bin/env python

localrules:
    normalize_genome_coverage,
    subtract_inputs,
    bedgraph_to_bigwig

# bam must be sorted by name for bedpe.
# Not doing this in the bowtie step is an artifact of samtools indexing required position-sorted bam for splitting species in workflows with spike-ins.
rule get_fragments:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-{{species}}.bam" if SISAMPLES else "alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam"
    output:
        f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-{{species}}-fragments.bedpe"
    threads:
        config["threads"]
    log:
        "logs/get_fragments/get_fragments_{sample}-{species}.log"
    shell: """
        rm -f .{wildcards.sample}_{wildcards.species}*.bam
        (samtools sort -n -T .{wildcards.sample}_{wildcards.species} -@ {threads} {input.bam} | \
         bedtools bamtobed -bedpe -i stdin > {output}) &> {log}
        """

rule midpoint_coverage:
    input:
        bedpe = lambda wc: f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-experimental-fragments.bedpe" if wc.counttype=="counts" else f"alignment/fragments/{{sample}}_{FACTOR}-chipseq-spikein-fragments.bedpe",
        fasta = lambda wc: os.path.abspath(config["genome"]["fasta"]) if wc.counttype=="counts" else config["spike_in"]["fasta"]
    output:
        f"coverage/{{counttype,counts|sicounts}}/{{sample}}_{FACTOR}-chipseq-{{counttype}}-midpoints.bedgraph"
    log:
        "logs/midpoint_coverage/midpoint_coverage_{sample}-{counttype}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} {{width=$6-$2}} {{(width % 2 != 0)? (mid=(width+1)/2+$2) : ((rand()<0.5)? (mid=width/2+$2) : (mid=width/2+$2+1))}} {{print $1, mid, mid+1, $7}}' {input.bedpe} | \
         sort -k1,1 -k2,2n | \
         bedtools genomecov -i stdin -g <(faidx {input.fasta} -i chromsizes) -bga | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule protection_coverage:
    input:
        bam = lambda wc: f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-experimental.bam" if wc.counttype=="counts" and SISAMPLES else f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam" if wc.counttype=="counts" else f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers-spikein.bam"
    output:
        f"coverage/{{counttype,counts|sicounts}}/{{sample}}_{FACTOR}-chipseq-{{counttype}}-protection.bedgraph"
    log:
        "logs/protection_coverage/protection_coverage_{sample}-{counttype}.log"
    shell: """
        (bedtools genomecov -ibam {input.bam} -bga -pc | \
         LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_{factor}-chipseq-counts-{readtype}.bedgraph",
        bam_experimental = "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam",
        bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam" if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_experimental = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-experimental.bam".format(sample=CHIPS[wc.sample]["control"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else [],
        input_bam_spikein = lambda wc: "alignment/{sample}_{factor}-chipseq-uniquemappers-spikein.bam".format(sample=CHIPS[wc.sample]["control"], factor=wc.factor) if wc.norm=="spikenorm" and wc.sample in CHIPS else []
    output:
        normalized = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{readtype}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        readtype="plus|minus|protection|midpoints"
    log:
        "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{readtype}-{factor}.log"
    run:
        if wildcards.norm=="libsizenorm" or wildcards.sample in INPUTS:
            shell("""
                  (awk -v norm_factor=$(samtools view -c {input.bam_experimental} | \
                                        paste -d "" - <(echo "/1000000") | bc -l) \
                   'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
                  """)
        else:
            shell("""
                  (awk -v norm_factor=$(paste -d "" \
                          <(samtools view -c {input.bam_spikein}) <(echo "*") \
                          <(samtools view -c {input.input_bam_experimental}) <(echo "/") \
                          <(samtools view -c {input.input_bam_spikein}) <(echo "/1000000") | bc -l) \
                          'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
                  """)

rule subtract_inputs:
    input:
        ip_sample = "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{readtype}.bedgraph",
        input_sample = lambda wc: f"coverage/{wc.norm}/{{sample}}_{FACTOR}-chipseq-{wc.norm}-{wc.readtype}.bedgraph".format(sample=CHIPS[wc.sample]["control"]),
    output:
        "coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{readtype}-input-subtracted.bedgraph",
    shell: """
        bedtools unionbedg -i {input.ip_sample} {input.input_sample} | \
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4-$5}}' > {output}
        """

rule bedgraph_to_bigwig:
    input:
        bedgraph = f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bedgraph",
        fasta = lambda wc: os.path.abspath(config["spike_in"]["fasta"]) if wc.norm=="sicounts" else os.path.abspath(config["genome"]["fasta"])
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw"
    log:
        "logs/bedgraph_to_bigwig/bedgraph_to_bigwig_{sample}-{readtype}-{norm}.log"
    shell: """
        (bedGraphToBigWig {input.bedgraph} <(faidx {input.fasta} -i chromsizes) {output}) &> {log}
        """

rule smoothed_midpoint_coverage:
    input:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw"
    output:
        f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}_smoothed.bw"
    params:
        bandwidth = config["smooth_bandwidth"]
    conda:
        "../envs/smooth_coverage.yaml"
    log:
        "logs/smoothed_midpoint_coverage/smoothed_midpoint_coverage_{sample}-{norm}-{readtype}.log"
    shell: """
        (python scripts/smooth_midpoint_coverage.py -b {params.bandwidth} -i {input} -o {output}) &> {log}
        """

