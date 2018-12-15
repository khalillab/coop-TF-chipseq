#!/usr/bin/env python

localrules:
    aggregate_read_numbers,
    plot_read_processing,
    plot_fragment_lengths

rule get_fragment_lengths:
    input:
        expand("alignment/{sample}_{factor}-chipseq-uniquemappers.bam", sample=SAMPLES, factor=FACTOR),
    output:
        f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.tsv"
    params:
        header = "\t".join(["fragsize"] + list(SAMPLES.keys()))
    threads: config["threads"]
    run:
        bam = input[0]
        shell("""samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' > {output}""")
        for bam in input[1:]:
            shell("""join -1 1 -2 2 -t $'\t' -e 0 -a 1 -a 2 {output} <(samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}') > qual_ctrl/fragment_length_distributions/.frag_length.temp; mv qual_ctrl/fragment_length_distributions/.frag_length.temp {output}""")
        shell("""sed -i "1i {params.header}" {output}""")

rule plot_fragment_lengths:
    input:
        table = f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.tsv"
    output:
        plot = f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.svg"
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/paired_end_fragment_length.R"

rule aggregate_read_numbers:
    input:
        adapter = expand("logs/clean_reads/clean_reads-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align_{sample}.log", sample=SAMPLES),
    output:
        f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing_summary.tsv"
    log:
        "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map\tpaired" > {output}) &> {log}""")
        for sample, adapter, align in zip(SAMPLES.keys(), input.adapter, input.align):
            shell("""
                  (grep -e "Total read pairs processed:" -e "Pairs written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}
                   grep -e "1 time" {align} | awk 'BEGIN{{sum=0; ORS="\t"}} {{sum+=$1}} END{{print sum}}' >> {output}
                   grep -e "exactly 1 time" {align} | awk 'BEGIN{{sum=0; ORS="\t"}} {{sum+=$1}} END{{print sum}}' >> {output}
                   grep -e "concordantly exactly 1 time" {align} | awk '{{print $1}}' >> {output}) &> {log}
                   """)

rule plot_read_processing:
    input:
        f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing_summary.tsv"
    output:
        surv_abs_out = f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing-survival-absolute.svg",
        surv_rel_out = f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing-survival-relative.svg",
        loss_out  = f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing-loss.svg",
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/processing_summary.R"

