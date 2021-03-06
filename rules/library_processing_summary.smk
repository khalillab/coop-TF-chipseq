#!/usr/bin/env python

localrules:
    aggregate_read_numbers,
    plot_read_processing,
    plot_fragment_lengths,
    build_spikein_counts_table,
    plot_spikein_pct,

rule get_fragment_lengths:
    input:
        expand("alignment/{sample}_{factor}-chipseq-noduplicates.bam", sample=SAMPLES, factor=FACTOR),
    output:
        f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.tsv"
    params:
        header = "\t".join(["fragsize"] + list(SAMPLES.keys()))
    threads: config["threads"]
    run:
        bam = input[0]
        shell("""samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $2, $1}}' > {output}""")
        for bam in input[1:]:
            shell("""join -1 1 -2 2 -t $'\t' -e 0 -a 1 -a 2 --nocheck-order {output} <(samtools view {bam} | cut -f9 | sed 's/-//g' | sort -k1,1n -S 80% --parallel {threads} | uniq -c | awk 'BEGIN{{OFS="\t"}}{{print $1, $2}}') > qual_ctrl/fragment_length_distributions/.frag_length.temp; mv qual_ctrl/fragment_length_distributions/.frag_length.temp {output}""")
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
        markdup = expand("logs/remove_duplicates/remove_duplicates_duplicate_stats_{sample}.log", sample=SAMPLES)
    output:
        f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing_summary.tsv"
    log:
        "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map\tno_dups" > {output}) &> {log}""")
        for sample, adapter, align, markdup in zip(SAMPLES.keys(), input.adapter, input.align, input.markdup):
            shell("""
                  (grep -e "Total read pairs processed:" -e "Pairs written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}
                   grep -e "1 time" {align} | awk 'BEGIN{{sum=0; ORS="\t"}} {{sum+=$1}} END{{print sum}}' >> {output}
                   grep -e "READ:" -e "WRITTEN:" {markdup} | cut -d ' ' -f2 | awk 'BEGIN{{ORS="\t"}} {{print $1/2}} END{{ORS="\\n"; print ""}}' >> {output}) &>> {log}
                   """)
                   # grep -e "exactly 1 time" {align} | awk 'BEGIN{{sum=0; ORS="\t"}} {{sum+=$1}} END{{print sum}}' >> {output}
                   # grep -e "concordantly exactly 1 time" {align} | awk '{{print $1}}' >> {output}) &> {log}

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

rule build_spikein_counts_table:
    input:
        ip_bam_experimental = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-experimental.bam", sample=get_samples(spikein=True, paired=True)),
        ip_bam_spikein = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-spikein.bam", sample=get_samples(spikein=True, paired=True)),
        input_bam_experimental = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-experimental.bam", sample=[v["control"] for k,v in get_samples(spikein=True, paired=True).items()]),
        input_bam_spikein = expand(f"alignment/{{sample}}_{FACTOR}-chipseq-noduplicates-spikein.bam", sample=[v["control"] for k,v in get_samples(spikein=True, paired=True).items()])
    output:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    params:
        groups = [v["group"] for k,v in get_samples(spikein=True, paired=True).items()]
    log:
        "logs/build_spikein_counts_table.log"
    run:
        shell("""(echo -e "sample\tgroup\ttotal_counts_input\texperimental_counts_input\tspikein_counts_input\ttotal_counts_IP\texperimental_counts_IP\tspikein_counts_IP" > {output}) &> {log} """)
        for sample, group, input_exp, input_si ,ip_exp, ip_si in zip(get_samples(spikein=True, paired=True).keys(), params.groups,
                                                                     input.input_bam_experimental, input.input_bam_spikein,
                                                                     input.ip_bam_experimental, input.ip_bam_spikein):
            shell("""(paste <(echo -e "{sample}\t{group}\t") \
                        <(samtools view -c {input_exp}) \
                        <(samtools view -c {input_si}) \
                        <(echo "") \
                        <(samtools view -c {ip_exp}) \
                        <(samtools view -c {ip_si}) | \
                        awk 'BEGIN{{FS=OFS="\t"}} {{$3=$4+$5; $6=$7+$8; print $0}}'>> {output}) &>> {log} """)

rule plot_spikein_pct:
    input:
        f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-counts.tsv"
    output:
        plot = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-plots-{{status}}.svg",
        stats = f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-stats-{{status}}.tsv"
    params:
        samplelist = lambda wc: get_samples(passing=(True if wc.status=="passing" else False), spikein=True, paired=True).keys(),
        conditions = conditiongroups_si if comparisons_si else [],
        controls = controlgroups_si if comparisons_si else []
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/spikein_abundance_chipseq.R"

