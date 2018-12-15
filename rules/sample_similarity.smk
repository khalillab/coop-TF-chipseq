#!/usr/bin/env python

rule map_to_windows:
    input:
        bg = f"coverage/libsizenorm/{{sample}}_{FACTOR}-chipseq-libsizenorm-midpoints.bedgraph",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        temp(f"qual_ctrl/scatter_plots/{FACTOR}_chipseq_{{sample}}-libsizenorm-midpoints-window-{{windowsize}}.bedgraph")
    log: "logs/map_to_windows/map_to_windows_{sample}-{windowsize}.log"
    shell: """
        (bedtools makewindows -g <(faidx {input.fasta} -i chromsizes) -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule join_window_counts:
    input:
        expand("qual_ctrl/scatter_plots/{factor}_chipseq_{sample}-libsizenorm-midpoints-window-{{windowsize}}.bedgraph", sample=SAMPLES, factor=FACTOR)
    output:
        f"qual_ctrl/scatter_plots/{FACTOR}_chipseq_union-bedgraph-libsizenorm-midpoint-window-{{windowsize}}-allsamples.tsv.gz"
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{windowsize}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}) &> {log}
        """

rule plot_scatter_plots:
    input:
        "qual_ctrl/scatter_plots/{factor}_chipseq_union-bedgraph-libsizenorm-midpoint-window-{windowsize}-allsamples.tsv.gz"
    output:
        "qual_ctrl/scatter_plots/{condition}-v-{control}/{status}/{condition}-v-{control}_{factor}_chipseq-libsizenorm-scatterplots-{status}-window-{windowsize}.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = lambda wc: get_samples(wc.status, [wc.condition, wc.control])
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/plot_scatter_plots.R"

