#!/usr/bin/env python

localrules:
    get_meme_sequences

#0. extend peak summit annotation to upstream and downstream distances
#1. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting poorly called peaks erroneously split into multiple peaks)
rule get_meme_sequences:
    input:
        peaks = "diff_exp/peaks/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-peaks-diffexp-results-{category}-{direction}-summits.bed",
        fasta = os.path.abspath(build_annotations(config["genome"]["fasta"]))
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa"
    params:
        upstr = config["motifs"]["meme-chip"]["upstream"],
        dnstr = config["motifs"]["meme-chip"]["downstream"]
    log:
        "logs/get_meme_sequences/get_meme_sequences_{condition}-v-{control}-{norm}-{category}-{direction}.log"
    shell: """
        (bedtools slop -l {params.upstr} -r {params.dnstr} -s -i {input.peaks} -g <(faidx {input.fasta} -i chromsizes) | \
        bedtools cluster -s -d 0 -i stdin | \
        bedtools groupby -g 7 -c 5 -o max -full -i stdin | \
        sort -k4,4V | \
        bedtools getfasta -name+ -fi {input.fasta} -bed stdin | \
        awk 'BEGIN{{FS=":|-"}} {{if ($1 ~ />/) {{print $1"::"$3":"$4+1"-"$5+1}} else {{print $0}}}}' \
        > {output}) &> {log}
        """

rule meme_chip:
    input:
        seq = "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}.fa",
        genome_fasta = os.path.abspath(build_annotations(config["genome"]["fasta"])),
        dbs = build_annotations("motifs/" + config["genome"]["name"] + "_all_dna_motifs.meme") if config["motifs"]["dna_motif_databases"] else [],
    output:
        "motifs/{condition}-v-{control}/{norm}/{category}/{condition}-v-{control}_tss-seq-{norm}-diffexp-results-{category}-{direction}-meme_chip/summary.tsv"
    params:
        db_command = "-db" if config["motifs"]["dna_motif_databases"] else [],
        meme_mode = config["motifs"]["meme-chip"]["meme-mode"],
        meme_nmotifs = config["motifs"]["meme-chip"]["meme-nmotifs"],
    log:
        "logs/meme_chip/meme_chip_{condition}-v-{control}-{norm}-{category}-{direction}.log"
    shell: """
        (meme-chip -oc motifs/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}/{wildcards.category}/{wildcards.condition}-v-{wildcards.control}_tss-seq-{wildcards.norm}-diffexp-results-{wildcards.category}-{wildcards.direction}-meme_chip {params.db_command} {input.dbs} -bfile <(fasta-get-markov {input.genome_fasta} -m 1) -order 1 -meme-mod {params.meme_mode} -meme-nmotifs {params.meme_nmotifs} -meme-p 1 -meme-norand -centrimo-local {input.seq}) &> {log}
        """

