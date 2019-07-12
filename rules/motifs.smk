#!/usr/bin/env python

localrules:
    get_meme_sequences

#0. extend peak summit annotation to upstream and downstream distances
#1. if multiple annotations overlap on same strand, keep the one that is the most significant (avoid multiple-counting poorly called peaks erroneously split into multiple peaks)
rule get_meme_sequences:
    input:
        peaks = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}-summits.bed",
        fasta = os.path.abspath(config["genome"]["fasta"])
    output:
        "motifs/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}.fa",
    params:
        search_dist = config["motifs"]["meme-chip"]["search-dist"],
    log:
        "logs/get_meme_sequences/get_meme_sequences_{condition}-v-{control}-{annotation}-{norm}-{direction}-{factor}.log"
    shell: """
        (bedtools slop -b {params.search_dist} -i {input.peaks} -g <(faidx {input.fasta} -i chromsizes) | \
         sort -k1,1 -k2,2n | \
        bedtools cluster -d 0 -i stdin | \
        bedtools groupby -g 7 -c 5 -o max -full -i stdin | \
        sort -k4,4V | \
        bedtools getfasta -name+ -fi {input.fasta} -bed stdin | \
        awk 'BEGIN{{FS=":|-"}} {{if ($1 ~ />/) {{print $1"::"$3":"$4+1"-"$5+1}} else {{print $0}}}}' \
        > {output}) &> {log}
        """

rule meme_chip:
    input:
        seq = "motifs/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}.fa",
        genome_fasta = os.path.abspath(config["genome"]["fasta"]),
        dbs = config["motifs"]["dna_motif_database"] if config["motifs"]["dna_motif_database"] else [],
    output:
        "motifs/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipseq-{norm}-{annotation}-diffbind-results-{direction}-meme_chip/summary.tsv",
    params:
        db_command = "-db" if config["motifs"]["dna_motif_database"] else [],
        meme_mode = config["motifs"]["meme-chip"]["meme-mode"],
        meme_nmotifs = config["motifs"]["meme-chip"]["meme-nmotifs"],
    log:
        "logs/meme_chip/meme_chip_{condition}-v-{control}-{annotation}-{norm}-{direction}-{factor}.log"
    shell: """
        (meme-chip -oc motifs/{wildcards.annotation}/{wildcards.condition}-v-{wildcards.control}/{wildcards.norm}/{wildcards.condition}-v-{wildcards.control}_{wildcards.factor}-chipseq-{wildcards.norm}-{wildcards.annotation}-diffbind-results-{wildcards.direction}-meme_chip {params.db_command} {input.dbs} -bfile <(fasta-get-markov {input.genome_fasta} -m 1) -order 1 -meme-mod {params.meme_mode} -meme-nmotifs {params.meme_nmotifs} -meme-p 1 -meme-norand -centrimo-local {input.seq}) &> {log}
        """

