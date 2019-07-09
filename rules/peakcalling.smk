#!/usr/bin/env python

localrules: combine_peaks

rule callpeaks_macs2:
    input:
        chip_bam = lambda wc: expand("alignment/{sample}_{factor}-chipseq-uniquemappers-{species}.bam", sample=[k for k,v in CHIPS_PASSING.items() if v["group"]==wc.group], factor=FACTOR, species=wc.species),
        input_bam = lambda wc: expand("alignment/{sample}_{factor}-chipseq-uniquemappers-{species}.bam", sample=[k for k,v in INPUTS_PASSING.items() if v["group"]==wc.group], factor=FACTOR, species=wc.species),
        fasta = lambda wc: os.path.abspath(config["genome"]["fasta"]) if wc.species=="experimental" else os.path.abspath(config["spike_in"]["fasta"]),
    output:
        tsv = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_peaks.xls",
        peaks = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_peaks.narrowPeak",
        summits = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_summits.bed",
        treat_bg = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_treat_pileup.bdg",
        cntrl_bg = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipseq_control_lambda.bdg"
    params:
        # bw = config["macs2"]["bw"],
        slocal = config["macs2"]["slocal"],
        llocal = config["macs2"]["llocal"],
        qscore = config["macs2"]["fdr"],
        # mfold_low = config["macs2"]["mfold_low"],
        # mfold_high = config["macs2"]["mfold_high"],
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/macs2/macs2_{group}-{species}-{factor}.log"
    shell: """
        (macs2 callpeak --treatment {input.chip_bam} --control {input.input_bam} --format BAMPE --name peakcalling/macs/{wildcards.group}/{wildcards.group}_{wildcards.species}-{wildcards.factor}-chipseq --SPMR --gsize $(faidx {input.fasta} -i chromsizes | awk '{{sum += $2}} END {{print sum}}') --qvalue {params.qscore} --slocal {params.slocal} --llocal {params.llocal} --keep-dup auto --bdg --call-summits) &> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.peaks}) &>> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.summits}) &>> {log}
        """
        # (macs2 callpeak --treatment {input.chip_bam} --control {input.input_bam} --format BAMPE --name peakcalling/macs/{wildcards.group}/{wildcards.group}_{wildcards.species}-{wildcards.factor}-chipseq --SPMR --bw {params.bw} --mfold {params.mfold_low} {params.mfold_high} --gsize $(faidx {input.fasta} -i chromsizes | awk '{{sum += $2}} END {{print sum}}') --qvalue {params.qscore} --slocal {params.slocal} --llocal {params.llocal} --keep-dup auto --bdg --call-summits) &> {log}

rule combine_peaks:
    input:
        cond = "peakcalling/macs/{condition}/{condition}_{species}-{factor}-chipseq_peaks.narrowPeak",
        ctrl = "peakcalling/macs/{control}/{control}_{species}-{factor}-chipseq_peaks.narrowPeak",
    output:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed"
    log:
        "logs/combine_peaks/combine_peaks-{condition}-v-{control}-{species}-{factor}.log"
    shell: """
        (bedtools multiinter -i {input} | \
         bedtools merge -i stdin | \
         sort -k1,1 -k2,2n > {output}) &> {log}
        """

