#!/usr/bin/env python

#   0) cut off the 5'-most base (left from A-tailing)
#   1) search for and remove the constant region of the TruSeq adapter in both reads
#   2) do quality trimming for both reads with --nextseq-trim 2-color quality trimming
rule clean_reads:
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["r1"],
        r2 = lambda wc: SAMPLES[wc.sample]["r2"],
    output:
        r1 = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.r1.fastq.gz",
        r2 = f"fastq/cleaned/{{sample}}_{FACTOR}-chipseq-cleaned.r2.fastq.gz",
        log = "logs/clean_reads/clean_reads-{sample}.log"
    params:
        qual_cutoff = config["cutadapt"]["qual_cutoff"],
        cut_5prime = config["cutadapt"]["cut_5prime"]
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["threads"]
    shell: """
        (cutadapt --cut={params.cut_5prime} -U {params.cut_5prime} --adapter=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --trim-n --cores={threads} --nextseq-trim={params.qual_cutoff} --minimum-length=6 --output={output.r1} --paired-output={output.r2} {input.r1} {input.r2}) &> {output.log}
        """

