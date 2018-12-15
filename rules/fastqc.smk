#!/usr/bin/env python

localrules:
    fastqc_aggregate

#fastQC on raw (demultiplexed) data
rule fastqc:
    input:
        fastq = lambda wc: (config["unmatched"][wc.readnumber] if wc.sample=="unmatched" else SAMPLES[wc.sample][wc.readnumber]) if wc.fqtype=="raw" else f"fastq/{wc.fqtype}/{wc.sample}_{FACTOR}-chipseq-{wc.fqtype}.{wc.readnumber}.fastq.gz"
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}_{readnumber}_fastqc-data-{fqtype}.txt",
    params:
        fname = lambda wc: (re.split('.fq|.fastq', os.path.split(config["unmatched"][wc.readnumber])[1])[0] if wc.sample=="unmatched" else re.split('.fq|.fastq', os.path.split(SAMPLES[wc.sample][wc.readnumber])[1])[0]) if wc.fqtype=="raw" else f"{wc.sample}_{FACTOR}-chipseq-{wc.fqtype}.{wc.readnumber}",
        adapter = lambda wc: {"r1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                              "r2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"}.get(wc.readnumber)
    threads:
        config["threads"]
    log:
        "logs/fastqc/fastqc_{fqtype}-{sample}-{readnumber}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}
        fastqc --adapters <(echo -e "adapter\t{params.adapter}") --nogroup --noextract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} {input.fastq}
        unzip -p qual_ctrl/fastqc/{wildcards.fqtype}/{params.fname}_fastqc.zip {params.fname}_fastqc/fastqc_data.txt > {output}) &> {log}
        """

fastqc_dict = {
        "per_base_qual":
        {   "title" : "Per base sequence quality",
            "fields": "base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus"
            } ,
        "per_tile_qual":
        {   "title" : "Per tile sequence quality",
            "fields": "tile\tbase\tmean\tsample\tstatus"
            },
        "per_seq_qual":
        {   "title" : "Per sequence quality scores",
            "fields": "quality\tcount\tsample\tstatus"
            },
        "per_base_seq_content":
        {   "title" : "Per base sequence content",
            "fields": "base\tg\ta\tt\tc\tsample\tstatus"
            },
        "per_seq_gc":
        {   "title" : "Per sequence GC content",
            "fields": "gc_content\tcount\tsample\tstatus"
            },
        "per_base_n":
        {   "title" : "Per base N content",
            "fields": "base\tn_count\tsample\tstatus"
            },
        "seq_length_dist":
        {   "title" : "Sequence Length Distribution",
            "fields": "length\tcount\tsample\tstatus"
            },
        "seq_duplication":
        {   "title" : "Total Deduplicated Percentage",
            "fields": "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus"
            },
        "adapter_content":
        {   "title" : "Adapter Content",
            "fields": "position\tpct\tsample\tstatus"
            }
        }

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}_{read}_fastqc-data-raw.txt", sample=(["unmatched"]if config["unmatched"]["r1"] and config["unmatched"]["r2"] else []) + list(SAMPLES.keys()), read=["r1", "r2"]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}_{read}_fastqc-data-cleaned.txt", sample=SAMPLES, read=["r1", "r2"]),
        aligned = expand("qual_ctrl/fastqc/aligned/{sample}_{read}_fastqc-data-aligned.txt", sample=SAMPLES, read=["r1", "r2"]),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}_{read}_fastqc-data-unaligned.txt", sample=SAMPLES, read=["r1", "r2"]),
    output:
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_quality.tsv',
        per_tile_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_tile_quality.tsv',
        per_seq_qual =  f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_quality.tsv',
        per_base_seq_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_sequence_content.tsv',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_gc.tsv',
        per_base_n = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_n.tsv',
        seq_length_dist = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_length_distribution.tsv',
        seq_duplication = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_duplication_levels.tsv',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-adapter_content.tsv',
    run:
        shell("rm -f {output}")
        for fastqc_metric, out_path in output.items():
            title = fastqc_dict[fastqc_metric]["title"]
            fields = fastqc_dict[fastqc_metric]["fields"]
            for read_status, read_status_data in input.items():
                sample_id_list = ["_".join(x) for x in itertools.product((["unmatched"]if config["unmatched"]["r1"] and config["unmatched"]["r2"] else []) + list(SAMPLES.keys()), ["r1", "r2"])] if read_status=="raw" else ["_".join(x) for x in itertools.product(SAMPLES.keys(), ["r1", "r2"])]
                for sample_id, fastqc_data in zip(sample_id_list, read_status_data):
                    if sample_id in ["unmatched_r1", "unmatched_r2"] and title=="Adapter Content":
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{title}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{m=$2;for(i=2;i<=NF-2;i++)if($i>m)m=$i; print $1, m, "{sample_id}", "{read_status}"}}' {fastqc_data} | tail -n +2 >> {out_path}""")
                    else:
                        shell("""awk 'BEGIN{{FS=OFS="\t"}} /{title}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, "{sample_id}", "{read_status}"}}' {fastqc_data} | tail -n +2 >> {out_path}""")
            shell("""sed -i "1i {fields}" {out_path}""")

rule plot_fastqc_summary:
    input:
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_quality.tsv',
        per_tile_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_tile_quality.tsv',
        per_seq_qual =  f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_quality.tsv',
        per_base_seq_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_sequence_content.tsv',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_gc.tsv',
        per_base_n = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_n.tsv',
        seq_length_dist = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_length_distribution.tsv',
        seq_duplication = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_duplication_levels.tsv',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-adapter_content.tsv',
    output:
        per_base_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_quality.svg',
        per_tile_qual = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_tile_quality.svg',
        per_seq_qual =  f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_quality.svg',
        per_base_seq_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_sequence_content.svg',
        per_seq_gc = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_sequence_gc.svg',
        per_base_n = f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_n.svg',
        seq_length_dist = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_length_distribution.svg',
        seq_duplication = f'qual_ctrl/fastqc/{FACTOR}_chipseq-sequence_duplication_levels.svg',
        adapter_content = f'qual_ctrl/fastqc/{FACTOR}_chipseq-adapter_content.svg',
    conda:
        "../envs/tidyverse.yaml"
    script:
        "../scripts/fastqc_summary.R"

