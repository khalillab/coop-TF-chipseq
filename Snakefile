#!/usr/bin/env python

import os
import re
import itertools
from math import log2

configfile: "config.yaml"

FACTOR = config["factor"]

INPUTS = config["igg_samples"]
INPUTS_PASSING = {k:v for k,v in INPUTS.items() if v["pass-qc"]}
CHIPS = config["chip_samples"]
CHIPS_PASSING = {k:v for k,v in CHIPS.items() if v["pass-qc"]}
SAMPLES = {**INPUTS, **CHIPS}
PASSING = {**INPUTS_PASSING, **CHIPS_PASSING}
GROUPS = set(v["group"] for (k,v) in CHIPS.items())

#groups which have at least one passing chip and input sample, so that they are valid for peakcalling
validgroups = set(v["group"] for k,v in CHIPS_PASSING.items() if v["control"] in INPUTS_PASSING)

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))

FIGURES = config["figures"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys()) + ["unmatched"]),
    group = "|".join(set(re.escape(v["group"]) for k,v in CHIPS.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + ["all"])),
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in FIGURES.keys()),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    status = "all|passing",
    norm = "counts|libsizenorm",
    strand = "midpoints|protection",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

def get_samples(status, groups):
    sampledict = {"all": SAMPLES,
                  "passing": PASSING}.get(status)
    if "all" in groups:
        return(list(sampledict.keys()))
    else:
        return([k for k,v in sampledict.items() if v["group"] in groups])

include: "rules/clean_reads.smk"
include: "rules/alignment.smk"
include: "rules/fastqc.smk"
include: "rules/library_processing_summary.smk"
include: "rules/peakcalling.smk"
include: "rules/genome_coverage.smk"
include: "rules/sample_similarity.smk"
include: "rules/datavis.smk"
# include: "rules/mnase-seq_differential_occupancy.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #fastqc
        f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_quality.svg',
        #library processing summaries
        f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.svg",
        f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing-loss.svg",
        #alignment
        expand("alignment/{sample}_{factor}-chipseq-uniquemappers.bam", sample=SAMPLES, factor=FACTOR),
        #peakcalling
        expand("peakcalling/macs/{group}/{group}_{factor}-chipseq_peaks.narrowPeak", group=GROUPS, factor=FACTOR),
        #genome coverage
        expand("coverage/{norm}/{sample}_{factor}-chipseq-{norm}-{readtype}.bw", norm=["counts","libsizenorm"], sample=SAMPLES, readtype=["protection", "midpoints", "midpoints_smoothed"], factor=FACTOR),
        #scatterplots
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}_chipseq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), factor=FACTOR, status=statuscheck(SAMPLES, PASSING), windowsize=config["scatterplot_binsizes"]),
        #datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{factor}}-chipseq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bygroup.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, status=statuscheck(SAMPLES, PASSING), readtype=["protection", "midpoints"], factor=FACTOR) if config["plot_figures"] else []

