#!/usr/bin/env python

import os
import re
import itertools
from math import log2

configfile: "config.yaml"

FACTOR = config["factor"]
FIGURES = config["figures"]

INPUTS = config["input_samples"]
CHIPS = config["chip_samples"]
SAMPLES = {**INPUTS, **CHIPS}

def get_samples(search_dict=CHIPS,
                passing=False,
                spikein=False,
                paired=False,
                groups=None):
    if passing:
        search_dict = {k:v for k,v in search_dict.items() if v["pass-qc"]}
    if spikein:
        search_dict = {k:v for k,v in search_dict.items() if v["spikein"]}
    if paired:
        search_dict = {k:v for k,v in search_dict.items() \
                if v["control"] in get_samples(search_dict=INPUTS,
                                               passing=passing,
                                               spikein=spikein,
                                               paired=False)}
    if groups and "all" not in groups:
        search_dict = {k:v for k,v in search_dict.items() if v["group"] in groups}
    return search_dict

def statuscheck(dict1, dict2):
    if dict1 == dict2:
        if len(dict1) == 0:
            return []
        else:
            return ["passing"]
    else:
        return ["passing", "all"]

def conditioncheck(conditionlist):
    if len(conditionlist) == 0:
        return []
    elif len(conditionlist) == 1:
        return conditionlist
    else:
        return conditionlist + ["all"]

SISAMPLES = get_samples(search_dict=SAMPLES, spikein=True)

allgroups = [v["group"] for k,v in get_samples(passing=True, paired=True).items()]
allgroups_si = [v["group"] for k,v in get_samples(passing=True, spikein=True, paired=True).items()]
#groups with >= 2 passing and paired samples, so that they are valid for peakcalling and diff binding
validgroups = set(z for z in allgroups if allgroups.count(z)>=2)
validgroups_si = set(z for z in allgroups_si if allgroups_si.count(z)>=2)

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups]))

comparisons_si =  config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"] if list(d.keys())[0] and list(d.values())[0] in validgroups_si]))

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys()) + ["unmatched"]),
    group = "|".join(set(re.escape(v["group"]) for k,v in CHIPS.items())),
    control = "|".join(set(re.escape(x) for x in controlgroups + (conditiongroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + (controlgroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in FIGURES.keys()),
    # annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    # annotation = "|".join(re.escape(x) for x in set(list(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES])) + list(config["differential_occupancy"]["annotations"].keys() if config["differential_occupancy"]["annotations"] else []) + ["peaks"])),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    readtype = "|".join(list(itertools.chain.from_iterable([[x, x+"-input-subtracted"] for x in ["plus", "minus", "midpoints", "protection"]]))),
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

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

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
        #fastqc
        f'qual_ctrl/fastqc/{FACTOR}_chipseq-per_base_quality.svg',
        #library processing summaries
        f"qual_ctrl/fragment_length_distributions/{FACTOR}_chipseq_fragment_length_distributions.svg",
        f"qual_ctrl/read_processing/{FACTOR}_chipseq_read_processing-loss.svg",
        expand(f"qual_ctrl/spikein/{FACTOR}-chipseq_spikein-plots-{{status}}.svg", status=statuscheck(get_samples(spikein=True, paired=True), get_samples(passing=True, spikein=True, paired=True))),
        #alignment
        expand(f"alignment/{{sample}}_{FACTOR}-chipseq-uniquemappers.bam", sample=SAMPLES),
        #peakcalling
        expand(f"peakcalling/sample_peaks/{{sample}}_experimental-{FACTOR}-chipseq_peaks.narrowPeak", sample=get_samples(passing=True, paired=True)),
        expand(f"peakcalling/sample_peaks/{{sample}}_spikein-{FACTOR}-chipseq_peaks.narrowPeak", sample=get_samples(passing=True, spikein=True, paired=True)),
        expand(f"peakcalling/{{group}}/{{group}}_experimental-{FACTOR}-chipseq-idrpeaks.narrowPeak", group=validgroups),
        expand(f"peakcalling/{{group}}/{{group}}_spikein-{FACTOR}-chipseq-idrpeaks.narrowPeak", group=validgroups_si),
        #genome coverage
        expand(f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw", norm=["counts","libsizenorm"], sample=SAMPLES, readtype=["protection", "midpoints", "midpoints_smoothed"]),
        expand(f"coverage/{{norm}}/{{sample}}_{FACTOR}-chipseq-{{norm}}-{{readtype}}.bw", norm=["sicounts","spikenorm"], sample=SISAMPLES, readtype=["protection", "midpoints", "midpoints_smoothed"]),
        expand(f"coverage/libsizenorm/{{sample}}_{FACTOR}-chipseq-libsizenorm-{{readtype}}.bw", sample=get_samples(paired=True), readtype=["protection-input-subtracted", "midpoints-input-subtracted", "midpoints-input-subtracted_smoothed"]),
        expand(f"coverage/spikenorm/{{sample}}_{FACTOR}-chipseq-spikenorm-{{readtype}}.bw", sample=get_samples(spikein=True, paired=True), readtype=["protection-input-subtracted", "midpoints-input-subtracted", "midpoints-input-subtracted_smoothed"]),
        #scatterplots
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}_chipseq-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), factor=FACTOR, status=statuscheck(SAMPLES, get_samples(search_dict=SAMPLES, passing=True)), windowsize=config["scatterplot_binsizes"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}_chipseq-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), factor=FACTOR, status=statuscheck(SISAMPLES, get_samples(search_dict=SISAMPLES, passing=True)), windowsize=config["scatterplot_binsizes"]) if comparisons_si else [],
        #datavis
        # expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{factor}}-chipseq_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, status=statuscheck(SAMPLES, get_samples(search_dict=SAMPLES, passing=True)), readtype=["midpoints", "midpoints-input-subtracted", "protection", "protection-input-subtracted"], factor=FACTOR) if config["plot_figures"] else [],
        # expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{readtype}}/{{factor}}-chipseq_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_{{readtype}}-heatmap-bysample.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, status=statuscheck(SISAMPLES, SIPASSING), readtype=["midpoints", "midpoints-input-subtracted", "protection", "protection-input-subtracted"], factor=FACTOR) if comparisons_si and config["plot_figures"] else [],
        # differential binding
        # expand(f"diff_binding/peaks/{{condition}}-v-{{control}}/{{condition}}-v-{{control}}_experimental-{FACTOR}-peaks.bed", zip, condition=conditiongroups, control=controlgroups)

