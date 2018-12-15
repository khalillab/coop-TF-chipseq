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
    # species = "experimental|spikein",
    read_status = "raw|cleaned|aligned|unaligned",
    figure = "|".join(re.escape(x) for x in FIGURES.keys()),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|libsizenorm",
    strand = "SENSE|ANTISENSE|plus|minus|midpoints|protection",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

# status_norm_sample_dict = {
#     "all":
#         {   "libsizenorm" : SAMPLES,
#             "spikenorm" : SISAMPLES
#         },
#     "passing":
#         {   "libsizenorm" : PASSING,
#             "spikenorm" : SIPASSING
#         }
#     }

# def get_samples(status, norm, groups):
#     if "all" in groups:
#         return(list(status_norm_sample_dict[status][norm].keys()))
#     else:
#         return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

include: "rules/clean_reads.smk"
include: "rules/alignment.smk"
# include: "rules/chip-seq_fastqc.smk"
# include: "rules/chip-seq_library_processing_summary.smk"
# include: "rules/chip-seq_peakcalling.smk"
# include: "rules/chip-seq_genome_coverage.smk"
# include: "rules/chip-seq_sample_similarity.smk"
# include: "rules/chip-seq_datavis.smk"
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
        expand("fastq/cleaned/{sample}_{factor}-chipseq-cleaned.{read}.fastq.gz", sample=SAMPLES, factor=FACTOR, read=["r1", "r2"])
