#!/bin/bash

snakemake -pr \
    -R `cat <(snakemake --lc --rerun-incomplete) \
            <(snakemake --li --rerun-incomplete) \
            <(snakemake --lp --rerun-incomplete) | sort -u` \
    --rerun-incomplete \
    --cores 1 \
    --use-conda \
    --conda-prefix ../conda
