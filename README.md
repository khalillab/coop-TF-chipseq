
# ChIP-seq analysis pipeline, [Bragdon *et. al.* 2022](https://doi.org/10.1101/2022.05.22.492993)

[Snakemake](https://snakemake.github.io/) workflow used to analyze ChIP-seq data for the 2022 publication [*Cooperative assembly confers regulatory specificity and long-term genetic circuit stability*](https://doi.org/10.1101/2022.05.22.492993). For the pipeline with raw data, see the Zenodo archive (coming).

## workflow summary
The workflow has the following major steps:

- adapter and quality trimming with [cutadapt](http://cutadapt.readthedocs.io/en/stable/guide.html)
- alignment with [bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
- summaries of quality statistics from [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- summaries of library processing statistics
- peakcalling with [MACS2](https://github.com/taoliu/MACS) and [IDR](https://github.com/nboley/idr)
- generation of coverage tracks representing fragment midpoints, and fragment protection
- library size and spike-in normalization of coverage
- genome-wide scatterplots and correlations
- data visualization

## requirements

### required software

- Unix-like operating system (tested on CentOS 7.2.1511)
- Git
- [conda](https://conda.io/docs/user-guide/install/index.html)

### required files

- Paired FASTQ files of demultiplexed ChIP-seq libraries.

- FASTA files of the experimental and spike-in genomes.

- [BED6](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) format annotation files:
    - optional: annotations for data visualization

## instructions

TODO
