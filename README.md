# SPLASH+

![Image of SPLASH+](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH_plus.png)

SPLASH+ is a new analytic method to detect a wide range of biological processes that diversify transcripts, including but not limited to RNA splicing, mutations, RNA editing, and V(D)J recombination inference directly on raw sequencing reads by integrating a micro-assembly and biological interpretation framework with the recently developed [SPLASH](https://doi.org/10.1016/j.cell.2023.10.028) algorithm. SPLASH is a unified reference-free algorithm that performs statistical inference directly on raw sequencing reads. SPLASH+ builds on SPLASH by utilizing new approaches to analyze SPLASH’s output, including a new, reference-free statistical approach for de novo assembly (being called as `Compactors`) as well as a framework for interpretation and annotation by assigning a meaningful biological class to each SPLASH's call. 

## How to run SPLASH+
SPLASH+ pipeline consists of 3 main steps:
1. `SPLASH`: to obtain sequences (`anchors`) that are followed by a set of sample-dependent diverse sequences (`targets`)
2. `Compactors`: for de novo local assembly of sequences called by SPLASH
3. `Biological interpretation`: to assign a biologically-relevant event (single base pair change, alternative splicing, ...) accounting for the observed sequence diversity.   

### 1- SPLASH
![Image of SPLASH](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH.png)
SPLASH can be run on an input set of fastq files by following the steps in https://github.com/salzman-lab/SPLASH. After running SPLASH, the output file will be a list of significant anchors (we refer to it as `anchors.txt` in this readme), where each anchor is associated with a set of statistically-significant sample-dependent target sequences. `anchors.txt` will then be used in the next step (Compactors) to perform a local de novo assembly and obtain extended sequences for each called anchor to facilitate and improve biological interpretation.  

### 2- Compactors
![Image of Compactor](https://github.com/salzman-lab/SPLASH-plus/blob/main/Compactor.png)
Compactors tests the sequence composition at each position to the right of each seed to evaluate whether the nucleotides presented at that position constitute noise or biological signal. This test is applied recursively on read sets, resulting in one or multiple assembled sequences (compactors) for each called anchor. The compactor step is implemented in a fully-containerized Nextflow pipeline (**nf-compactors**) with minimal installation requirements. 

Compactors need two input files:
1. `anchors.txt`: a single column file containing the list of significant anchors from SPLASH
2. `samplesheet.csv`: each line in this file provides the path to an input FASTQ file used for running SPLASH.

After running Compactors, two output files will be generated:
1. `compactor_summary.tsv`: Contains the resulting assembled sequences (compactors) for each significant significant anchor. This file will then be used in the next step for biological interpretation.
2. `sample_specificity.tsv`: reporting the supporting read counts for each compactor in input samples.

#### Quick Start for running Compactors pipeline:

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/). You can also use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines.

3. Create your `--fastq_samplesheet`, and run the pipeline. [The FASTQ samplesheet should be of this format](https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv). `anchors_file` can be a any TSV presenting seeds or anchors in a column called `anchor`.

   ```console
   nextflow run salzmanlab/compactors \
       -r main \
       -latest \
       -profile test,YOURPROFILE \
       --fastq_samplesheet samplesheet.csv \
       --anchors_file anchors.txt \
       --outdir <OUTDIR>
   ```
### 3- Biological interpretation




## Contact
Please contact Roozbeh Dehghannasiri (rdehghan@stanford.edu).

## Citation
Dehghannasiri*, R., Henderson*, G., Bierman, R., Chaung, K., Baharav, T., Wang, P., and Salzman, J. [Unsupervised reference-free inference reveals unrecognized regulated transcriptomic complexity in human single cells](https://www.biorxiv.org/content/10.1101/2022.12.06.519414v2), bioRxiv, (2023).
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).



