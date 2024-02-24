# SPLASH+

![Image of SPLASH+](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH_plus.png)

SPLASH+ is a new analytic method to detect a wide range of biological processes that diversify transcripts, including but not limited to RNA splicing, mutations, RNA editing, and V(D)J recombination inference directly on raw sequencing reads by integrating a micro-assembly and biological interpretation framework with the recently developed [SPLASH](https://doi.org/10.1016/j.cell.2023.10.028) algorithm. SPLASH is a unified reference-free algorithm that performs statistical inference directly on raw sequencing reads. SPLASH+ builds on SPLASH by utilizing new approaches to analyze SPLASH’s output, including a new, reference-free statistical approach for de novo assembly (being called as `Compactors`) as well as a framework for interpretation and annotation by assigning a meaningful biological class to each SPLASH's call. 

## How to run SPLASH+
SPLASH+ pipeline consists of 3 main steps:
1. `SPLASH`: to obtain sequences (`anchors`) that are followed by a set of sample-dependent diverse sequences (`targets`)
2. `Compactors`: for de novo local assembly of sequences called by SPLASH
3. `Biological interpretation`: to assign a biologically-relevant event (single base pair change, alternative splicing, ...) accounting for the observed sequence diversity.   

## Running SPLASH
![Image of SPLASH](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH.png)

## Running Compactors
![Image of Compactor](https://github.com/salzman-lab/SPLASH-plus/blob/main/Compactor.png)
**nf-compactors** is a pipeline designed to run the seed-based assembly tool, compactors. 

Given a set of seeds (i.e. 'anchors' from SPLASH) and a FASTQ list, compactors analyzes the sequence composition directly to the right of the seeds across the raw FASTQs. Compactors tests the sequence composition at each position to the right of each seed to evaluate whether the nucleotides presented at that position constitute noise or biological signal. This test is applied recursively on read sets, resulting in one or multiple assembled sequences (compactors) for each seed. The resulting assembled sequences are reported with their supporting read counts in `compactor_summary.tsv` and their read counts in individual samples are presented in `sample_specificity.tsv`.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules).


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.10.3`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

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

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.





## Contact
Please contact Roozbeh Dehghannasiri (rdehghan@stanford.edu).

## Citation
Dehghannasiri*, R., Henderson*, G., Bierman, R., Chaung, K., Baharav, T., Wang, P., and Salzman, J. [Unsupervised reference-free inference reveals unrecognized regulated transcriptomic complexity in human single cells](https://www.biorxiv.org/content/10.1101/2022.12.06.519414v2), bioRxiv, (2023).
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).



