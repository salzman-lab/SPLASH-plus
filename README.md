# SPLASH+

![Image of SPLASH+](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH_plus.png)
SPLASH+ is a new analytic method to detect a wide range of biological processes that diversify transcripts, including but not limited to RNA splicing, mutations, RNA editing, and V(D)J recombination inference directly on raw sequencing reads by integrating a micro-assembly and biological interpretation framework with the recently developed [SPLASH](https://doi.org/10.1016/j.cell.2023.10.028) algorithm. SPLASH is a unified reference-free algorithm that performs statistical inference directly on raw sequencing reads. SPLASH+ builds on SPLASH by utilizing new approaches to analyze SPLASH’s output, including a new, reference-free statistical approach for de novo assembly (being called as `Compactors`) as well as a framework for interpretation and annotation by assigning a meaningful biological class to each SPLASH's call. 

## How to run SPLASH+
SPLASH+ pipeline consists of 3 main steps:
1. [Running SPLASH](https://github.com/salzman-lab/SPLASH-plus/blob/main/README.md#1--splash): to obtain sequences (`anchors`) that are followed by a set of sample-dependent diverse sequences (`targets`)
2. [Running Compactors](https://github.com/salzman-lab/SPLASH-plus/blob/main/README.md#2--compactors): for de novo local assembly of sequences called by SPLASH
3. [Running Biological Interpretation](https://github.com/salzman-lab/SPLASH-plus/blob/main/README.md#3--biological-interpretation): to assign a biologically relevant event (single base pair change, alternative splicing, ...) accounting for the observed sequence diversity.   

### 1- SPLASH
![Image of SPLASH](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH.png)
SPLASH can be run on an input set of FASTQ files by following the steps in https://github.com/salzman-lab/SPLASH. After running SPLASH, the output file will be a list of significant anchors (we refer to it as `anchors.txt` in this readme), where each anchor is associated with a set of statistically significant sample-dependent target sequences. `anchors.txt` will then be used in the next step (Compactors) to perform a local de novo assembly and obtain extended sequences for each called anchor to facilitate and improve biological interpretation.  

### 2- Compactors
![Image of Compactor](https://github.com/salzman-lab/SPLASH-plus/blob/main/Compactor.png)
Compactors analyze the sequence composition at each position to the right of each seed to evaluate whether the nucleotides presented at that position constitute noise or biological signal. This test is applied recursively on read sets, resulting in one or multiple assembled sequences (compactors) for each called anchor. The compactor step is implemented in a fully containerized Nextflow pipeline (**nf-compactors**) with minimal installation requirements. 

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
For biological interpretation of called anchors (obtained from step 1) using their assembled compactors (obtained from step 2), we provide a script [SPLASH_plus_classification.R](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH_plus_classification.R) to categorize anchors into biologically meaningful events. Currently, we consider 6 different categories: Single base pair changes, alternative splicing, internal splicing (such as insertions, and deletions), 3'UTR, Centromere, and Repeats. The script needs the following inputs:

- `directory`:  Directory for writing output files 
- `compactor_file`: path to the compactors file `compactor_summary.tsv` generated from the compactors step
- `STAR_executable`: path to [STAR](https://github.com/alexdobin/STAR) executable file
- `samtools_executable`: path to [Samtools](https://www.htslib.org/) executable file
- `bedtools_executable`: path to [bedtools](https://bedtools.readthedocs.io/en/latest/) executable file
- `bowtie2_executable`: path to [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) executable file
- `STAR_reference`: path to STAR index files for the reference genome
- `annotated_splice_juncs`: path to the file containing annotated splice junctions from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `annotated_exon_boundaries`: path to the file containing annotated exon boundaries from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `gene_coordinates`: path to the file containing gene coordinates from the reference transcriptome (can be either downloaded or generated from `SPLASH_build.R`)
- `centromere_annotation_file`: (optional) path to the centromere annotation file
- `repeats_annotation_file`: (optional) path to annotation file for repetitive elements
- `UTR_annotation_file`: (optional) path to UTR annotation file
 
The script will generate a file `classified_anchors.tsv` in the same directory specified by the `directory` input argument. The file contains significant anchors along with their compactors, biological classification, and alignment information.

## SPLASH+ output file description

#### Building index and annotation files needed for running classification script 
To be able to run `SPLASH_plus_classification.R` for a reference assembly, you need STAR index for reference genome and three annotation files (`annotated_splice_juncs`, `annotated_exon_boundaries`, `gene_coordinates`) for annotated splice junctions, exons, and genes in the reference transcriptome. 
To build these files, you should obtain a fasta file for the reference genome and a gtf file for the transcriptome annotation. You can then perform the following two steps (note that fasta and gtf files should be from the same assembly as they need to have consistent coordinates, chr names for accurate annotating of anchors):
- **STAR index**: You can use default parameters to build [STAR](https://github.com/alexdobin/STAR) index: 
`STAR --runThreadN 4 --runMode genomeGenerate --genomeDir STAR_index_files --genomeFastaFiles $fasta file$ --sjdbGTFfile $gtf file$`
- **Annotation files**: the three files for annotated exon boundaries, annotated splice junctions, and gene coordinates can be built by running a script we have provided [SPLASH_plus_build.R](https://github.com/salzman-lab/SPLASH-plus/blob/main/SPLASH_plus_build.R). `SPLASH_plus_build.R` needs 3 inputs:
  - `$gtf_file$`: absolute path to the gtf file,
  - `$hisat2_directory$`: directory containing HISAT2 codes downloaded from [HISAT2 repository](https://github.com/DaehwanKimLab/hisat2), the script assumes that there are two python scripts at: `$hisat2_directory$/extract_exons.py` and `$hisat2_directory$/extract_splice_sites.py`),
  - `$outfile_name$`: the name used for the annotation files that script will generate.

  The `SPLASH_plus_build.R` can be run using the following command:  
  `Rscript SPLASH_build.R $gtf_file$ $hisat2_directory$ $outfile_name$`
  If the script finishes successfully, it will generate 3 output annotation files in the same directory as the script: 
  - `$outfile_name$_known_splice_sites.txt` for annotated splice sites (can be used as `annotated_splice_juncs` input for `SPLASH_plus_classification.R`) 
  - `$outfile_name$_exon_coordinates.bed` for annotated exon boundaries (can be used as `annotated_exon_boundaries` input for `SPLASH_plus_classification.R`)
  - `$outfile_name$_genes.bed` for annotated gene coordinates (can be used as `gene_coordinates` input for `SPLASH_plus_classification.R`)

#### Downloading pre-built annotation files for human and mouse genomes:
The human files were built for both [T2T assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_009914755.1/) and [GRCh38 assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/). The mouse files were built based on [mm39 assembly](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/). The annotation files can be downloaded using the following links:
- **Human (T2T)**:
   - `annotated_splice_juncs`: https://drive.google.com/file/d/1owlOQyP1z4cyFvYcAAA-qQmc-K6jGbs9/view?usp=share_link
   - `annotated_exon_boundaries`: https://drive.google.com/file/d/1R-4-ICDAzmIBgQmlOF22nNrCWoSgrmHi/view?usp=share_link
   - `gene_coordinates`: https://drive.google.com/file/d/1L0A7iGXEYiOsPQ0QiJayKPybJ79ZDi2F/view?usp=sharing
- **Human (GRCh38)**:
   - `annotated_splice_juncs`: https://drive.google.com/file/d/1izVHy1m-ddlNgJtFKfWcHdtkc_Y5bHHP/view?usp=sharing
   - `annotated_exon_boundaries`: https://drive.google.com/file/d/1oK6OgQnFFVvybBo0EZ5aIyeoZLAtMyZF/view?usp=sharing
   - `gene_coordinates`: https://drive.google.com/file/d/1REfnl9ZNYcsb-1jSurDHcsL7QFJ00JEp/view?usp=sharing
 - **Mouse (mm39)**:
   - `annotated_splice_juncs`: https://drive.google.com/file/d/1iJhf421nMRDC0uCo_0jh7Nkns8NAieTE/view?usp=sharing
   - `annotated_exon_boundaries`: https://drive.google.com/file/d/1npE0rkxhsDtJk3FeMdfuZwc5Elfuk4bq/view?usp=sharing
   - `gene_coordinates`: https://drive.google.com/file/d/1V8By-yq7AmgXY-XDhipgjjsamL0ghhJa/view?usp=sharing

## Contact
Please contact Roozbeh Dehghannasiri (rdehghan@stanford.edu).

## Citation
Dehghannasiri*, R., Henderson*, G., Bierman, R., Chaung, K., Baharav, T., Wang, P., and Salzman, J. [Unsupervised reference-free inference reveals unrecognized regulated transcriptomic complexity in human single cells](https://www.biorxiv.org/content/10.1101/2022.12.06.519414v2), bioRxiv, (2023).
> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).



