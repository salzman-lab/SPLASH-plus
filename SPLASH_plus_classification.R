# the script takes the compactor_summary.tsv file and then classifies each anchor to splicing, 3UTR, mutation events  
# the only input argument is the directory for the compactor_summary.tsv file 

if (!require("stringdist")) {
  install.packages("stringdist", dependencies = TRUE)
  library(stringdist)
}
if (!require("Biostrings")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Biostrings")
  library(Biostrings)
}
if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}
if (!require("GenomicAlignments")) {
  if (!requireNamespace("BiocManager", quietly = TRUE))``
  install.packages("BiocManager")
  BiocManager::install("GenomicAlignments")
  library(GenomicAlignments)
}
if (!require("data.table")) {
  install.packages("data.table", dependencies = TRUE)
  library(data.table)
}

riffle <- function(a, b) {   # this function is used to interleave the elements of two vectors into a vector to create a fasta file
  seqmlab <- seq(length = length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}

##################################################
############ INPUT ARGUMENTS #####################
##################################################
args <- commandArgs(trailingOnly = TRUE)
SPLASH_directory = args[1]                  # Directory for writing output files
compactor_file = args[2]                    # path to the compactors file
STAR_executable = args[3]                   # path to STAR executable file
Samtools_executable = args[4]               # path to Samtools executable file
bedtools_executable = args[5]               # path to Bedtools executable file
Bowtie2_executable = args[6]                # path to Bowtie2 executable file
STAR_reference = args[7]                    # path to STAR index files for the reference genome
annotated_splice_juncs_file = args[8]       # path to the file containing annotated splice junctions
annotated_exon_boundaries_file = args[9]    # path to the file containing annotated exon boundaries
gene_coords_file = args[10]                 # path to the file containing gene coordinates
centromere_annotation_file = args[11]       # path to centromere annotation file 
repeats_annotation_file = args[12]          # path to repeat annotation file 
UTR_annotation_file = args[13]              # path to UTR annotation file 
##################################################
##################################################
##################################################

######### read in compactors ##############
compactors_dt = fread(compactor_file)
###########################################


compactors_dt = compactors_dt[!duplicated(paste(anchor, compactor_majority, sep = "--"))]
compactors_dt[,num_compactor_per_anchor:=length(unique(compactor_majority)), by = anchor]
compactors_dt = compactors_dt[num_compactor_per_anchor > 1]           # Discard anchors with only one compactor
compactors_dt[,compactor_order:=1:num_compactor_per_anchor, by = anchor]
anchor_rank_dt = data.table(unique(compactors_dt$anchor), 1:length(unique(compactors_dt$anchor)))
names(anchor_rank_dt) = c("anchor", "anchor_index")
compactors_dt = merge(compactors_dt, anchor_rank_dt, all.x = TRUE, all.y = FALSE, by.x = "anchor", by.y = "anchor")
compactors_dt_high_rank = compactors_dt[compactor_order<3]
compactors_dt_high_rank[,length_compactor:=nchar(compactor_valid), by = compactor_valid] # As I want to compute distances between the shorter compactors, I need to first consider only that part that is shared between the thwo valid compactors
compactors_dt_high_rank[,min_length_compactor:=min(length_compactor), by = anchor]
compactors_dt_high_rank[,compactor_valid_min_length:=substr(compactor_valid, 1, min_length_compactor), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[,ex_anchor_compactor_valid_min_length :=substr(compactor_valid_min_length, 28, nchar(compactor_valid_min_length)), by = compactor_valid_min_length] # here I try to remove the anchor part of the compactor to compute distance only between the varying part of the compactor

compactors_dt_high_rank_shifted = compactors_dt_high_rank[, data.table::shift(.SD, 1, NA, "lead", TRUE), .SDcols = 1:ncol(compactors_dt_high_rank)]
compactors_dt_high_rank = cbind(compactors_dt_high_rank, compactors_dt_high_rank_shifted[,list(anchor_lead_1, ex_anchor_compactor_valid_min_length_lead_1)])
compactors_dt_high_rank = compactors_dt_high_rank[anchor_lead_1 == anchor]
compactors_dt_high_rank[,lev_dist:=stringdist(ex_anchor_compactor_valid_min_length, ex_anchor_compactor_valid_min_length_lead_1, method = "lv"), by = 1:nrow(compactors_dt_high_rank)] # Levenstein distance between the top two compacotrs of each anchor
compactors_dt_high_rank[,ham_dist:=stringdist(ex_anchor_compactor_valid_min_length, ex_anchor_compactor_valid_min_length_lead_1, method = "hamming"), by = 1:nrow(compactors_dt_high_rank)] # Hamming distance between the top two compactors of each anchor
compactors_dt_high_rank[,lev_operations:=as.character(attributes(adist(ex_anchor_compactor_valid_min_length, ex_anchor_compactor_valid_min_length_lead_1, count=T))$trafos), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[,run_length_D:=max(rle( strsplit(lev_operations,"")[[1]] == 'D')$lengths[which( rle( strsplit(lev_operations,"")[[1]] == 'D')$values=="TRUE")]), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[,run_length_I:=max(rle( strsplit(lev_operations,"")[[1]] == 'I')$lengths[which( rle( strsplit(lev_operations,"")[[1]] == 'I')$values=="TRUE")]), by = 1:nrow(compactors_dt_high_rank)]
compactors_dt_high_rank[is.na(run_length_I), run_length_I:=0]
compactors_dt_high_rank[is.na(run_length_D), run_length_D:=0]
compactors_dt = merge(compactors_dt, unique(compactors_dt_high_rank[, list(anchor, ham_dist, lev_dist, lev_operations, run_length_D, run_length_I)]), all.x = TRUE, all.y = FALSE, by.x = "anchor", by.y = "anchor")

compactors_dt[, anchor_event:=""] # anchor event defines the corresponding event for the anchor

###########################################################################################
# Anchors with the same mu_lev and mu_ham are classified as Base pair change anchors ######
compactors_dt[ham_dist==lev_dist, anchor_event:= paste("Base_pair_change_", ham_dist, sep = "")]
###########################################################################################
###########################################################################################


###########################################################################
############# STAR alignment of the compactors to the genome ##############
###########################################################################

### I first add ranks for the anchors
compactors_dt[,compactor_index:= paste(">", anchor_index, "_", compactor_order, sep = ""), by = 1:nrow(compactors_dt)] # compactor_index is the unique identifier for each compactor sequence

compactors_fasta = riffle(compactors_dt$compactor_index, compactors_dt$compactor_majority)
compactors_fasta  = data.table(compactors_fasta) # the fasta file containg all compactor sequences
write.table(compactors_fasta, paste(SPLASH_directory, "compactors_fasta.fa", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
system(paste(STAR_executable, " --runThreadN 4 --genomeDir  ", STAR_reference, " --readFilesIn ", SPLASH_directory,"compactors_fasta.fa", " --outFileNamePrefix ", SPLASH_directory, "STAR_alignment/Compactors", " --twopassMode Basic --alignIntronMax 1000000 --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 0 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10 --outSAMtype SAM --chimOutType SeparateSAMold --outSAMunmapped None --clip3pAdapterSeq AAAAAAAAA --outSAMattributes NH HI AS nM NM ", sep = ""))

alignment_info_compactors = fread(paste(SPLASH_directory, "STAR_alignment/CompactorsAligned.out.sam", sep = ""), header = FALSE, skip = "NH:") ## now grabbing alignment information for compactor sequences after running STAR
alignment_info_compactors[,V1:= paste(">",V1,sep = ""), by = V1]
alignment_info_compactors[,num_alignments:=.N, by = V1]
alignment_info_compactors[,STAR_num_mismatches:=as.numeric(strsplit(V16, split = ":")[[1]][3]), by = V16]

compactors_dt[,c("STAR_flag", "STAR_chr", "STAR_coord", "STAR_CIGAR", "STAR_num_alignments", "STAR_num_mismatches"):=NULL]
compactors_dt = merge(compactors_dt, alignment_info_compactors[!duplicated(V1), list(V1, V2, V3, V4, V6, num_alignments, STAR_num_mismatches)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "V1")
setnames(compactors_dt,c("V2", "V3", "V4", "V6", "num_alignments"),c("STAR_flag", "STAR_chr", "STAR_coord", "STAR_CIGAR", "STAR_num_alignments"))

### the flags for the STAR alignment status of each compactor
compactors_dt[,compactor_index:=gsub(">", "", compactor_index), by = compactor_index]
compactors_dt[,is.aligned_STAR:=0]   # specifies whether compactor was mapped by STAR
compactors_dt[,is.STAR_chimeric:=0]  # specifies whether STAR reports a chimeric alignment for the compactor
compactors_dt[,is.STAR_SJ:=0]        # specifies whether STAR reports a gapped alignment (i.e., a splice junction for the compactor)

num_chimeric_alignments = data.table()
num_chimeric_alignments = fread(paste(SPLASH_directory, "STAR_alignment/CompactorsLog.final.out", sep = ""), sep = "|", skip = 35)
num_chimeric_alignments[,V2:=gsub("\t", "", V2), by = V2]
if (num_chimeric_alignments$V2[1]!=0){
  chimeric_alignment_info_compactors = fread(paste(SPLASH_directory,"STAR_alignment/CompactorsChimeric.out.sam", sep = ""),header=FALSE,skip="NH:")
  compactors_dt[compactor_index%in%chimeric_alignment_info_compactors$V1, is.STAR_chimeric:=1]
}
compactors_dt[STAR_CIGAR%like%"N", is.STAR_SJ:=1] # the compactors with non NA STAR alignment are flagged as mapped by STAR
compactors_dt[!is.na(STAR_chr) & (is.STAR_chimeric == 0), is.aligned_STAR:=1] # the compactors with non NA STAR alignment are flagged as mapped by STAR

##########################################################################################
##########################################################################################
##########################################################################################


##############################################################################################
########### assigning gene name to each compactor based on its alignment to genome ###########
##############################################################################################
compactors_dt[,compactor_gene:=NULL]
system(paste(Samtools_executable, " view -S -b ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.sam > ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.bam", sep = ""))
system(paste(bedtools_executable, " bamtobed -split -i ",SPLASH_directory,"STAR_alignment/CompactorsAligned.out.bam | sed '/^chr/!d' | sort -k1,1 -k2,2n > ", SPLASH_directory,"STAR_alignment/called_exons.bed", sep = ""))
system(paste(bedtools_executable, " intersect -a ",SPLASH_directory,"STAR_alignment/called_exons.bed -b ", gene_coords_file, " -wb -loj | cut -f 4,10   | ", bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory, "STAR_alignment/compactor_genes.txt", sep = ""))
compactor_genes = fread(paste(SPLASH_directory, "STAR_alignment/compactor_genes.txt", sep = ""), sep = "\t", header = FALSE)
names(compactor_genes) = c("compactor_index", "compactor_gene")
compactors_dt = merge(compactors_dt, compactor_genes[!duplicated(compactor_index)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "compactor_index")
compactors_dt[compactor_gene == ".", compactor_gene:=NA]
compactors_dt[,num_compactor_gene_anchor := length(unique(compactor_gene)), by = anchor]  #number of unique compactor gene names for each anchor (is useful for distinguishing TE-like events in which different compactors might map to multiple genes)
##############################################################################################
##############################################################################################
##############################################################################################

####################################################################################################################
######################### extracting splice junctions for compactors with splice alignment #########################
####################################################################################################################
# I first write sam file that has only those splice alignments that are the best alignmentfor a compactor
alignment_info_compactors = fread(paste(SPLASH_directory, "STAR_alignment/CompactorsAligned.out.sam", sep = ""), header = FALSE, skip = "NH:")
alignment_info_compactors = alignment_info_compactors[!duplicated(V1)][V6%like%"N"]
write.table(alignment_info_compactors, paste(SPLASH_directory, "STAR_alignment/top_splice_alignments.out.sam", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
system(paste(Samtools_executable, " view -H ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.bam > ", SPLASH_directory, "STAR_alignment/sam_header.txt", sep = ""))
# now I need to concatenate the header to the new sam file and then convert the sam file to a bam file and then run bamtobed split to get the extracted splice junctions
system(paste("cat ", SPLASH_directory, "STAR_alignment/sam_header.txt ", SPLASH_directory, "STAR_alignment/top_splice_alignments.out.sam > ", SPLASH_directory, "STAR_alignment/top_splice_alignments_with_header.out.sam" , sep = ""))
system(paste(Samtools_executable, " view -S -b ", SPLASH_directory, "STAR_alignment/top_splice_alignments_with_header.out.sam > ", SPLASH_directory, "STAR_alignment/top_splice_alignments_with_header.out.bam", sep = ""))
system(paste(bedtools_executable, " bamtobed -split -i ", SPLASH_directory, "STAR_alignment/top_splice_alignments_with_header.out.bam > ", SPLASH_directory, "STAR_alignment/extracted_splice_junction.bed", sep = ""))
splice_junctions = fread(paste(SPLASH_directory, "STAR_alignment/extracted_splice_junction.bed", sep = ""), header = FALSE) # read in the extracted splice junctions and then by creating a shifted version of them and then appending them to the original data table columnwise I extract splice junctions as the 2nd coord from the first line and 1st coord from the second line + 1
splice_junctions_shifted = splice_junctions[, shift(.SD, 1, NA, "lead", TRUE), .SDcols=1:6]
splice_junctions = cbind(splice_junctions, splice_junctions_shifted)
splice_junctions = splice_junctions[V4==V4_lead_1] # to subset to only those with the same compactor_index
splice_junctions[,splice_junc:= paste(V1, ":", V3, ":", V2_lead_1+1, sep = ""), by = 1:nrow(splice_junctions)]
splice_junctions[,all_splice_juncs := paste(splice_junc, collapse = "--"), by = V4]

#######################################
##### annotating AS splice sites ######
known_splice_sites = fread(annotated_splice_juncs_file)

# looking at the splice sites with at least two distinct junctions
known_splice_sites[,chr_V2:= paste(V1, V2, sep = ":"), by = paste(V1, V2)]
known_splice_sites[,chr_V3:= paste(V1, V3, sep = ":"), by = paste(V1, V3)]
known_splice_sites[,num_uniq_V2_for_V3:=length(unique(chr_V2)), by = chr_V3] # number of partner splice cites for each V2
known_splice_sites[,num_uniq_V3_for_V2:=length(unique(chr_V3)), by = chr_V2] # number of partner splice cites for each V3
alt_v2 = known_splice_sites[num_uniq_V3_for_V2 > 1]$V2 # the V2 coordinates that have more than one splice site partner
alt_v3 = known_splice_sites[num_uniq_V2_for_V3 > 1]$V3 # the V3 coordinates that have more than one splice site partner
total = c(alt_v2, alt_v3, alt_v2-1, alt_v3-1, alt_v2+1, alt_v3+1) # I concatenate all coordinates with their +-1 counterparts
splice_junctions[,SSA_AS_annot := 0]
splice_junctions[,SSB_AS_annot := 0]
splice_junctions[V3%in%total, SSA_AS_annot := 1]
splice_junctions[V2_lead_1%in%total, SSB_AS_annot := 1]
splice_junctions[,SS_AS_annot:= paste(SSA_AS_annot, ":", SSB_AS_annot, sep = ""), by = 1:nrow(splice_junctions)] # the flag that specifies whether each 5' and 3' splice site is known to be involved in annotated alternative splicing
splice_junctions[,all_SS_AS_annot:= paste(SS_AS_annot, collapse = "--"), by = V4]
#######################################
#######################################

##############################################################################
### now check to see if splice site is an annotated exon boundary ############
known_exon_boundaries = fread(annotated_exon_boundaries_file, header = TRUE, sep = "\t")
known_exon_boundaries = known_exon_boundaries[!duplicated(paste(V1, V2, V3))]
total = c(known_exon_boundaries$chr_V2, known_exon_boundaries$chr_V2_1, known_exon_boundaries$chr_V2_2, known_exon_boundaries$chr_V3, known_exon_boundaries$chr_V3_1, known_exon_boundaries$chr_V3_2)
splice_junctions[,SSA_annot:=0]
splice_junctions[,SSB_annot:=0]
splice_junctions[paste(V1, V3, sep = "")%in%total, SSA_annot:=1]
splice_junctions[paste(V1, V2_lead_1, sep = "")%in%total, SSB_annot:=1]
splice_junctions[,SS_annot := paste(SSA_annot, ":", SSB_annot, sep = ""), by = 1:nrow(splice_junctions)] # the flag that shows whether the 5' and 3' SS of each splice alignment is annotated as an exon boundary
splice_junctions[,all_SS_annot := paste(SS_annot, collapse = "--"), by = V4]
##############################################################################
##############################################################################

compactors_dt[,c("all_splice_juncs", "all_SS_AS_annot", "all_SS_annot") := NULL]
splice_junctions = splice_junctions[!duplicated(V4)]
compactors_dt = merge(compactors_dt, splice_junctions[, list(V4, all_splice_juncs, all_SS_AS_annot, all_SS_annot)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "V4")
##################################################################################################################################
##################################################################################################################################
##################################################################################################################################



##########################################################################################################
######## Bowtie2 alignment of anchors to compactors to remove duplicate anchors ##########################
##########################################################################################################
anchors_fasta = riffle(paste(">", compactors_dt[!duplicated(anchor)]$anchor_index, sep = ""), compactors_dt[!duplicated(anchor)]$anchor)
write.table(anchors_fasta, paste(SPLASH_directory, "anchors_fasta.fa", sep = ""), quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
system(paste("mkdir ",SPLASH_directory, "Bowtie_alignment", sep = "")) # this creates the directory for BOWTIE alignment and index files
system(paste(Bowtie2_executable, "-build ", SPLASH_directory, "compactors_fasta.fa ", SPLASH_directory, "Bowtie_alignment/compactors_bowtie2_index", sep = "")) # making the index files for compactor sequences
#aligning anchors to compactors
system(paste(Bowtie2_executable, " -f -k 10 -x ", SPLASH_directory, "Bowtie_alignment/compactors_bowtie2_index ", "-U ", SPLASH_directory, "anchors_fasta.fa -S ", SPLASH_directory, "Bowtie_alignment/anchors_to_compactors_bowtie.sam", sep = ""))

Bowtie_alignment_info_compactors = fread(paste(SPLASH_directory, "Bowtie_alignment/anchors_to_compactors_bowtie.sam", sep = ""), header = FALSE, fill = TRUE)
Bowtie_alignment_info_compactors = Bowtie_alignment_info_compactors[V12 == "AS:i:0"]
# now I find the corresponding anchor index for all the compactors to which each anchor has been mapped
# and then those that have been mapped to a compactor with higher rank will be flagged as duplicate
Bowtie_alignment_info_compactors[,aligned_compactor_anchor_index := as.numeric(strsplit(V3, split = "_")[[1]][1]), by = V3]
Bowtie_alignment_info_compactors$aligned_compactor_anchor_index = as.numeric(Bowtie_alignment_info_compactors$aligned_compactor_anchor_index)
Bowtie_alignment_info_compactors[,min_aligned_compactor_anchor_index := min(aligned_compactor_anchor_index), by = V1]
Bowtie_alignment_info_compactors[, duplicate_anchor:=0]
Bowtie_alignment_info_compactors[min_aligned_compactor_anchor_index<V1, duplicate_anchor:=1]
Bowtie_alignment_info_compactors[duplicate_anchor==1 & aligned_compactor_anchor_index==min_aligned_compactor_anchor_index, high_ranked_compactors := paste(V3,collapse = ":"), by = V1]
Bowtie_alignment_info_compactors = Bowtie_alignment_info_compactors[duplicate_anchor==1 & !is.na(high_ranked_compactors)]
Bowtie_alignment_info_compactors = Bowtie_alignment_info_compactors[!duplicated(V1)]
compactors_dt = merge(compactors_dt,Bowtie_alignment_info_compactors[,list(V1,duplicate_anchor, high_ranked_compactors)], all.x = TRUE, all.y = FALSE, by.x = "anchor_index", by.y = "V1")
compactors_dt[is.na(duplicate_anchor), duplicate_anchor := 0]
##########################################################################################
##########################################################################################
##########################################################################################


#############################################################################################
################ intersecting with centromere, repeat, and UTR annotation databases #########
#############################################################################################

##### centromere annotation
system(paste(Samtools_executable, " view -S -b ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.sam > ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.bam", sep = ""))
system(paste(bedtools_executable, " bamtobed -split -i ", SPLASH_directory, "STAR_alignment/CompactorsAligned.out.bam | sed '/^chr/!d' | sort -k1,1 -k2,2n > ", SPLASH_directory, "STAR_alignment/called_exons.bed", sep = ""))
system(paste(bedtools_executable, " intersect -a ", SPLASH_directory, "STAR_alignment/called_exons.bed -b ", centromere_annotation_file, " -wb -loj | cut -f 4,10   | sort -nrk1| ", bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory,"STAR_alignment/centromere_annotation.txt", sep = ""))
centromere_annotation = fread(paste(SPLASH_directory, "STAR_alignment/centromere_annotation.txt", sep = ""), sep = "\t", header = FALSE)
names(centromere_annotation) = c("compactor_index", "centromere_annotation")
compactors_dt = merge(compactors_dt,centromere_annotation[!duplicated(compactor_index)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "compactor_index")
compactors_dt[, centromere_annotation:=gsub(".,", "", centromere_annotation, fixed = TRUE), by = centromere_annotation]

##### repeat annotation
system(paste(bedtools_executable, " intersect -a ", SPLASH_directory, "STAR_alignment/called_exons.bed -b ", repeats_annotation_file, " -wb -loj | cut -f 4,10   | sort -nrk1| ", bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory,"STAR_alignment/repeat_annotation.txt", sep = ""))
repeat_annotation = fread(paste(SPLASH_directory, "STAR_alignment/repeat_annotation.txt", sep = ""), sep = "\t", header = FALSE)
names(repeat_annotation) = c("compactor_index", "repeat_annotation")
compactors_dt = merge(compactors_dt, repeat_annotation[!duplicated(compactor_index)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "compactor_index")
compactors_dt[,repeat_annotation:=gsub(".,", "", repeat_annotation, fixed = TRUE), by = repeat_annotation]

##### UTR annotation
system(paste(bedtools_executable, " intersect -a ", SPLASH_directory, "STAR_alignment/called_exons.bed -b ", UTR_annotation_file, " -wb -loj | cut -f 4,10   | sort -nrk1| ", bedtools_executable, " groupby -g 1 -c 2 -o distinct  > ", SPLASH_directory,"STAR_alignment/UTR_annotation.txt", sep = ""))
UTR_annotation = fread(paste(SPLASH_directory, "STAR_alignment/UTR_annotation.txt", sep = ""), sep = "\t", header = FALSE)
names(UTR_annotation) = c("compactor_index", "UTR_annotation")
compactors_dt = merge(compactors_dt,UTR_annotation[!duplicated(compactor_index)], all.x = TRUE, all.y = FALSE, by.x = "compactor_index", by.y = "compactor_index")
compactors_dt[,UTR_annotation:=gsub(".,", "", UTR_annotation, fixed = TRUE), by = UTR_annotation]
##########################################################################################################
##########################################################################################################
##########################################################################################################


#############################################################################
########## Assigning biological classification to each anchor  ##############
#############################################################################

## Each anchor gets one of this categories: Base_pair_change, Splicing, Internal_splicing, 3UTR, Repeat, Centromere 
compactors_dt[,is_STAR_SJ_in_top_two:=0]
compactors_dt[(compactor_order < 3), is_STAR_SJ_in_top_two := sum(is.STAR_SJ), by = anchor]
compactors_dt[,is_STAR_SJ_in_top_two := max(is_STAR_SJ_in_top_two), by = anchor]
compactors_dt[,is_aligned_STAR_in_top_two := 0]
compactors_dt[(compactor_order < 3), is_aligned_STAR_in_top_two := sum(is.aligned_STAR), by = anchor]
compactors_dt[,is_aligned_STAR_in_top_two := max(is_aligned_STAR_in_top_two), by = anchor]
compactors_dt[ is_STAR_SJ_in_top_two>0 & (ham_dist!=lev_dist | ham_dist > 5) & num_compactor_gene_anchor == 1, anchor_event := "Splicing"] # if both top compactors are mapped as SJ by STAR we keep event as splicing
compactors_dt[(lev_dist < run_length_D + run_length_I+1) & (run_length_D > 1 | run_length_I>1), anchor_event := "Internal_splicing"]
compactors_dt[ ,help:=2*max(run_length_D, run_length_I)+1,by = 1 : nrow(compactors_dt)]
compactors_dt[ lev_dist < help & (run_length_D > 1 | run_length_I > 1), anchor_event := "Internal_splicing"]
compactors_dt[, help:=NULL]
compactors_dt[!is.na(STAR_CIGAR), num_STAR_match := sum(explodeCigarOpLengths(STAR_CIGAR, ops = c("M"))[[1]]), by = STAR_CIGAR]

repeat_anchors = unique(compactors_dt[anchor_event == ""  & repeat_annotation != "." & (compactor_order == 1 | compactor_order == 2)]$anchor_index)
UTR_anchors = unique(compactors_dt[anchor_event == ""  & UTR_annotation!="." & (compactor_order == 1 | compactor_order == 2)]$anchor_index)
centromere_anchors = unique(compactors_dt[anchor_event == "" & centromere_annotation != "." & (compactor_order == 1 | compactor_order == 2)]$anchor_index)

## below we assign Repeat, Centormere, and 3UTR classes to those anchors that had not been already classified and mapped to an annotated element
## if an anchor maps to an element in more than one class, the class is assigned according to this priority: 3UTR > Centromere > Repeat
compactors_dt[anchor_index %in% repeat_anchors, anchor_event := "Repeat"] 
compactors_dt[anchor_index %in% centromere_anchors, anchor_event := "Centromere"]
compactors_dt[anchor_index %in% UTR_anchors, anchor_event := "3UTR"]
##################################################################################
##################################################################################

### writing out the final output file containing classified anchors and their compactors
write.table(compactors_dt, paste(directory, "classified_anchors.tsv", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)