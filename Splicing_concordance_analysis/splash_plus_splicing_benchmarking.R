#### this script takes the splicing results for SPLASH+, leafcutter and SpliZ on Tabula Sapiens data and computes the concordance of the calls by each method for 
### each tissue and finally generates the concordance heatmap Figure 3B in the paper (the last code block)
library(data.table)
library(ggplot2)

##### read in leafcutter output files for each donor tissue ######
leafcutter_TSP1_blood = fread("Leafcutter_TSP1-Blood.txt")
leafcutter_TSP1_lung = fread("Leafcutter_TSP1-Lung.txt")
leafcutter_TSP1_muscle = fread("Leafcutter_TSP1-Muscle.txt")
leafcutter_TSP2_lung = fread("Leafcutter_TSP2-Lung.txt")
leafcutter_TSP2_muscle = fread("Leafcutter_TSP2-Muscle.txt")
leafcutter_TSP7_blood = fread("Leafcutter_TSP7-Blood.txt")
leafcutter_TSP4_muscle = fread("Leafcutter_TSP4-Muscle.txt")
###################################################################

## read in splash plus files, each file contains all the splicing genes called by splash plus for all the donors containing that tissue #### 
SPLASH_plus_muscle = fread("SPLASH_plus_muscle.txt")
SPLASH_plus_lung = fread("SPLASH_plus_lung.txt")
SPLASH_plus_blood = fread("SPLASH_plus_blood.txt")
#################################################################################

### read in SpliZ files for each donor, each file contains the splicing genes called by spliz for across all tissues within each donor ############
spliz_TSP1 = fread("SpliZ_TSP1_SS2.tsv")
spliz_TSP2 = fread("SpliZ_TSP2_SS2.tsv")
spliz_TSP4 = fread("SpliZ_TSP4_SS2.tsv")
spliz_TSP7 = fread("SpliZ_TSP7_SS2.tsv")
##################################################################################


leafcutter_TSP1_blood[,donor:="TSP1"]
leafcutter_TSP1_blood[,tissue:="Blood"]
leafcutter_TSP1_lung[,donor:="TSP1"]
leafcutter_TSP1_lung[,tissue:="Lung"]
leafcutter_TSP1_muscle[,donor:="TSP1"]
leafcutter_TSP1_muscle[,tissue:="Muscle"]
leafcutter_TSP2_lung[,donor:="TSP2"]
leafcutter_TSP2_lung[,tissue:="Lung"]
leafcutter_TSP2_muscle[,donor:="TSP2"]
leafcutter_TSP2_muscle[,tissue:="Muscle"]
leafcutter_TSP7_blood[,donor:="TSP7"]
leafcutter_TSP7_blood[,tissue:="Blood"]
leafcutter_TSP4_muscle[,donor:="TSP4"]
leafcutter_TSP4_muscle[,tissue:="Muscle"]
leafcutter = rbind(leafcutter_TSP1_blood, leafcutter_TSP1_lung, leafcutter_TSP1_muscle, leafcutter_TSP2_lung, leafcutter_TSP2_muscle, leafcutter_TSP7_blood, leafcutter_TSP4_muscle)
leafcutter_called = leafcutter[abs(logef)>1.5 & p.adjust<0.05] # intron clusters with p-val <0.05 and effectsize >1.5 are called


spliz_TSP1[,donor:="TSP1"]
spliz_TSP1 = spliz_TSP1[perm_pval_adj_scZ<0.05]
spliz_TSP1 = unique(spliz_TSP1[,list(geneR1A_uniq,tissue,donor)]) 
names(spliz_TSP1) = c("gene","grouping_level_1","donor")
spliz_TSP2[,donor:="TSP2"]
spliz_TSP2 = spliz_TSP2[perm_pval_adj_scZ<0.05]
spliz_TSP2 = unique(spliz_TSP2[,list(geneR1A_uniq,tissue,donor)]) 
names(spliz_TSP2) = c("gene","grouping_level_1","donor")
spliz_TSP4[,donor:="TSP4"]
spliz_TSP4 = spliz_TSP4[scZ_pval<0.05]
spliz_TSP7[,donor:="TSP7"]
spliz = rbind(spliz_TSP1,spliz_TSP2,unique(spliz_TSP4[,list(gene,grouping_level_1,donor)]),unique(spliz_TSP7[,list(gene,grouping_level_1,donor)]))

#################################
#################################
### COMPARISON FOR MUSCLE #######
#################################
#################################
Leafcutter_TSP1 = data.table(unique(leafcutter_called[tissue=="Muscle" & donor=="TSP1"]$genes))
Leafcutter_TSP1 = Leafcutter_TSP1[!is.na(V1)] # I want to remove introns from unknown gene names
Leafcutter_TSP2 = data.table(unique(leafcutter_called[tissue=="Muscle" & donor=="TSP2"]$genes))
Leafcutter_TSP2 = Leafcutter_TSP2[!is.na(V1)] # I want to remove introns from unknown gene names
Leafcutter_TSP4 = data.table(unique(leafcutter_called[tissue=="Muscle" & donor=="TSP4"]$genes))
Leafcutter_TSP4 = Leafcutter_TSP4[!is.na(V1)] # I want to remove introns from unknown gene names
SPLASH_plus_TSP1 = data.table(SPLASH_plus_muscle[donor=="TSP1"][!duplicated(Gene)]$Gene)
SPLASH_plus_TSP2 = data.table(SPLASH_plus_muscle[donor=="TSP2"][!duplicated(Gene)]$Gene)
SPLASH_plus_TSP4 = data.table(SPLASH_plus_muscle[donor=="TSP4"][!duplicated(Gene)]$Gene)

Leafcutter_TSP1_transformed <- Leafcutter_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP1)] # I do this to get separate names for introns with multiple genes with overlapping names so that for an intron with annotated name "Gene1,Gene2", I get "Gene1" and "Gene2"  
Leafcutter_TSP2_transformed <- Leafcutter_TSP2[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP2)]
Leafcutter_TSP4_transformed <- Leafcutter_TSP4[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP4)]
SPLASH_plus_TSP1_transformed <- SPLASH_plus_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP1)]
SPLASH_plus_TSP2_transformed <- SPLASH_plus_TSP2[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP2)]
SPLASH_plus_TSP4_transformed <- SPLASH_plus_TSP4[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP4)]

Leafcutter_TSP1_transformed[,found_in_Leafcutter_TSP2:=0]
Leafcutter_TSP1_transformed[,found_in_Leafcutter_TSP4:=0]
Leafcutter_TSP1_transformed[,found_in_SPLASH_plus_TSP1:=0]
Leafcutter_TSP1_transformed[gene%in%Leafcutter_TSP2_transformed$gene,found_in_Leafcutter_TSP2:=1]
Leafcutter_TSP1_transformed[gene%in%Leafcutter_TSP4_transformed$gene,found_in_Leafcutter_TSP4:=1]
Leafcutter_TSP1_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
Leafcutter_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP2+found_in_Leafcutter_TSP4+found_in_SPLASH_plus_TSP1]
setorder(Leafcutter_TSP1_transformed, nrow, -sum_found, gene)
Leafcutter_TSP1_transformed = Leafcutter_TSP1_transformed[!duplicated(nrow)]

Leafcutter_TSP2_transformed[,found_in_Leafcutter_TSP1:=0]
Leafcutter_TSP2_transformed[,found_in_Leafcutter_TSP4:=0]
Leafcutter_TSP2_transformed[,found_in_SPLASH_plus_TSP2:=0]
Leafcutter_TSP2_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
Leafcutter_TSP2_transformed[gene%in%Leafcutter_TSP4_transformed$gene,found_in_Leafcutter_TSP4:=1]
Leafcutter_TSP2_transformed[gene%in%SPLASH_plus_TSP2_transformed$gene,found_in_SPLASH_plus_TSP2:=1]
Leafcutter_TSP2_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_Leafcutter_TSP4+found_in_SPLASH_plus_TSP2]
setorder(Leafcutter_TSP2_transformed, nrow, -sum_found, gene)
Leafcutter_TSP2_transformed = Leafcutter_TSP2_transformed[!duplicated(nrow)]

Leafcutter_TSP4_transformed[,found_in_Leafcutter_TSP1:=0]
Leafcutter_TSP4_transformed[,found_in_Leafcutter_TSP2:=0]
Leafcutter_TSP4_transformed[,found_in_SPLASH_plus_TSP4:=0]
Leafcutter_TSP4_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
Leafcutter_TSP4_transformed[gene%in%Leafcutter_TSP2_transformed$gene,found_in_Leafcutter_TSP2:=1]
Leafcutter_TSP4_transformed[gene%in%SPLASH_plus_TSP4_transformed$gene,found_in_SPLASH_plus_TSP4:=1]
Leafcutter_TSP4_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_Leafcutter_TSP2+found_in_SPLASH_plus_TSP4]
setorder(Leafcutter_TSP4_transformed, nrow, -sum_found, gene)
Leafcutter_TSP4_transformed = Leafcutter_TSP4_transformed[!duplicated(nrow)]

SPLASH_plus_TSP1_transformed[,found_in_Leafcutter_TSP1:=0]
SPLASH_plus_TSP1_transformed[,found_in_SPLASH_plus_TSP2:=0]
SPLASH_plus_TSP1_transformed[,found_in_SPLASH_plus_TSP4:=0]
SPLASH_plus_TSP1_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
SPLASH_plus_TSP1_transformed[gene%in%SPLASH_plus_TSP2_transformed$gene,found_in_SPLASH_plus_TSP2:=1]
SPLASH_plus_TSP1_transformed[gene%in%SPLASH_plus_TSP4_transformed$gene,found_in_SPLASH_plus_TSP4:=1]
SPLASH_plus_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_SPLASH_plus_TSP2+found_in_SPLASH_plus_TSP4]
setorder(SPLASH_plus_TSP1_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP1_transformed = SPLASH_plus_TSP1_transformed[!duplicated(nrow)]

SPLASH_plus_TSP2_transformed[,found_in_Leafcutter_TSP2:=0]
SPLASH_plus_TSP2_transformed[,found_in_SPLASH_plus_TSP1:=0]
SPLASH_plus_TSP2_transformed[,found_in_SPLASH_plus_TSP4:=0]
SPLASH_plus_TSP2_transformed[gene%in%Leafcutter_TSP2_transformed$gene,found_in_Leafcutter_TSP2:=1]
SPLASH_plus_TSP2_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
SPLASH_plus_TSP2_transformed[gene%in%SPLASH_plus_TSP4_transformed$gene,found_in_SPLASH_plus_TSP4:=1]
SPLASH_plus_TSP2_transformed[,sum_found:=found_in_Leafcutter_TSP2+found_in_SPLASH_plus_TSP1+found_in_SPLASH_plus_TSP4]
setorder(SPLASH_plus_TSP2_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP2_transformed = SPLASH_plus_TSP2_transformed[!duplicated(nrow)]

SPLASH_plus_TSP4_transformed[,found_in_Leafcutter_TSP4:=0]
SPLASH_plus_TSP4_transformed[,found_in_SPLASH_plus_TSP1:=0]
SPLASH_plus_TSP4_transformed[,found_in_SPLASH_plus_TSP2:=0]
SPLASH_plus_TSP4_transformed[gene%in%Leafcutter_TSP4_transformed$gene,found_in_Leafcutter_TSP4:=1]
SPLASH_plus_TSP4_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
SPLASH_plus_TSP4_transformed[gene%in%SPLASH_plus_TSP2_transformed$gene,found_in_SPLASH_plus_TSP2:=1]
SPLASH_plus_TSP4_transformed[,sum_found:=found_in_Leafcutter_TSP4+found_in_SPLASH_plus_TSP1+found_in_SPLASH_plus_TSP2]
setorder(SPLASH_plus_TSP4_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP4_transformed = SPLASH_plus_TSP4_transformed[!duplicated(nrow)]

spliz_muscle = spliz[grouping_level_1=="Muscle"]
muscle_splicing_genes_all_methods = list(Leafcutter_TSP1 =  unique(Leafcutter_TSP1_transformed$gene), Leafcutter_TSP2 = unique(Leafcutter_TSP2_transformed$gene), Leafcutter_TSP4 = unique(Leafcutter_TSP4_transformed$gene), SPLASH_plus_TSP1 = unique(SPLASH_plus_TSP1_transformed$gene), SPLASH_plus_TSP2 = unique(SPLASH_plus_TSP2_transformed$gene), SPLASH_plus_TSP4 = unique(SPLASH_plus_TSP4_transformed$gene), SpliZ_TSP1 = spliz_muscle[!gene%like%"unknown" & donor=="TSP1"]$gene, SpliZ_TSP2 = spliz_muscle[!gene%like%"unknown" & donor=="TSP2"]$gene,SpliZ_TSP4 = spliz_muscle[!gene%like%"unknown" & donor=="TSP4"]$gene)

#################################
#################################
### COMPARISON FOR LUNG #########
#################################
#################################
Leafcutter_TSP1 = data.table(unique(leafcutter_called[tissue=="Lung" & donor=="TSP1"]$genes))
Leafcutter_TSP1 = Leafcutter_TSP1[!is.na(V1)] # I want to remove unknown gene names
Leafcutter_TSP2 = data.table(unique(leafcutter_called[tissue=="Lung" & donor=="TSP2"]$genes))
Leafcutter_TSP2 = Leafcutter_TSP2[!is.na(V1)] # I want to remove unknown gene names
SPLASH_plus_TSP1 = data.table(SPLASH_plus_lung[donor=="TSP1"][!duplicated(Gene)]$Gene)
SPLASH_plus_TSP2 = data.table(SPLASH_plus_lung[donor=="TSP2"][!duplicated(Gene)]$Gene)

Leafcutter_TSP1_transformed <- Leafcutter_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP1)]
Leafcutter_TSP2_transformed <- Leafcutter_TSP2[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP2)]
SPLASH_plus_TSP1_transformed <- SPLASH_plus_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP1)]
SPLASH_plus_TSP2_transformed <- SPLASH_plus_TSP2[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP2)]

Leafcutter_TSP1_transformed[,found_in_Leafcutter_TSP2:=0]
Leafcutter_TSP1_transformed[,found_in_SPLASH_plus_TSP1:=0]
Leafcutter_TSP1_transformed[gene%in%Leafcutter_TSP2_transformed$gene,found_in_Leafcutter_TSP2:=1]
Leafcutter_TSP1_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
Leafcutter_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP2+found_in_SPLASH_plus_TSP1]
setorder(Leafcutter_TSP1_transformed, nrow, -sum_found, gene)
Leafcutter_TSP1_transformed = Leafcutter_TSP1_transformed[!duplicated(nrow)]

Leafcutter_TSP2_transformed[,found_in_Leafcutter_TSP1:=0]
Leafcutter_TSP2_transformed[,found_in_SPLASH_plus_TSP2:=0]
Leafcutter_TSP2_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
Leafcutter_TSP2_transformed[gene%in%SPLASH_plus_TSP2_transformed$gene,found_in_SPLASH_plus_TSP2:=1]
Leafcutter_TSP2_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_SPLASH_plus_TSP2]
setorder(Leafcutter_TSP2_transformed, nrow, -sum_found, gene)
Leafcutter_TSP2_transformed = Leafcutter_TSP2_transformed[!duplicated(nrow)]

SPLASH_plus_TSP1_transformed[,found_in_Leafcutter_TSP1:=0]
SPLASH_plus_TSP1_transformed[,found_in_SPLASH_plus_TSP2:=0]
SPLASH_plus_TSP1_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
SPLASH_plus_TSP1_transformed[gene%in%SPLASH_plus_TSP2_transformed$gene,found_in_SPLASH_plus_TSP2:=1]
SPLASH_plus_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_SPLASH_plus_TSP2]
setorder(SPLASH_plus_TSP1_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP1_transformed = SPLASH_plus_TSP1_transformed[!duplicated(nrow)]

SPLASH_plus_TSP2_transformed[,found_in_Leafcutter_TSP2:=0]
SPLASH_plus_TSP2_transformed[,found_in_SPLASH_plus_TSP1:=0]
SPLASH_plus_TSP2_transformed[gene%in%Leafcutter_TSP2_transformed$gene,found_in_Leafcutter_TSP2:=1]
SPLASH_plus_TSP2_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
SPLASH_plus_TSP2_transformed[,sum_found:=found_in_Leafcutter_TSP2+found_in_SPLASH_plus_TSP1]
setorder(SPLASH_plus_TSP2_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP2_transformed = SPLASH_plus_TSP2_transformed[!duplicated(nrow)]

spliz_lung = spliz[grouping_level_1=="Lung"]
lung_splicing_genes_all_methods = list(Leafcutter_TSP1 =  unique(Leafcutter_TSP1_transformed$gene), Leafcutter_TSP2 = unique(Leafcutter_TSP2_transformed$gene), SPLASH_plus_TSP1 = unique(SPLASH_plus_TSP1_transformed$gene), SPLASH_plus_TSP2 = unique(SPLASH_plus_TSP2_transformed$gene), SpliZ_TSP1 = spliz_lung[!gene%like%"unknown" & donor=="TSP1"]$gene, SpliZ_TSP2 = spliz_lung[!gene%like%"unknown" & donor=="TSP2"]$gene)

#################################
#################################
### COMPARISON FOR BLOOD  #######
#################################
#################################
spliz_blood = spliz[grouping_level_1=="Blood"]

Leafcutter_TSP1 = data.table(unique(leafcutter_called[tissue=="Blood" & donor=="TSP1"]$genes))
Leafcutter_TSP1 = Leafcutter_TSP1[!is.na(V1)] # I want to remove unknown gene names
Leafcutter_TSP7 = data.table(unique(leafcutter_called[tissue=="Blood" & donor=="TSP7"]$genes))
Leafcutter_TSP7 = Leafcutter_TSP7[!is.na(V1)] # I want to remove unknown gene names
SPLASH_plus_TSP1 = data.table(SPLASH_plus_blood[donor=="TSP1"][!duplicated(Gene)]$Gene)
SPLASH_plus_TSP7 = data.table(SPLASH_plus_blood[donor=="TSP7"][!duplicated(Gene)]$Gene)

Leafcutter_TSP1_transformed <- Leafcutter_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP1)]
Leafcutter_TSP7_transformed <- Leafcutter_TSP7[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(Leafcutter_TSP7)]
SPLASH_plus_TSP1_transformed <- SPLASH_plus_TSP1[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP1)]
SPLASH_plus_TSP7_transformed <- SPLASH_plus_TSP7[, .(gene = unlist(lapply(strsplit(V1, ","), sort))), by = 1:nrow(SPLASH_plus_TSP7)]

Leafcutter_TSP1_transformed[,found_in_Leafcutter_TSP7:=0]
Leafcutter_TSP1_transformed[,found_in_SPLASH_plus_TSP1:=0]
Leafcutter_TSP1_transformed[gene%in%Leafcutter_TSP7_transformed$gene,found_in_Leafcutter_TSP7:=1]
Leafcutter_TSP1_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
Leafcutter_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP7+found_in_SPLASH_plus_TSP1]
setorder(Leafcutter_TSP1_transformed, nrow, -sum_found, gene)
Leafcutter_TSP1_transformed = Leafcutter_TSP1_transformed[!duplicated(nrow)]

Leafcutter_TSP7_transformed[,found_in_Leafcutter_TSP1:=0]
Leafcutter_TSP7_transformed[,found_in_SPLASH_plus_TSP7:=0]
Leafcutter_TSP7_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
Leafcutter_TSP7_transformed[gene%in%SPLASH_plus_TSP7_transformed$gene,found_in_SPLASH_plus_TSP7:=1]
Leafcutter_TSP7_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_SPLASH_plus_TSP7]
setorder(Leafcutter_TSP7_transformed, nrow, -sum_found, gene)
Leafcutter_TSP7_transformed = Leafcutter_TSP7_transformed[!duplicated(nrow)]

SPLASH_plus_TSP1_transformed[,found_in_Leafcutter_TSP1:=0]
SPLASH_plus_TSP1_transformed[,found_in_SPLASH_plus_TSP7:=0]
SPLASH_plus_TSP1_transformed[gene%in%Leafcutter_TSP1_transformed$gene,found_in_Leafcutter_TSP1:=1]
SPLASH_plus_TSP1_transformed[gene%in%SPLASH_plus_TSP7_transformed$gene,found_in_SPLASH_plus_TSP7:=1]
SPLASH_plus_TSP1_transformed[,sum_found:=found_in_Leafcutter_TSP1+found_in_SPLASH_plus_TSP7]
setorder(SPLASH_plus_TSP1_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP1_transformed = SPLASH_plus_TSP1_transformed[!duplicated(nrow)]

SPLASH_plus_TSP7_transformed[,found_in_Leafcutter_TSP7:=0]
SPLASH_plus_TSP7_transformed[,found_in_SPLASH_plus_TSP1:=0]
SPLASH_plus_TSP7_transformed[gene%in%Leafcutter_TSP7_transformed$gene,found_in_Leafcutter_TSP7:=1]
SPLASH_plus_TSP7_transformed[gene%in%SPLASH_plus_TSP1_transformed$gene,found_in_SPLASH_plus_TSP1:=1]
SPLASH_plus_TSP7_transformed[,sum_found:=found_in_Leafcutter_TSP7+found_in_SPLASH_plus_TSP1]
setorder(SPLASH_plus_TSP7_transformed, nrow, -sum_found, gene)
SPLASH_plus_TSP7_transformed = SPLASH_plus_TSP7_transformed[!duplicated(nrow)]
blood_splicing_genes_all_methods = list(Leafcutter_TSP1 =  unique(Leafcutter_TSP1_transformed$gene), Leafcutter_TSP7 = unique(Leafcutter_TSP7_transformed$gene), SPLASH_plus_TSP1 = unique(SPLASH_plus_TSP1_transformed$gene), SPLASH_plus_TSP7 = unique(SPLASH_plus_TSP7_transformed$gene), SpliZ_TSP1 = spliz_blood[!gene%like%"unknown" & donor=="TSP1"]$gene, SpliZ_TSP7 = spliz_blood[!gene%like%"unknown" & donor=="TSP7"]$gene)



########################################################################
############# PLOT THE CONCORDANCE HEATMAP IN FIGURE 3B  ###############
########################################################################
method = c("Leafcutter","SpliZ","SPLASH+")
Intersection = c("(Lung) Donors 1,2")
value = c(length(intersect(lung_splicing_genes_all_methods$Leafcutter_TSP1,lung_splicing_genes_all_methods$Leafcutter_TSP2))/length(unique(c(lung_splicing_genes_all_methods$Leafcutter_TSP1,lung_splicing_genes_all_methods$Leafcutter_TSP2))), length(intersect(lung_splicing_genes_all_methods$SpliZ_TSP1,lung_splicing_genes_all_methods$SpliZ_TSP2))/length(unique(c(lung_splicing_genes_all_methods$SpliZ_TSP1,lung_splicing_genes_all_methods$SpliZ_TSP2))), length(intersect(lung_splicing_genes_all_methods$SPLASH_plus_TSP1,lung_splicing_genes_all_methods$SPLASH_plus_TSP2))/length(unique(c(lung_splicing_genes_all_methods$SPLASH_plus_TSP1,lung_splicing_genes_all_methods$SPLASH_plus_TSP2))))
dt_lung = data.table(method,Intersection,value)

method = c("Leafcutter","SpliZ","SPLASH+")
Intersection = c("(Blood) Donors 1,7")
value = c(length(intersect(blood_splicing_genes_all_methods$Leafcutter_TSP1,blood_splicing_genes_all_methods$Leafcutter_TSP7))/length(unique(c(blood_splicing_genes_all_methods$Leafcutter_TSP1,blood_splicing_genes_all_methods$Leafcutter_TSP7))), length(intersect(blood_splicing_genes_all_methods$SpliZ_TSP1,blood_splicing_genes_all_methods$SpliZ_TSP7))/length(unique(c(blood_splicing_genes_all_methods$SpliZ_TSP1,blood_splicing_genes_all_methods$SpliZ_TSP7))), length(intersect(blood_splicing_genes_all_methods$SPLASH_plus_TSP1,blood_splicing_genes_all_methods$SPLASH_plus_TSP7))/length(unique(c(blood_splicing_genes_all_methods$SPLASH_plus_TSP1,blood_splicing_genes_all_methods$SPLASH_plus_TSP7))))
dt_blood = data.table(method,Intersection,value)

method = c("Leafcutter","SpliZ","SPLASH+")
Intersection = c("(Muscle) Donors 1,2","(Muscle) Donors 1,2","(Muscle) Donors 1,2","(Muscle) Donors 1,4","(Muscle) Donors 1,4","(Muscle) Donors 1,4","(Muscle) Donors 2,4","(Muscle) Donors 2,4","(Muscle) Donors 2,4","(Muscle) Donors 1,2,4","(Muscle) Donors 1,2,4","(Muscle) Donors 1,2,4")
value = c(length(intersect(muscle_splicing_genes_all_methods$Leafcutter_TSP1,muscle_splicing_genes_all_methods$Leafcutter_TSP2))/length(unique(c(muscle_splicing_genes_all_methods$Leafcutter_TSP1,muscle_splicing_genes_all_methods$Leafcutter_TSP2))), length(intersect(muscle_splicing_genes_all_methods$SpliZ_TSP1,muscle_splicing_genes_all_methods$SpliZ_TSP2))/length(unique(c(muscle_splicing_genes_all_methods$SpliZ_TSP1,muscle_splicing_genes_all_methods$SpliZ_TSP2))), length(intersect(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,muscle_splicing_genes_all_methods$SPLASH_plus_TSP2))/length(unique(c(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,muscle_splicing_genes_all_methods$SPLASH_plus_TSP2))), length(intersect(muscle_splicing_genes_all_methods$Leafcutter_TSP1,muscle_splicing_genes_all_methods$Leafcutter_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$Leafcutter_TSP1,muscle_splicing_genes_all_methods$Leafcutter_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SpliZ_TSP1,muscle_splicing_genes_all_methods$SpliZ_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$SpliZ_TSP1,muscle_splicing_genes_all_methods$SpliZ_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4))), length(intersect(muscle_splicing_genes_all_methods$Leafcutter_TSP2,muscle_splicing_genes_all_methods$Leafcutter_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$Leafcutter_TSP2,muscle_splicing_genes_all_methods$Leafcutter_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SpliZ_TSP2,muscle_splicing_genes_all_methods$SpliZ_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$SpliZ_TSP2,muscle_splicing_genes_all_methods$SpliZ_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SPLASH_plus_TSP2,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4))/length(unique(c(muscle_splicing_genes_all_methods$SPLASH_plus_TSP2,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4))),length(intersect(muscle_splicing_genes_all_methods$Leafcutter_TSP1, intersect(muscle_splicing_genes_all_methods$Leafcutter_TSP2,muscle_splicing_genes_all_methods$Leafcutter_TSP4)))/length(unique(c(muscle_splicing_genes_all_methods$Leafcutter_TSP1,muscle_splicing_genes_all_methods$Leafcutter_TSP2,muscle_splicing_genes_all_methods$Leafcutter_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SpliZ_TSP1, intersect(muscle_splicing_genes_all_methods$SpliZ_TSP2,muscle_splicing_genes_all_methods$SpliZ_TSP4)))/length(unique(c(muscle_splicing_genes_all_methods$SpliZ_TSP1,muscle_splicing_genes_all_methods$SpliZ_TSP2,muscle_splicing_genes_all_methods$SpliZ_TSP4))), length(intersect(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,intersect(muscle_splicing_genes_all_methods$SPLASH_plus_TSP2,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4)))/length(unique(c(muscle_splicing_genes_all_methods$SPLASH_plus_TSP1,muscle_splicing_genes_all_methods$SPLASH_plus_TSP2,muscle_splicing_genes_all_methods$SPLASH_plus_TSP4)))          )
dt_muscle = data.table(method,Intersection,value)
dt=rbind(dt_lung,dt_blood,dt_muscle)
dt$method=factor(dt$method,levels = c("SPLASH+","Leafcutter","SpliZ"))
ggplot(dt, aes(Intersection,method)) + geom_tile(aes(fill = value))+scale_fill_gradient(low = "white", high = "blue") 
########################################################################################################################
########################################################################################################################