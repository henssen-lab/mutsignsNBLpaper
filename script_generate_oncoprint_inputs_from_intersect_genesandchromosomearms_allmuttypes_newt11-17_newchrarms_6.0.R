#### SCRIPT TO SELECT MUTATED GENES AND REGIONS OF INTEREST IN NBL
#### FROM THE INTERSECT FILES
#### GENERATE ONCOPRINT INPUTS
#### INCLUDING TRANSLOC T(11;17)

snv_FINAL = NULL
indels_FINAL = NULL
svs_FINAL = NULL
svs_close_FINAL = NULL
cnvs_FINAL = NULL

snv_FINAL_oncohist = NULL
indels_FINAL_oncohist = NULL
svs_FINAL_oncohist = NULL
svs_close_FINAL_oncohist = NULL

#Open genes list
gene_list_essential <- read.table("/fast/users/rodrigue_c/work/refs/KNOWNgenes_NBL_essentials_v2.0.txt")
#gene_list_essential <- read.table("/fast/users/rodrigue_c/work/refs/KNOWNgenes_NBL_ALL.txt")
gene_list_oncohist <- read.table("/fast/users/rodrigue_c/work/refs/KNOWNgenes_oncohistones.txt")
#gene_list_all <- read.table("/fast/users/rodrigue_c/work/refs/KNOWNgenes_NBL_all+histones.txt")

#Open intersect files
snv <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/vep_analysis/results_vep_format_filtered_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_SNVS_mutect_PASS_21.04.07.txt.bedformat.genecodev19.onlygenes.protein_coding.vepinfo")
indels <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/vep_analysis/results_vep_format_filtered_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_INDELS_mutect2_PASS_newinfostrands_21.04.07.txt.genecodev19.onlygenes.protein_coding.vepinfo")
#svs_bkps_entiresv <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/intersect_mutations_ALLgenes_gencodev19/results_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt.bedformat.genecodev19.onlygenes.protein_coding.bkps+entiresv")
svs_bkps <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/intersect_mutations_ALLgenes_gencodev19/results_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_intraSVs+translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt.bedformat.genecodev19.onlygenes.protein_coding")
svs_close_bkps <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/intersect_mutations_ALLgenes_gencodev19/results_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_AROUND20KBP_BKPS_intraSVs+translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.04.28.txt.bedformat.genecodev19.onlygenes.protein_coding")
cnvs <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/intersect_mutations_ALLgenes_gencodev19/results_21.04.07/RESULTS_all_patients_berlin+peifer_cohorts_CNV_GAINS+LOSSES+AMP+HOMDELS_ASCAT_martin_results_cutoffploidy+cutoff1.25_ALLpatients_21.04.21.txt.bedformat.genecodev19.onlygenes.protein_coding")
cnv_arms <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/intersect_mutations_ALLgenes_gencodev19/results_21.04.07_chromosomearms/RESULTS_all_patients_berlin+peifer_cohorts_CNV_GAINS+LOSSES+AMP+HOMDELS_ASCAT_martin_results_cutoffploidy+cutoff1.25_ALLpatients_21.04.21.txt.bedformat.genomicregionsNBL.overlap60")

#Grep mutations associated with genes essential
for (i in 1:dim(gene_list_essential)[1]){
  snv_grep <- unique(snv[snv$V4==gene_list_essential[i,1],c(4,5,7,12)])
  indels_grep <- unique(indels[indels$V4==gene_list_essential[i,1],c(4,5,7,12)])
  svs_grep <- unique(svs_bkps[svs_bkps$V5==gene_list_essential[i,1],c(5,6,8)])
  svs_close_grep <- unique(svs_close_bkps[svs_close_bkps$V5==gene_list_essential[i,1],c(5,6,8)])
  cnvs_grep <- unique(cnvs[cnvs$V5==gene_list_essential[i,1],c(5,6,8)])
  snv_FINAL <- rbind(snv_FINAL, snv_grep)
  indels_FINAL <- rbind(indels_FINAL, indels_grep)
  svs_FINAL <- rbind(svs_FINAL, svs_grep)
  svs_close_FINAL <- rbind(svs_close_FINAL, svs_close_grep)
  cnvs_FINAL <- rbind(cnvs_FINAL, cnvs_grep)
}

#Grep mutations associated with genes from oncohistones
for (i in 1:dim(gene_list_oncohist)[1]){
  snv_grep <- unique(snv[snv$V4==gene_list_oncohist[i,1],c(4,5,7,12)])
  indels_grep <- unique(indels[indels$V4==gene_list_oncohist[i,1],c(4,5,7,12)])
  svs_grep <- unique(svs_bkps[svs_bkps$V5==gene_list_oncohist[i,1],c(5,6,8)])
  svs_close_grep <- unique(svs_close_bkps[svs_close_bkps$V5==gene_list_oncohist[i,1],c(5,6,8)])
  snv_FINAL_oncohist <- rbind(snv_FINAL_oncohist, snv_grep)
  indels_FINAL_oncohist <- rbind(indels_FINAL_oncohist, indels_grep)
  svs_FINAL_oncohist <- rbind(svs_FINAL_oncohist, svs_grep)
  svs_close_FINAL_oncohist <- rbind(svs_close_FINAL_oncohist, svs_close_grep)
}

#Format cnv chromosome arms
cnv_arms$V11 <- paste(cnv_arms$V5,"_",cnv_arms$V8, sep="")
cnv_arms_FINAL <- cnv_arms[,c(11,6,8,9)]

#Grep transloc
sv_all <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_intraSVs+translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt")

#SV translocations t(11;17) and t(1;17)
sv_all <- sv_all[sv_all$V7=="translocation",]
t11_17_1 <- data.frame(patient_id=sv_all[sv_all$V2=="11"&sv_all$V3>=53700000&sv_all$V3<=135006516&sv_all$V4=="17"&sv_all$V5>=24000000&sv_all$V5<=81195210,1])
t11_17_2 <- data.frame(patient_id=sv_all[sv_all$V2=="17"&sv_all$V3>=24000000&sv_all$V3<=81195210&sv_all$V4=="11"&sv_all$V5>=53700000&sv_all$V5<=135006516,1])
t11_17 <- unique(rbind(t11_17_1,t11_17_2))

t11_17$mut_type <- "translocation"
t11_17$Gene.Symbol <- "t(11;17)"

t11_17_FINAL <- t11_17[,c(3,1,2)]
names(t11_17_FINAL) <- c("Gene.Symbol", "sample_id", "mut_type")

#Get everything under the same format c("Gene.Symbol", "sample_id", "mut_type", "conseq")
svs_FINAL$conseq <- "."
svs_close_FINAL$conseq <- "."
cnvs_FINAL$conseq <- "."
t11_17_FINAL$conseq <- "."
if (dim(svs_FINAL_oncohist)[1]>0){
  svs_FINAL_oncohist$conseq <- "."
}
if (dim(svs_close_FINAL_oncohist)[1]>0){
  svs_close_FINAL_oncohist$conseq <- "."
}

#Save the files for the oncoprint
write.table(snv_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SNVS_oncoprint_NBLgenesessential.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(indels_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/INDELS_oncoprint_NBLgenesessential.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(svs_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_oncoprint_NBLgenesessential.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(svs_close_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_around20kbp_oncoprint_NBLgenesessential.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(cnvs_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/CNVS_oncoprint_NBLgenesessential.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(cnv_arms_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/CNVS_oncoprint_NBLchromarms_v2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(snv_FINAL_oncohist, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SNVS_oncoprint_NBLoncohistones.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(indels_FINAL_oncohist, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/INDELS_oncoprint_NBLoncohistones.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(svs_FINAL_oncohist, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_oncoprint_NBLoncohistones.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(svs_close_FINAL_oncohist, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_around20kbp_oncoprint_NBLoncohistones.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(t11_17_FINAL, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/Transloc11_17_oncoprint_NBLchromarms_v2.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

