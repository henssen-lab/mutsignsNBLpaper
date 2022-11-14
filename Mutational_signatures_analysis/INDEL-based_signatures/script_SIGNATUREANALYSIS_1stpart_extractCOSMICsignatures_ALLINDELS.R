#############################################
###        INDEL SIGNATURE ANALYSIS       ###
### 1ST SCRIPT: EXTRACT COSMIC SIGNATURES ###
#############################################


#libraries
library("YAPSA")

#data load
data(sigs_pcawg)
data(cutoffs_pcawg)


###### FORMAT FILES ######

#open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
patients_type_file <- patients_type_file[,c(1,3)]
names(patients_type_file) <- c("sample_id", "risk_group")

#open indels list and format
indels_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_21.10.13_ALLINDELS_perpatient_perriskgroup/RESULTS_all_patients_berlin+peifer_cohorts_INDELS_mutect2_PASS_newinfostrands_21.01.26_simplified.txt", sep="\t")
names(indels_file) <- c("sample_id", "chrom", "pos", "ref", "alt")
indels_file_PT <- merge(indels_file, patients_type_file, by="sample_id")
names(indels_file_PT) <- c("sample_id", "chrom", "pos", "ref", "alt", "risk_group")

#exclude patients
indels_file_PT <- indels_file_PT[grep(indels_file_PT$sample_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54",invert=TRUE),]

#final indels list
indels_file_PT_final <- indels_file_PT[,c(2:5,1)]
names(indels_file_PT_final) <- c("CHROM", "POS", "REF", "ALT", "PID")

#format COSMIC signatures
cosmic_signs <- read.delim("/data/gpfs-1/users/rodrigue_c/work/refs/COSMIC_v3.2_INDELsignatures_ID_GRCh37.txt",sep="\t")
rownames(cosmic_signs) <- rownames(PCAWG_SP_ID_sigs_df)
cosmic_signs <- cosmic_signs[,-1]


###### ANALYSIS ######

#data preprocessing
vcf_like_indel_trans_df <- translate_to_hg19(indels_file_PT_final,"CHROM") #add "chr" to the chromosome field
mutational_cataloge_indel_df <- create_indel_mutation_catalogue_from_df(in_dat = vcf_like_indel_trans_df, in_signature_df = cosmic_signs) #generate the mutation counts for each feature
save(mutational_cataloge_indel_df, file="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/ALLINDELS_COSMIC3.2_signatures_mutcatalogforanalysis_22.06.13.RData")

#prepare parameters for the analysis
#cutoff df
current_cutoff_pid_vector <- cutoffPCAWG_ID_WGS_Pid_df[3,] #For Indel signatures the optimal cost factor was identified to be 3, thus the third line in cutoffPCAWG_ID_WGS_Pid_df can be chosen for analysis.
current_cutoff_pid_vector <- unlist(append(current_cutoff_pid_vector, 0, 14)) #we don't know the PCAWG cutoff for ID15 and ID18 - we include them with cutoff=0
current_cutoff_pid_vector <- unlist(append(current_cutoff_pid_vector, 0, 17))
names(current_cutoff_pid_vector) <- names(cosmic_signs)
#aetiology df
current_sigInd_df <- PCAWG_SP_ID_sigInd_df[1:14,]
current_sigInd_df$sig <- as.character(current_sigInd_df$sig)
current_sigInd_df[15,]<-c("ID15",15,"grey","unknown") #same, ID15 and ID18 are not included in their table, so we include them with the aetiology from COSMIC.
current_sigInd_df[16,]<-c("ID16",16,"pink","unknown")
current_sigInd_df[17,]<-c("ID17",17,"black","TOP2A")
current_sigInd_df[18,]<-c("ID18",18,"mediumblue","Colibactin exposure")


### SUPERVISED MUTATIONAL SIGNATURE ANALYSIS
current_LCDlistsList <- LCD_complex_cutoff_combined(in_mutation_catalogue_df = mutational_cataloge_indel_df, in_signatures_df = cosmic_signs, in_cutoff_vector = current_cutoff_pid_vector, in_filename = NULL, in_method = "abs", in_sig_ind_df = current_sigInd_df)

#Get the absolute exposures for all patients in our cohort
cohort_exp <- t(current_LCDlistsList$cohort$exposures)

#Filter signatures with exposures <5% of the total exposures in our cohort
def_signs = NULL

for (s in 1:dim(cohort_exp)[2]){
  if ((sum(cohort_exp[,s])*100/sum(cohort_exp))>=5){
    def_signs <- c(def_signs, colnames(cohort_exp)[s])
  }
}


### GET EXPOSURE FOR THE FINAL SIGNATURES PASSING THE FILTERS

#Re-generate the params for the analysis (only for filtered signatures)
cosmic_signs_FILT <- cosmic_signs[names(cosmic_signs) %in% def_signs]
current_cutoff_pid_vector_FILT <- current_cutoff_pid_vector[names(current_cutoff_pid_vector) %in% def_signs]
current_sigInd_df_FILT <- current_sigInd_df[current_sigInd_df$sig %in% def_signs,]

#Analysis
current_LCDlistsList_FILT <- LCD_complex_cutoff_combined(in_mutation_catalogue_df = mutational_cataloge_indel_df, in_signatures_df = cosmic_signs_FILT, in_cutoff_vector = current_cutoff_pid_vector_FILT, in_filename = NULL, in_method = "abs", in_sig_ind_df = current_sigInd_df_FILT)

# Retrieve absolute exposures (results)
current_LCDlistsList_FILT_abs <- data.frame(t(current_LCDlistsList_FILT$cohort$exposures))
current_LCDlistsList_FILT_abs$sample_id <- rownames(current_LCDlistsList_FILT_abs)
abs_exp <- merge(current_LCDlistsList_FILT_abs, patients_type_file, by="sample_id")
#write results
write.table(abs_exp, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_COUNT_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

# Retrieve relative exposures (results)
current_LCDlistsList_FILT_freq <- data.frame(t(current_LCDlistsList_FILT$cohort$norm_exposures))
current_LCDlistsList_FILT_freq$sample_id <- rownames(current_LCDlistsList_FILT_freq)
abs_freq <- merge(current_LCDlistsList_FILT_freq, patients_type_file, by="sample_id")
#write results
write.table(abs_freq, "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

