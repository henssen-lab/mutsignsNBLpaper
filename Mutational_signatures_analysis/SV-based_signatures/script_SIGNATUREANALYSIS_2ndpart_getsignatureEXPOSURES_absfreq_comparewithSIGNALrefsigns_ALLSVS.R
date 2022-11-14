##################################################################
###                 SV SIGNATURE ANALYSIS                      ###
### 2ND SCRIPT: GET EXPOSURES - COMPARE DE NOVO WITH REF SIGNS ###
##################################################################

#libraries
library(Palimpsest)
library(bedr)
library(RCircos)
library(stringr)

#declare directories
datadir <- "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis"
resdir <- file.path(datadir,paste("SV_signatures_ALLpatients_2022.06.15/",sep=""))


###### OPEN INPUT FILES ######

#load data
SV_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_intraSVs+translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt",sep="\t")

#open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
patients_type_file <- patients_type_file[,c(1,3)]
names(patients_type_file) <- c("sample_id", "risk_group")
patients_type_file <- patients_type_file[grep(patients_type_file$sample_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54",invert=TRUE),]

SV_file_merge <- merge(SV_file,patients_type_file,by.x="V1", by.y="sample_id")
names(SV_file_merge) <- c("sample_id", "chr1", "pos1", "chr2", "pos2", "info", "svtype", "strand", "length", "risk_group")

#exclude patients
SV_file_merge <- SV_file_merge[grep(SV_file_merge$sample_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54",invert=TRUE),]
SV_data <- SV_file_merge

#preprocess SV Data
SV_data <- SV_data[,c(1,7,2,3,4,5)]


###### FORMAT INPUT FILES ######

#format SV data
for (i in 1:dim(SV_data)[1]){
  if (SV_data[i,2]=="deletion"){
    SV_data[i,2] <- "DEL"
  } else if (SV_data[i,2]=="insertion"){
    SV_data[i,2] <- "INS"
  } else if (SV_data[i,2]=="translocation"){
    SV_data[i,2] <- "BND"
  } else if (SV_data[i,2]=="inversion"){
    SV_data[i,2] <- "INV"
  } else if (SV_data[i,2]=="tandem-duplication"){
    SV_data[i,2] <- "DUP"
  }
  SV_data[i,3] <- paste("chr",SV_data[i,3],sep="")
  SV_data[i,5] <- paste("chr",SV_data[i,5],sep="")
  SV_data[i,7] <- "."
  SV_data[i,8] <- "."
  SV_data[i,9] <- "."
}

#Get rid of all insertions
SV_data <- SV_data[SV_data$svtype!="INS",] #SV signatures do not take into account insertions

#New header for palimpsest
names(SV_data) <- c("Sample", "Type", "CHROM_1", "POS_1", "CHROM_2", "POS_2", "Tumor_Varcount", "Tumor_Depth", "Normal_Depth")


###### ANALYSIS ######

# Preprocess SV inputs and annotate for further analysis:
SV.vcf <- preprocessInput_sv(input_data =  SV_data, resdir = resdir)
SV_input <- palimpsest_input(vcf = SV.vcf,Type = "SV")

#Load from previous step
load("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/SV_denovo_signatures.RData")
SV_denovo_sigs <- denovo_signatures

# define list of colors for visualizing mutational signatures. Selecting default colors
SV_cols <- signature_colour_generator(rownames(SV_denovo_sigs))

# Calculate contribution of signatures in each sample:
SVsignatures_exp <- deconvolution_fit_SV(vcf = SV.vcf,input_data = SV_input$mut_nums,input_signatures = SV_denovo_sigs,sig_cols = SV_cols,resdir = resdir, gv="hg19")

#Filter signatures with exposures <5% of the total exposures in our cohort
def_signs = NULL

for (s in 1:dim(SVsignatures_exp$sig_num)[2]){
  if ((sum(SVsignatures_exp$sig_num[,s])*100/sum(SVsignatures_exp$sig_num))>=5){
    def_signs <- c(def_signs, colnames(SVsignatures_exp$sig_num)[s])
  }
}

SVsignatures_exp_FILT <- SVsignatures_exp$sig_num[names(SVsignatures_exp$sig_num) %in% def_signs]


###### SV REF SIGNATURE COMPARISON ########

# select reference signatures
SV_ref_signs <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/compare_SVsign_denovo_vs_REF/REF_SVsign_degasperi_et_al_v2.txt", sep="\t", stringsAsFactors=FALSE, header=TRUE)

#match columns between both types of signatures
SV_denovo_sigs_order <- SV_denovo_sigs[,c(8,12,14,10,4,20,24,26,22,16,32,36,38,34,28,2,7,11,13,9,3,19,23,25,21,15,31,35,37,33,27,1)]
colnames(SV_denovo_sigs_order) <- names(SV_ref_signs)

# Compare the de novo signatures with published SV signatures
SV_cosine_match <- compare_results(reference_sigs = SV_ref_signs, extraction_1 = SV_denovo_sigs_order)
#write results
write.table(SV_cosine_match, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/denovosigns_vs_SIGNALrefsigns_cosinedistance.txt", sep="\t", row.names=FALSE, quote=FALSE)
SV_cosine_similarities <- deconvolution_compare(SV_denovo_sigs_order, SV_ref_signs) #cosine matrix for plotting
#write results
write.table(SV_cosine_match, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/denovosigns_vs_SIGNALrefsigns_cosinedistance_MATRIXforplots.txt", sep="\t", row.names=FALSE, quote=FALSE)


###### GET SIGNATURE EXPOSURE ########

#COUNTS
SVsignatures_exp_abs <- data.frame(SVsignatures_exp$sig_num)
SVsignatures_exp_abs$sample_id <- rownames(SVsignatures_exp_abs)
abs_exp <- merge(SVsignatures_exp_abs, patients_type_file, by="sample_id", all=TRUE)
abs_exp[is.na(abs_exp)] <- 0
write.table(abs_exp, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/exposure_SV_denovo+SIGNAL_signature_absexp_per_patient_COUNT_22.06.13.txt",sep="\t", row.names=FALSE, quote=FALSE)
#FREQ
SVsignatures_exp_freq <- data.frame(SVsignatures_exp$sig_props)
SVsignatures_exp_freq$sample_id <- rownames(SVsignatures_exp_freq)
freq_exp <- merge(SVsignatures_exp_freq, patients_type_file, by="sample_id", all=TRUE)
freq_exp[is.na(freq_exp)] <- 0
write.table(freq_exp, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/exposure_SV_denovo+SIGNAL_signature_absexp_per_patient_FREQ_22.06.13.txt",sep="\t", row.names=FALSE, quote=FALSE)

# Estimate the probability of each SV being due to each mutational signature - in case we want to use this info
SV_origins <- signature_origins(input = SV.vcf,Type = "SV",signature_contribution = SVsignatures_exp,input_signatures = SV_denovo_sigs)
write.table(data.frame(SV_origins), "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/table_SVmutsignature_origin_foreachSV_info_22.06.13.txt",sep="\t", row.names=FALSE, quote=FALSE)
