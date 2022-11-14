################################################################
#### SCRIPT CORRELATE ALL SIGNATURES BETWEEN THEM (PEARSON) ####
################################################################


#libraries
library("Hmisc")
library("corrplot")
library("stringr")

## Prepare data frame
DT_ALL_PATIENTS = NULL


####### OPEN INPUT FILES #######

#SNV signatures
SNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/exposure_SNV_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt")
SNVsign <- SNVsign[,-dim(SNVsign)[2]]
SNVsign_names <- str_sort(names(SNVsign[-1]),numeric=TRUE)

#INDELS signatures
INDsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt")
INDsign <- INDsign[,-dim(INDsign)[2]]
INDsign_names <- str_sort(names(INDsign[-1]),numeric=TRUE)

#SV signatures
SVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/exposure_SV_denovo+SIGNAL_signature_absexp_per_patient_FREQ_22.06.13.txt")
SVsign <- SVsign[,-dim(SVsign)[2]]
SVsign_names <- str_sort(names(SVsign[-1]),numeric=TRUE)

#CNV signatures
CNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/cnv_signatures_analysis/cnv_signature_analysis_NEWfilters_22.06.15/exposure_CNA_denovo_signature_freqexp_per_patient_FREQ_22.06.13_REFITFREQ.txt")
CNVsign <- CNVsign[,-dim(CNVsign)[2]]
CNVsign_names <- str_sort(names(CNVsign[-1]),numeric=TRUE)


####### CORRELATION TEST #######

#Final table
DT_ALL_PATIENTS2 <- merge(SNVsign, INDsign, by="sample_id")
DT_ALL_PATIENTS1 <- merge(DT_ALL_PATIENTS2, SVsign, by="sample_id", all = TRUE)
DT_ALL_PATIENTS <- merge(DT_ALL_PATIENTS1, CNVsign, by="sample_id", all = TRUE)

#convert format sample id to rownames
rownames(DT_ALL_PATIENTS) <- DT_ALL_PATIENTS[,1]
DT_ALL_PATIENTS <- DT_ALL_PATIENTS[,-1]

#order signatures
NAMES_all_signs <- c(SNVsign_names, INDsign_names, CNVsign_names, SVsign_names)
DT_ALL_PATIENTS <- DT_ALL_PATIENTS[NAMES_all_signs]

#correlation matrix
CORMAT_data <- cor(DT_ALL_PATIENTS,method="pearson",use="na.or.complete")
signatures <- rownames(CORMAT_data)
save_CORMAT_data <- cbind(signatures, data.frame(CORMAT_data))
write.table(save_CORMAT_data, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sign_coocurrence/plots_exposurebased_NEWresults_postREVISION_22.06.17/correlationmatrix_pearsoncorr_basedonexposure_newFREQexp_ONLYsignatures_SNV+INDEL+SV+CNV_ALLPATIENTS_postREVISION_22.06.17.txt", sep="\t", quote=FALSE, row.names=FALSE)

#test pvalue correlation
CORMAT_data_test <- cor.mtest(DT_ALL_PATIENTS, conf.level = 0.95, method = "pearson", exact=TRUE)


####### PLOT RESULTS CORRELATION TEST #######

#PLOT
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/sign_coocurrence/plots_exposurebased_NEWresults_postREVISION_22.06.17/plots/PLOT_correlationmatrix_basedonexposuredata_ONLYsignatures_SNV+INDEL+SV+CNVnewsignIDS_ALLPATIENTS_postREVISION_22.06.17.pdf", width = 10, height = 10)
corrplot(CORMAT_data, p.mat = CORMAT_data_test$p, order = 'original', insig = 'blank', sig.level = 0.05, pch.cex = 0.9, pch.col = 'grey20', method="color", type="lower", col=col1(30))
dev.off()
