##############################################################################################
#### SCRIPT CORRELATE ALL SIGNATURES WITH COMPLEX REARRANGEMENTS AND GENES THEM (PEARSON) ####
##############################################################################################


#libraries
library("Hmisc")
library("corrplot")
library("stringr")
library("ClustImpute")

## Prepare data frame
DT_ALL_PATIENTS = NULL


#### OPEN INPUT FILES ####

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
CNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis//cnv_signatures_analysis/cnv_signature_analysis_NEWfilters_22.06.15/exposure_CNA_denovo_signature_freqexp_per_patient_FREQ_22.06.13_REFITFREQ.txt")
CNVsign <- CNVsign[,-dim(CNVsign)[2]]
CNVsign_names <- str_sort(names(CNVsign[-1]),numeric=TRUE)

#Genes NBL mutated
NBLgenes <- read.delim("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/inputs/ALLpatients_allmutations_NBLgenesessential_v2.txt")
NBLgenes <- NBLgenes[,c(dim(NBLgenes)[2],1:(dim(NBLgenes)[2]-1))]
NBLgenes_names <- names(NBLgenes[-1])

#Genes DNA repair
DNAgenes <- read.delim("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/inputs/ALLpatients_SOMATICmutations_DNAgenesrepair+PROfound_genesHRR.txt")
DNAgenes <- DNAgenes[,c(dim(DNAgenes)[2],1:(dim(DNAgenes)[2]-1))]
DNAgenes_names <- names(DNAgenes[-1])

#Complex rearrangements
comprearrang <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_FILTEREDATLEAST2CALLERS+AAnormal/MATRIX_complexrearrangements_types_COUNTS_perpatient_114patients_RISKGROUPinfo_FILTEREDATLEAST2CALLERS+AAnormal.txt", sep="\t", stringsAsFactors=FALSE)
comprearrang <- comprearrang[,-dim(comprearrang)[2]]
comprearrang_names <- names(comprearrang[-1])

#### CREATE MATRIX ####

#Transform NA counts into 0 for no presence of mutations in genes
genes_T1 <- merge(SNVsign, NBLgenes, by="sample_id",all = TRUE)
genes_T2 <- merge(genes_T1, DNAgenes, by="sample_id",all = TRUE)
genes_T2 <- genes_T2[,-c(2:5)]
genes_T2[is.na(genes_T2)] <- 0
genes_T2[,-1] <- scale(genes_T2[,-1])

#Final table
DT_ALL_PATIENTS1 <- merge(SNVsign, INDsign, by="sample_id", all = TRUE)
DT_ALL_PATIENTS2 <- merge(DT_ALL_PATIENTS1, SVsign, by="sample_id", all = TRUE)
DT_ALL_PATIENTS3 <- merge(DT_ALL_PATIENTS2, CNVsign, by="sample_id", all = TRUE)
DT_ALL_PATIENTS4 <- merge(DT_ALL_PATIENTS3, genes_T2, by="sample_id", all = TRUE)
DT_ALL_PATIENTS <- merge(DT_ALL_PATIENTS4, comprearrang, by="sample_id", all = TRUE)

#convert format sample id to rownames
rownames(DT_ALL_PATIENTS) <- DT_ALL_PATIENTS[,1]
DT_ALL_PATIENTS <- DT_ALL_PATIENTS[,-1]

#DT_ALL_PATIENTS[is.na(DT_ALL_PATIENTS)] <- 0
df_cluster_input_N <- ClustImpute(DT_ALL_PATIENTS,nr_cluster=3, nr_iter=25, c_steps=25, n_end=10) #In case we want to impute the missing values
DT_ALL_PATIENTS <- df_cluster_input_N$complete_data

#order signatures
NAMES_all_signs <- c(SNVsign_names, INDsign_names, CNVsign_names, SVsign_names, NBLgenes_names, DNAgenes_names, comprearrang_names)
DT_ALL_PATIENTS <- DT_ALL_PATIENTS[NAMES_all_signs]
DT_ALL_PATIENTS <- scale(DT_ALL_PATIENTS)

#correlation matrix
CORMAT_data <- cor(DT_ALL_PATIENTS,method="spearman",use="na.or.complete")
signatures <- rownames(CORMAT_data)
save_CORMAT_data <- cbind(signatures, data.frame(CORMAT_data))
write.table(save_CORMAT_data, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/analysis_correlation_postREVISION_22.07.04/correlationmatrix_spearmancorr_basedonexposure_newFREQexp_signatures+comprearrang+genes+PROfound_genesHRR_ALLPATIENTS_postREVISION_22.11.2.txt", sep="\t", quote=FALSE, row.names=FALSE)

#test pvalue correlation
CORMAT_data_test <- cor.mtest(DT_ALL_PATIENTS, conf.level = 0.95, method = "spearman", exact=FALSE)


###PLOT SIGN + COMPLEX REARRANGEMENTS VS GENES

a <- CORMAT_data[c(24:(length(NBLgenes_names)+length(DNAgenes_names)+23)),c(1:23,(length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]
b <- CORMAT_data_test$p[c(24:(length(NBLgenes_names)+length(DNAgenes_names)+23)),c(1:23,(length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]

#get rid of empty rows
coeff = NULL
for (i in 1:dim(b)[1]){
    if (min(b[i,])>=0.05){
        coeff <- c(coeff, i)
    }
}

a_1 <- a[-coeff,]
b_1 <- b[-coeff,]

#PLOT
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/analysis_correlation_postREVISION_22.07.04/plots/PLOT_correlationmatrix_basedonexposuredata_signatures+genes+complexrearrang_ALLPATIENTS_PASS_NEWEXPOSURE_postREVISION_22.11.2_PROfound_genesHRR.pdf", width = 7, height = 10)
corrplot(a_1, p.mat = b_1, insig = 'blank',sig.level = 0.05, pch.cex = 0.9, pch.col = 'grey20', method="color", col=col1(30),diag = TRUE)
dev.off()


###PLOT COMPLEX REARRANGEMENTS VS COMPLEX REARRANGEMENTS

a <- CORMAT_data[c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data))),c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]
b <- CORMAT_data_test$p[c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data))),c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]

#PLOT
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/analysis_correlation_postREVISION_22.07.04/plots/PLOT_correlationmatrix_basedonexposuredata_complexrearrang_ALLPATIENTS_PASS_NEWEXPOSURE_postREVISION_22.09.27.pdf", width = 7, height = 10)
corrplot(a, p.mat = b, insig = 'blank', sig.level = 0.05, pch.cex = 0.9, pch.col = 'grey20', method="color", type="lower", col=col1(30))
dev.off()


###PLOT SIGN VS COMPLEX REARRANGEMENTS

a <- CORMAT_data[c(1:23),c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]
b <- CORMAT_data_test$p[c(1:23),c((length(NBLgenes_names)+length(DNAgenes_names)+23+1):length(rownames(CORMAT_data)))]

#PLOT
col1 <- colorRampPalette(c("#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#FFFFFF", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"))

pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/comprearrang_genes_signs_coocurrence/analysis_correlation_postREVISION_22.07.04/plots/PLOT_correlationmatrix_basedonexposuredata_signatures+complexrearrang_ALLPATIENTS_PASS_NEWEXPOSURE_postREVISION_22.09.27.pdf", width = 7, height = 10)
corrplot(a, p.mat = b, insig = 'blank', sig.level = 0.05, pch.cex = 0.9, pch.col = 'grey20', method="color", col=col1(30),diag = TRUE)
dev.off()
