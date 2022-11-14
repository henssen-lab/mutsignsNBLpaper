##############################################################
###                   SIGNATURE ANALYSIS                   ###
### 3RD SCRIPT: EXTRACT AND COMPARE WITH COSMIC SIGNATURES ###
##############################################################


#libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(wesanderson)


###### OPEN INPUT FILES ######

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
names(patients_type_file) <- c("sample_id", "cohort", "risk_group")
patients_type_file <- patients_type_file[,c(1,3)]

###load previous object x_snvs_muttype_counts and Cancer.analysis
load("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_mutation_type_counts_TRInucleotidecontext_22.06.13.RData")
load("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_denovo_mutSignatures_analysis_results_absexp_3signs_22.06.13.RData")


###### SIGNATURE ANALYSIS ######

### DE-NOVO Signature extraction
# Retrieve signatures (results)
deNovo.signs <- Cancer.analysis$Results$signatures
# Retrieve exposures (results)
deNovo.exp <- Cancer.analysis$Results$exposures
deNovo.exp_DF <- data.frame(t(coerceObj(x = deNovo.exp, to = "data.frame")))
deNovo.exp_DF$sample_id <- rownames(deNovo.exp_DF)
#Retrieve risk group
merge_deNovo.exp_DF <- merge(deNovo.exp_DF, patients_type_file, by="sample_id")


###### COMPARE DENOVO SIGNATURES WITH KNOWN COSMIC ONES ######

# Get signatures from data (imported as data.frame)
cosmic_signs_DF <- read.delim("/data/gpfs-1/users/rodrigue_c/work/refs/COSMIC_v3.2_SBS_GRCh37.txt",sep="\t")
rownames(cosmic_signs_DF) <- cosmic_signs_DF$Type
#Convert it to mutSignatures object
cosmixSigs <- cosmic_signs_DF %>% dplyr::select(starts_with("SBS")) %>% as.mutation.signatures()
# Compare de-novo signatures with selected COSMIC signatures
match_signs <- matchSignatures(mutSign = deNovo.signs , reference = cosmixSigs, threshold = 0.5, plot = TRUE)
save(match_signs, file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_denovo_vs_COSMIC3.2_signature_comparison_3denovo+absexp_22.06.13.RData")

#Save table with distance values of de-novo signatures vs COSMIC signatures: TOP similar signatures (cosine distance<0.25)
distance_DF <- match_signs$distanceDataFrame[order(match_signs$distanceDataFrame[3]),]
TOP_similar_signs <- distance_DF[distance_DF$dist<0.15,]
write.table(TOP_similar_signs,"/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_table_comparison_3denovo+absexp_vs_COSMIC3.2_all114patients_500perms_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

#Plot cosine distances between de novo and COSMIC signatures
pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/plots/Plot_ALLSNVS_NBL_114patients_Peifer+Berlin_signaturecomparison_3denovo+absexp_vs_COSMIC3.2_22.06.13.pdf", width = 10, height = 4)
match_signs$plot
dev.off()


###### ESTIMATE EXPOSURES TO COSMIC SIGNATURES ######

### FILTER ALL SIGNATURES WITH <5% OF EXPOSURE IN OUR COHORT

#Retrieve the known COSMIC signatures I'm interested in.
#The algorithm is way more accurate if ONLY AND ALL the relevant signatures are used (i.e., the signatures we are reasonably expecting to be operative in the tumor samples being analyzed)
cosmixSigs_NBL <- cosmic_signs_DF %>% dplyr::select(ends_with(unique(TOP_similar_signs$refSign))) %>% as.mutation.signatures()
cosmic_signs_INDATA <- resolveMutSignatures(mutCountData = x_snvs_muttype_counts, signFreqData = cosmixSigs_NBL)

# Retrieve exposures (results)
cosmic_exp <- cosmic_signs_INDATA$results$count.result
cosmic_exp_DF_filter <- data.frame(t(coerceObj(x = cosmic_exp, to = "data.frame")))

#Retrieve signatures with exposures <5% of the total exposures in our cohort
def_signs = NULL

for (s in 1:dim(cosmic_exp_DF_filter)[2]){
  if ((sum(cosmic_exp_DF_filter[,s])*100/sum(cosmic_exp_DF_filter))>=5){
    def_signs <- c(def_signs, names(cosmic_exp_DF_filter)[s])
  }
}

### GET EXPOSURE FOR THE FINAL SIGNATURES PASSING THE FILTERS

cosmixSigs_NBL_FINAL <- cosmic_signs_DF %>% dplyr::select(ends_with(unique(def_signs))) %>% as.mutation.signatures()
cosmic_signs_INDATA_FINAL <- resolveMutSignatures(mutCountData = x_snvs_muttype_counts, signFreqData = cosmixSigs_NBL_FINAL)

# Retrieve exposures (results)
cosmic_exp_FINAL <- cosmic_signs_INDATA_FINAL$results$count.result
cosmic_exp_DF_FINAL <- data.frame(t(coerceObj(x = cosmic_exp_FINAL, to = "data.frame")))
cosmic_exp_DF_FINAL$sample_id <- rownames(cosmic_exp_DF_FINAL)
#Retrieve risk group
merge_cosmic_exp_DF_FINAL <- merge(cosmic_exp_DF_FINAL, patients_type_file, by="sample_id")
#write results
write.table(merge_cosmic_exp_DF_FINAL, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/exposure_SNV_COSMIC_signature_absexp_per_patient_COUNT_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

# Retrieve relative exposures per sample (results)
cosmic_exp_FINAL2 <- cosmic_signs_INDATA_FINAL$results$freq.result
cosmic_exp_DF_FINAL2 <- data.frame(t(coerceObj(x = cosmic_exp_FINAL2, to = "data.frame")))
cosmic_exp_DF_FINAL2$sample_id <- rownames(cosmic_exp_DF_FINAL2)
#Retrieve risk group
merge_cosmic_exp_DF_FINAL2 <- merge(cosmic_exp_DF_FINAL2, patients_type_file, by="sample_id")
#write results
write.table(merge_cosmic_exp_DF_FINAL2, "/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/exposure_SNV_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

