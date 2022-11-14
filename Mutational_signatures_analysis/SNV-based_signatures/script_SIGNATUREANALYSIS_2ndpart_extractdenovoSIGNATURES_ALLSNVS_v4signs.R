##############################################
###           SIGNATURE ANALYSIS           ###
### 2ND SCRIPT: EXTRACT DE NOVO SIGNATURES ###
##############################################


#libraries
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
patients_type_file <- patients_type_file[,c(1:3)]

###load previous object x_snvs_muttype_counts
load("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_mutation_type_counts_TRInucleotidecontext_22.06.13.RData")


###### SIGNATURE ANALYSIS ######

### DE-NOVO Signature extraction
num.sign <- 4 #how many signatures should we extract
# Define parameters for the non-negative matrix factorization procedure (recommended 500-1000 iterations)
Cancer.params <- mutSignatures::setMutClusterParams(num.sign, num_totIterations = 500, num_parallelCores = 16, algorithm= "brunet", approach = "counts")
Cancer.analysis <- decipherMutationalProcesses(input = x_snvs_muttype_counts, params = Cancer.params)
save(Cancer.analysis, file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_denovo_mutSignatures_analysis_results_absexp_4signs_22.06.13.RData")
# Retrieve signatures (results)
deNovo.signs <- Cancer.analysis$Results$signatures
# Retrieve exposures (results)
deNovo.exp <- Cancer.analysis$Results$exposures
deNovo.exp_DF <- data.frame(t(coerceObj(x = deNovo.exp, to = "data.frame")))
deNovo.exp_DF$sample_id <- rownames(deNovo.exp_DF)
#Retrieve risk group
merge_deNovo.exp_DF <- merge(deNovo.exp_DF, patients_type_file, by="sample_id")
write.table(merge_deNovo.exp_DF,"/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_table_absEXPOSURE_denovo_4signs_all114patients_500perms_22.06.13.txt",row.names=FALSE,sep="\t",quote=FALSE)

