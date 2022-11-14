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
num.sign <- 7 #how many signatures should we extract
pdf(file=paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/plot_silhouette_500perms_k",num.sign,".pdf",sep=""))
# Define parameters for the non-negative matrix factorization procedure (recommended 500-1000 iterations)
Cancer.params <- mutSignatures::setMutClusterParams(num.sign, num_totIterations = 500, num_parallelCores = 32, algorithm= "brunet", approach = "counts")
Cancer.analysis <- decipherMutationalProcesses(input = x_snvs_muttype_counts, params = Cancer.params)
dev.off()

