##############################################
###          SV SIGNATURE ANALYSIS         ###
### 1ST SCRIPT: DE NOVO SIGNATURE ANALYSIS ###
##############################################

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
SV_denovo_sigs <- NMF_Extraction(input_matrices =  SV_input, range_of_sigs = 1:10, nrun = 500, resdir = resdir)

