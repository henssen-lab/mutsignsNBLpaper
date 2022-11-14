###############################################################################
#### SCRIPT TO GENERATE MATRIX OF COMPLEX REARRANGEMENTS (COUNTS AND FREQ) ####
###############################################################################


#libraries
library(stringr)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE INPUT FILES #############

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v3.0.txt",sep="")

#Format
patients_type_file <- patients_type_file[,c(1,3)]
names(patients_type_file) <- c("sample_id", "risk_group")

#Exclude patients
patients_type_file <- patients_type_file[!grepl(patients_type_file$sample_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54"),]


############# OPEN AND PREPARE MUTATION FILES #############

#complex rearrangements counts
complex_svs <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/allCOUNTS_complexrearrangements_jabba+AA+MC_collapseequal_allpatients.txt",sep=""), sep="\t")
complex_svs <- complex_svs[!grepl(complex_svs$patient_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54"),]


### prepare input MATRIX
#empty matrix
MATRIX_complex_svs <- matrix(0,length(unique(complex_svs$patient_id)),length(unique(complex_svs$type)))
rownames(MATRIX_complex_svs) <- unique(complex_svs$patient_id)
colnames(MATRIX_complex_svs) <- sort(unique(complex_svs$type))
MATRIX_complex_svs <- as.data.frame(MATRIX_complex_svs)

#fill matrix
for (i in 1:dim(MATRIX_complex_svs)[1]){
	patient <- rownames(MATRIX_complex_svs)[i]
	results <- complex_svs[complex_svs$patient_id==patient,]
	if (dim(results[grep(results[,1],pattern="bfb"),])[1]>0){
		num <- results[grep(results[,1],pattern="bfb"),2]
		MATRIX_complex_svs$bfb[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="chromoplexy"),])[1]>0){
		num <- results[grep(results[,1],pattern="chromoplexy"),2]
		MATRIX_complex_svs$chromoplexy[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="chromothripsis"),])[1]>0){
		num <- results[grep(results[,1],pattern="chromothripsis"),2]
		MATRIX_complex_svs$chromothripsis[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="Clustered-rearrangement"),])[1]>0){
		num <- results[grep(results[,1],pattern="Clustered-rearrangement"),2]
		MATRIX_complex_svs$`Clustered-rearrangement`[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="Complex-non-cyclic"),])[1]>0){
		num <- results[grep(results[,1],pattern="Complex-non-cyclic"),2]
		MATRIX_complex_svs$`Complex-non-cyclic`[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="ecDNA"),])[1]>0){
		num <- results[grep(results[,1],pattern="ecDNA"),2]
		MATRIX_complex_svs$ecDNA[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="pyrgo"),])[1]>0){
		num <- results[grep(results[,1],pattern="pyrgo"),2]
		MATRIX_complex_svs$pyrgo[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="rigma"),])[1]>0){
		num <- results[grep(results[,1],pattern="rigma"),2]
		MATRIX_complex_svs$rigma[i] <- num
	}
	if (dim(results[grep(results[,1],pattern="tic"),])[1]>0){
		num <- results[grep(results[,1],pattern="tic"),2]
		MATRIX_complex_svs$tic[i] <- num
	}
}

MATRIX_complex_svs$sample_id <- rownames(MATRIX_complex_svs)

#merge with risk group info and all patients
#end up with a matrix for the 114 patients
MATRIX_complex_svs_M <- merge(MATRIX_complex_svs, patients_type_file, by="sample_id", all=TRUE)
MATRIX_complex_svs_M[is.na(MATRIX_complex_svs_M)] <- 0

#Generate frequency matrix
MATRIX_complex_svs_M_FREQ <- MATRIX_complex_svs_M

for (l in 1:dim(MATRIX_complex_svs_M_FREQ)[1]){
	MATRIX_complex_svs_M_FREQ[l,2:10] <-MATRIX_complex_svs_M_FREQ[l,2:10]/sum(MATRIX_complex_svs_M_FREQ[l,2:10])
}

MATRIX_complex_svs_M_FREQ[is.na(MATRIX_complex_svs_M_FREQ)] <- 0

### save matrix results
#counts
write.table(MATRIX_complex_svs_M,paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_COUNTS_perpatient_114patients_RISKGROUPinfo_",typeofrun,AAversion,".txt",sep=""),row.names=FALSE, quote=FALSE, sep="\t")
#freq
write.table(MATRIX_complex_svs_M_FREQ,paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_FREQ_perpatient_114patients_RISKGROUPinfo_",typeofrun,AAversion,".txt", sep=""),row.names=FALSE, quote=FALSE, sep="\t")

