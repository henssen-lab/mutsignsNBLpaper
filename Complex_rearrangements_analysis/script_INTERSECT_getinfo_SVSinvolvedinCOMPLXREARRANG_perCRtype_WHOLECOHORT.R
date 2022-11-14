#####################################################################################
###### SCRIPT INTERSECT SVS AND COMPLEX REARRANGEMENTS 2ND PLOT REARRANG TYPES ######
#####################################################################################

### Library
library("bedr")
library("reshape2")

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal

### Initialize DFs
num_total_SVs <- as.numeric(system("cat /fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/*/*.ALLSVS.jabbainput.TIER1+2.bedpe | grep -vP \"start|NBL47|NBL49|NBL50|NBL53|NBL54|CB2044\" | wc -l", intern=TRUE))

all_patients_DF = NULL
FILES_all = NULL
ALL_coords_CR = NULL


############# OPEN AND PREPARE INPUT FILES #############

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v3.0.txt",sep="")

#Format
patients_type_file <- patients_type_file[,c(1,3)]
names(patients_type_file) <- c("sample_id", "risk_group")

#Exclude patients
patients_type_file <- patients_type_file[grep(patients_type_file$sample_id,pattern="CB2044|NBL47|NBL53|NBL54|NBL49|NBL50",invert=TRUE),]


############# GENERATE 2 SUMMARY TABLES : 1. ALL REGIONS OF CR 2. PROPORTION OF SVS INVOLVED IN CR #############

#List of directories
list_d <- grep(dir(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords",sep="")), pattern='CB2044|NBL47|NBL53|NBL54|NBL49|NBL50', invert=TRUE, value=TRUE)

#Open files for each patient
for (d in 1:length(list_d)){
    list_f <- grep(dir(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", list_d[d], sep="")), pattern="Linear-amplification", invert=TRUE, value=TRUE)
    if (length(list_f)>0){
        all_lines_CR = NULL
        patient_id <- unlist(strsplit(list_f, "_"))[1]
        patient_type <- patients_type_file[patients_type_file$sample_id == patient_id,2]
        #Complex rearrangements coord
        for (f in 1:length(list_f)){
            file <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", list_d[d], "/", list_f[f], sep=""), sep="")
            CR_type <- unlist(strsplit(list_f[f], "_"))[2]
            ##
            #### LOOK AT THE INTERSECTION BETWEEN SVS AND CR PER PATIENT
            #### GET PROPORTION OF SVS INVOLVED IN COMPLEX REARRANGEMENTS
            #Format complex rearrangements coords into bedr format
            file_bedr <- paste(file$chr,":",file$start,"-",file$end,sep="")
            file_bedr.sort <- bedr.sort.region(file_bedr, method = "natural")
            #SVs
            file_SV <- read.table(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/",patient_id,"/",patient_id,".ALLSVS.jabbainput.TIER1+2.bedpe", sep=""), sep="\t")
            #Format SV coords into bedr format
            #Split each SV in 2 to also include translocations
            #I accept a window of +/-200bp around the bkps
            file_SV_bedr1 <- paste("chr",file_SV$V1,":",file_SV$V2,"-",file_SV$V3,sep="")
            file_SV_bedr2 <- paste("chr",file_SV$V4,":",file_SV$V5,"-",file_SV$V6,sep="")
            file_SV_bedr <- c(file_SV_bedr1, file_SV_bedr2)
            file_SV_bedr.sort <- bedr.sort.region(file_SV_bedr, method = "natural")
            interesect_file_SV <- unique(file_SV_bedr.sort[in.region(file_SV_bedr.sort, file_bedr.sort)])
            ### PROPORTION OF SVS
            #Get the proportion of SVs involved in complex rearrangements
            if (length(interesect_file_SV)[1]>0){
                for (x in 1:length(interesect_file_SV)[1]){
                    pos1 <- unlist(strsplit(interesect_file_SV[x],split="-"))[2]
                    lines_pos1 <- file_SV[with(file_SV, grepl(as.numeric(pos1), paste(V2,V3,V5,V6))),] #Retrieve original SV.
                    all_lines_CR <- unique(rbind(all_lines_CR,lines_pos1))
                }
                num_SVs_patient <- dim(all_lines_CR)[1]
                DF_patient <- data.frame(patient_id, patient_type,CR_type,num_SVs_patient)
                all_patients_DF <- rbind(all_patients_DF, DF_patient)
            } else {
                num_SVs_patient <- 0
                DF_patient <- data.frame(patient_id, patient_type,CR_type,num_SVs_patient)
                all_patients_DF <- rbind(all_patients_DF, DF_patient)
            }
        }
    }
}

### GET FINAL TABLE
matrix_CR_RG <- as.data.frame.matrix(table(all_patients_DF$patient_type, all_patients_DF$CR_type))

for (i in 1:dim(matrix_CR_RG)[1]){
    for (x in 1:dim(matrix_CR_RG)[2]){
        matrix_CR_RG[i,x] <- sum(all_patients_DF[all_patients_DF$patient_type==rownames(matrix_CR_RG)[i]&all_patients_DF$CR_type==colnames(matrix_CR_RG)[x],4])
    }
}

#COUNTS FOR TOTAL SAMPLES, not risk groups
#final data frame
df_CR_RG <- data.frame(colSums(matrix_CR_RG))
df_CR_RG[rownames(df_CR_RG)=="bfb",] <- 2 #Adding 2 BFB found through visual inspection
df_CR_RG$type_CR <- rownames(df_CR_RG)
df_CR_RG <- df_CR_RG[,c(2,1)]
names(df_CR_RG) <- c("type_CR", "num_SVs_CR")
rownames(df_CR_RG) <- NULL

#calculate percentage of SVs involved in each rearrangement type
for (p in 1:dim(df_CR_RG)[1]){
    df_CR_RG$percent_SVs_CR[p] <- (df_CR_RG$num_SVs_CR[p]*100)/num_total_SVs
}

### SAVE RESULTS

#Save proportion of svs
write.table(df_CR_RG, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/INTERSECT_SVS_complexrearrangements/COUNTS_SVS_implicated_complexrearrangements_PERCOMPLXREARRANGTYPE_WHOLECOHORT.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
