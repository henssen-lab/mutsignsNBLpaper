#############################################################
###### SCRIPT INTERSECT SVS AND COMPLEX REARRANGEMENTS ######
#############################################################

### Library
library("bedr")

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal

### Initialize DFs
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
        FILES_all = NULL
        for (f in 1:length(list_f)){
            file <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", list_d[d], "/", list_f[f], sep=""), sep="")
            CR_type <- unlist(strsplit(list_f[f], "_"))[2]
            file$patient_id <- patient_id
            file$patient_type <- patient_type
            file$CR_type <- CR_type
            FILES_all <- rbind(FILES_all, file)
        }
        #### TABLE ALL REGIONS OF CR FOR DOWNSTREAM ANALYSES
        #Save all rearrangements wih all coordinates for all patients in a single table
        ALL_coords_CR <- rbind(ALL_coords_CR, FILES_all)
        ##
        #### LOOK AT THE INTERSECTION BETWEEN SVS AND CR PER PATIENT
        #### GET PROPORTION OF SVS INVOLVED IN COMPLEX REARRANGEMENTS
        #Format complex rearrangements coords into bedr format
        file_bedr <- paste(FILES_all$chr,":",FILES_all$start,"-",FILES_all$end,sep="")
        file_bedr.sort <- bedr.sort.region(file_bedr, method = "natural")
        #SVs
        file_SV <- read.table(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/",patient_id,"/",patient_id,".ALLSVS.jabbainput.TIER1+2.bedpe", sep=""), sep="\t")
        #Format SV coords into bedr format
        #Split each SV in 2 to also include translocations
        file_SV_bedr1 <- paste("chr",file_SV$V1,":",file_SV$V2,"-",file_SV$V3,sep="")
        file_SV_bedr2 <- paste("chr",file_SV$V4,":",file_SV$V5,"-",file_SV$V6,sep="")
        file_SV_bedr <- c(file_SV_bedr1, file_SV_bedr2)
        file_SV_bedr.sort <- bedr.sort.region(file_SV_bedr, method = "natural")
        interesect_file_SV <- unique(file_SV_bedr.sort[in.region(file_SV_bedr.sort, file_bedr.sort)])
        #Total number of SVs
        num_total_SVs <- dim(file_SV)[1]
        ### PROPORTION OF SVS
        #Get the proportion of SVs involved in complex rearrangements
        if (length(interesect_file_SV)[1]>0){
            for (x in 1:length(interesect_file_SV)[1]){
                pos1 <- unlist(strsplit(interesect_file_SV[x],split="-"))[2]
                lines_pos1 <- file_SV[with(file_SV, grepl(pos1, paste(V2,V3,V5,V6))),]
                all_lines_CR <- unique(rbind(all_lines_CR,lines_pos1))
            }
            num_SVs_patient <- dim(all_lines_CR)[1]
            proportion_SVs_patient <- (num_SVs_patient/num_total_SVs)*100
            DF_patient <- data.frame(patient_id, patient_type, proportion_SVs_patient)
            all_patients_DF <- rbind(all_patients_DF, DF_patient)
        } else {
            num_SVs_patient <- 0
            proportion_SVs_patient <- (num_SVs_patient/num_total_SVs)*100
            DF_patient <- data.frame(patient_id, patient_type, proportion_SVs_patient)
            all_patients_DF <- rbind(all_patients_DF, DF_patient)
        }
    }
}

### SAVE RESULTS

#Save proportion of svs
write.table(all_patients_DF, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/INTERSECT_SVS_complexrearrangements/COUNTS_SVS_implicated_complexrearrangements_PERPATIENT_PER_RISK_GROUP.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)

#Save all regions in complex rearrangements
write.table(ALL_coords_CR, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/table_ALLregions_ALLcomplexrearrangements_114patients_infopatient_inforiskgroup.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
