#########################################################
##### SCRIPT TO COUNT NUM OF COMPLEX REARRANGEMENTS #####
#####     JABBA + AA + MC RESULTS COLLAPSED         #####
#########################################################

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal

##Initialize final table
ALL_PATIENTS_count <- NULL

##Set working directory with all folders (one for each patient) containing the info about rearrangements
setwd(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords",sep=""))

#Patients with rearrangements
patient_f <- dir(, pattern ="B")


#### GET COUNTS ####

#Retrieve counts for each patient
for (c in 1:length(patient_f)){
  patient_id <- patient_f[c]
  patient_R <- data.frame(dir(paste(getwd(),"/",patient_id, sep="")))
  names(patient_R) <- "filename"
  for (r in 1:dim(patient_R)[1]){
    type <- paste(unlist(strsplit( unlist(strsplit(patient_R$filename[r], "_"))[2]," ")),collapse="-")
    if (type =="dm" || type == "cpxdm" ){
      type <- "ecDNA"
    } else if (type == "BFB"){
      type <- "bfb"
    }
    patient_R$type[r] <- type
  }
  patient_T <- table(patient_R$type)
  patient_FINAL <- data.frame(patient_T)
  names(patient_FINAL) <- c("type", "Freq")
  patient_FINAL$patient_id <- patient_id
  ALL_PATIENTS_count <- rbind(ALL_PATIENTS_count, patient_FINAL)
}

#At the moment get rid of the Linear amplifications (don't take this as complex rearrangements)
ALL_PATIENTS_count <- ALL_PATIENTS_count[ALL_PATIENTS_count$type!="Linear-amplification",]

#Save table
write.table(ALL_PATIENTS_count, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/allCOUNTS_complexrearrangements_jabba+AA+MC_collapseequal_allpatients.txt",sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
