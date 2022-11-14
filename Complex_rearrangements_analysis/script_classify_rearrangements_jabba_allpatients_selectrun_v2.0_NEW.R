##################################################################################
#### SCRIPT RETRIEVE COMPLEX REARRANGEMENTS FROM JABBA (JABBA CLASSIFICATION) ####
##################################################################################

### Libraries
library(gGnome)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal

## Prepare data frame
DT_JABBA_ALL_PATIENTS = NULL


#### CLASSIFY JABBA RESULTS ####

## Open patient's file
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")


## Retrieve type of rearrangements called by jabba for each patient in patients_type_file
for (p in 1:dim(patients_type_file)[1]){
  patient_id = patients_type_file[p,1]
  ## Open jabba results file
  if (file.exists(paste("/fast/scratch/users/rodrigue_c/projects/jabba/",patient_id,"/run_",typeofrun,"SVs/jabba.simple.gg.rds", sep=""))){
    gg.jabba = gG(jabba = paste("/fast/scratch/users/rodrigue_c/projects/jabba/",patient_id,"/run_",typeofrun,"SVs/jabba.simple.gg.rds", sep=""))
    ## Classify all the events
    gg.jabba = events(gg.jabba)
    if (dim(data.frame(gg.jabba$meta$events))[1]>0){
      dt.jabba <- data.frame(gg.jabba$meta$events)
      dt.jabba <- dt.jabba[dt.jabba$type!="del",]
      dt.jabba <- dt.jabba[dt.jabba$type!="dup",]
      dt.jabba <- dt.jabba[dt.jabba$type!="inv",]
      dt.jabba <- dt.jabba[dt.jabba$type!="invdup",]
      dt.jabba <- dt.jabba[dt.jabba$type!="tra",]
      dt.jabba.count <- data.frame(table(dt.jabba$type))
      if (dim(dt.jabba)[1]>0){
        names(dt.jabba.count) <- c("type", "Freq")
        dt.jabba.count$patient_id <- patient_id
      }
      write.table(dt.jabba.count, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/complex_rearrang_counts/jabba_run_",typeofrun,"_all_rearrangements_classificationTABLE_patient",patient_id,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      write.table(dt.jabba, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/complex_rearrang_info/jabba_run_FILTEREDATLEAST2CALLERS_all_COMPLEX_rearrangementsINFO_patient",patient_id,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      print(paste(patient_id, " completed",sep=""))
    } else if (dim(data.frame(gg.jabba$meta$events))[1]==0){
      columns <- c("Freq")
      dt.jabba <- data.frame(matrix(nrow = 0, ncol = length(columns)))
      colnames(dt.jabba) = columns
      write.table(dt.jabba, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/complex_rearrang_counts/jabba_run_",typeofrun,"_all_rearrangements_classificationTABLE_patient",patient_id,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      write.table(dt.jabba, paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/complex_rearrang_info/jabba_run_FILTEREDATLEAST2CALLERS_all_COMPLEX_rearrangementsINFO_patient",patient_id,".txt",sep=""), sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
      print(paste(patient_id, " completed",sep=""))
    }
  }
}
