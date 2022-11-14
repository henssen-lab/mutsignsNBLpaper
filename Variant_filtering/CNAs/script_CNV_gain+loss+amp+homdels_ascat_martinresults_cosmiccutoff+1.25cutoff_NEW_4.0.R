###############################################################
#### SCRIPT TO FILTER ALL CNVS FROM ASCAT (MARTIN RESULTS) ####
###############################################################


############# DECLARE LISTS #############

MUT_all_PATIENTS = NULL


############# OPEN AND PREPARE MUTATION FILES #############

#open file with patient_id, risk group and cohort
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")

patients_type_file <- patients_type_file[patients_type_file$V1!="NBL31",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL36",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL61",]

#Create factor with patient id
patients_list <- patients_type_file[,1]

#open CNA files
cna_martin_results <- read.delim("/fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/ascat_results_hg19_martin/ALL.cnseg.tsv",sep="")

#Change format and retrieve ploidy status
for (i in 1:length(patients_list)){
    patient_id <- patients_list[i]
    res_new2 <- cna_martin_results[cna_martin_results$sample==patient_id,]
    ploidy_file <- read.delim(paste("/fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/ascat_results_hg19_martin/CN+purity_ascat_martin_perpatient_nologRinfo/",patient_id,".purity.tsv", sep=""),sep="")
    if (dim(res_new2)[1]>0){ 
      res_new2$ploidy <- ploidy_file$ploidy
      res_new2$length <- abs(as.numeric(res_new2[,4])-as.numeric(res_new2[,3]))
      res_new2$info <- paste("ASCAT=",res_new2[,5],"/",res_new2[,6],",",res_new2$mean_logr,",",res_new2$ploidy,",cnv,",res_new2$length, sep="")
      res_new2$flag <- "."
      for (a in 1:dim(res_new2)[1]){
        if (res_new2$ploidy[a]<=2.7 & (as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))>=5 & log2((as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))/as.numeric(res_new2$ploidy[a]))>1.25){
          res_new2$mut_type[a] <- "amplification"
        } else if(res_new2$ploidy[a]>2.7 & (as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))>=9 & log2((as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))/as.numeric(res_new2$ploidy[a]))>1.25){
        res_new2$mut_type[a] <- "amplification"
      } else if(res_new2$ploidy[a]<=2.7 & (as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))==0){
        res_new2$mut_type[a] <- "homozygous_loss"
      } else if(res_new2$ploidy[a]>2.7 & (as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))<(as.numeric(res_new2$ploidy[a])-2.7)){
        res_new2$mut_type[a] <- "homozygous_loss"
      } else if(log2((as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))/as.numeric(res_new2$ploidy[a]))<(-0.3)){
        res_new2$mut_type[a] <- "loss"
      } else if(log2((as.numeric(res_new2[a,5])+as.numeric(res_new2[a,6]))/as.numeric(res_new2$ploidy[a]))>0.3){
        res_new2$mut_type[a] <- "gain"
      } else {
        res_new2$mut_type[a] <- "neutral"
      }
      }
      res_new_cna <- res_new2[res_new2$mut_type!="neutral",]
      if (dim(res_new_cna)[1]>0){
        res_new_format <- res_new_cna[,c(1,2,3,2,4,10,12,11,9)]
        MUT_all_PATIENTS <- rbind(MUT_all_PATIENTS, res_new_format) 
      }
  }
}

write.table(MUT_all_PATIENTS, "RESULTS_all_patients_berlin+peifer_cohorts_CNV_GAINS+LOSSES+AMP+HOMDELS_ASCAT_martin_results_cutoffploidy+cutoff1.25_ALLpatients_21.04.21.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
