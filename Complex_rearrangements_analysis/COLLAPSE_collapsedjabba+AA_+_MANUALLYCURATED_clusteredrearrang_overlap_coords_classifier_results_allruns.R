#################################################################################
###### COLLAPSE JABBA+AA INFO WITH MANUAL CURATED CLUSTERED REARRANGEMENTS ######
#################################################################################


###Libraries
suppressMessages(library(GenomicRanges))


### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal
AAversion2 <- args[2] # select between forced, or normal
patient_id <- args[3]


##### PREPARE FILES AND WD

#files
list_GRcoords_MC = NULL
list_CRtypes_MC = NULL
list_GRcoords_JBAA = NULL
list_CRtypes_JBAA = NULL
indx_MC = NULL
indx_JBAA = NULL
LIST_ALL <- list(vector(mode = "list", length = 1))
L <- 0

#working directory
wd <- "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results"
dir.create(file.path(paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,sep="")), showWarnings = FALSE)


##### OPEN INPUT FILES #####

#Open Manual classification results for the patient
MC_file <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/manual_SVs_classification_clusters/CLUSTERED_SVs_10bkps_10Mbwindow_ALLpatients_SVsATLEAST2CALLERS_coordinatesofclusters.txt",sep="\t")
MC_file <- MC_file[MC_file$patient_id==patient_id,]

#Open JABBA+Amplicon Architect amplicons of the same patient
list_files_JBAA <- dir(paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,sep=""), pattern= patient_id)
list_files_JBAA_filt <- list_files_JBAA #since I'm interested in reclassifying all SVs, I don't filter unknown amplicons


##### GENERATE GRANGES LISTS FOR THE 2 ALGORITHMS #####

#### JABBA+AMPLICON ARCHITECT COLLAPSED FILES
if (exists("list_files_JBAA_filt")){
  L <- length(list_files_JBAA_filt)
  if (length(list_files_JBAA_filt)>0){
    #Iterate for each rearrangement file, listing the coordinates
    for (a in 1:length(list_files_JBAA_filt)){
      typeAA <- unlist(strsplit(list_files_JBAA_filt[a],"_"))[2]
      #Open AA info
      AA_file <- read.delim(paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,"/", list_files_JBAA_filt[a],sep=""))
      names(AA_file) <- c("chr", "start", "end")
      GR_coords_AA <- makeGRangesFromDataFrame(AA_file, keep.extra.columns=TRUE, ignore.strand=TRUE) #GRanges object to do the intersect
      list_GRcoords_JBAA <- c(list_GRcoords_JBAA, GR_coords_AA) #same lists but for AA
      list_CRtypes_JBAA <- c(list_CRtypes_JBAA, typeAA)
    }
  }
}

#### MANUAL CURATION FILES
if (exists("MC_file")){
  if (dim(MC_file)[1]>0){
    #Iterate for each rearrangement file, listing the coordinates
    for (a in 1:dim(MC_file)[1]){
      typeMC <- "Clustered-rearrangement"
      #Open AA info
      MC_file_2 <- MC_file[a,c(1:3)]
      names(MC_file_2) <- c("chr", "start", "end")
      GR_coords_MC <- makeGRangesFromDataFrame(MC_file_2, keep.extra.columns=TRUE, ignore.strand=TRUE) #GRanges object to do the intersect
      list_GRcoords_MC <- c(list_GRcoords_MC, GR_coords_MC) #same lists but for AA
      list_CRtypes_MC <- c(list_CRtypes_MC, typeMC)
    }
  }
}


##### INTERSECT THE AMPLICONS DETECTED BY THE DIFFERENT ALGORITHMS #####

#If we only have results from 1 algorithm
if (is.null(list_GRcoords_MC) && !is.null(list_GRcoords_JBAA)){
  stop(paste(patient_id,": Manual Curation list is empty - NO CHANGES",sep="")) #break script in case there are no JABBA and AA results
} else if (!is.null(list_GRcoords_MC) && !is.null(list_GRcoords_JBAA)){
  indx_MC = NULL
  indx_JBAA = NULL
  for (z in 1:length(list_GRcoords_MC)){
    for (z_1 in 1:length(list_GRcoords_JBAA)){
      INTERSECT <- suppressWarnings(intersect(list_GRcoords_JBAA[[z_1]], list_GRcoords_MC[[z]]))
      if ((length(INTERSECT)>0)){
        indx_MC <- c(indx_MC,z)
        indx_JBAA <- c(indx_JBAA,z_1)
      }
    }
  }
} else if (is.null(list_GRcoords_MC) && is.null(list_GRcoords_JBAA)){
  stop(paste(patient_id,": Both lists are empty - NO RESULTS",sep="")) #break script in case there are no JABBA and AA results
}


#We substract the CR that have been merged
if (!is.null(indx_MC) && !is.null(indx_JBAA)){
  list_GRcoords_MC <- list_GRcoords_MC[-indx_MC]
  list_CRtypes_MC <- list_CRtypes_MC[-indx_MC]
  list_GRcoords_JBAA <- list_GRcoords_JBAA[-indx_JBAA]
  list_CRtypes_JBAA <- list_CRtypes_JBAA[-indx_JBAA]
}

#We fish the rest of the CR that are not merged
if (length(list_GRcoords_MC)>0 && length(list_GRcoords_JBAA)==0){
  for (z in 1:length(list_GRcoords_MC)){
    DF <- data.frame(list_GRcoords_MC[[z]])[,c(1:3)]
    names(DF) <- c("chr", "start", "end")
    DF$MC_type <- list_CRtypes_MC[[z]]
    DF$JBAA_type <- "none"
    LIST_ALL <- c(LIST_ALL, list(DF))
  }
} else if (length(list_GRcoords_MC)==0 && length(list_GRcoords_JBAA)>0){
  print("nada que hacer")
} else if (length(list_GRcoords_MC)>0 && length(list_GRcoords_JBAA)>0){
  for (z in 1:length(list_GRcoords_MC)){
    DF <- data.frame(list_GRcoords_MC[[z]])[,c(1:3)]
    names(DF) <- c("chr", "start", "end")
    DF$MC_type <- list_CRtypes_MC[[z]]
    DF$JBAA_type <- "none"
    LIST_ALL <- c(LIST_ALL, list(DF))
  }
}

#Remove first entry on LIST_ALL which is empty
LIST_ALL <- LIST_ALL[-1]


##### WRITE RESULTS #####
#Iterate through the list
if (length(LIST_ALL)>0){
  for (l in 1:length(LIST_ALL)){
    DF <- as.data.frame(LIST_ALL[[l]])
    type_DF <- unique(c(DF$MC_type))
    write.table(DF[,c(1:3)], paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id, "/", patient_id, "_", type_DF, "_",(L+l),".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  }
}

print(paste(patient_id, ": Completed", sep=""))
