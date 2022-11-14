#############################################################
###### COLLAPSE JABBA+AA INFO WITH SHATTERSEEK RESULTS ######
#############################################################


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
list_GRcoords_SS = NULL
list_CRtypes_SS = NULL
list_GRcoords_JBAA = NULL
list_CRtypes_JBAA = NULL
indx_SS = NULL
indx_JBAA = NULL
LIST_ALL <- list(vector(mode = "list", length = 1))
L <- 0

#working directory
wd <- "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results"
dir.create(file.path(paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,sep="")), showWarnings = FALSE)


##### OPEN INPUT FILES #####

#Open Manual classification results for the patient
SS_file <- dir("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/shatterseek_chromothripsis/results_ATLEAST2CALLERS", pattern= patient_id)

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

#### SHATTERSEEK FILES
if (exists("SS_file")){
  if (length(SS_file)[1]>0){
    #Iterate for each rearrangement file, listing the coordinates
    for (a in 1:length(SS_file)[1]){
      typeSS <- "chromothripsis"
      #Open AA info
      SS_file_2 <- read.delim(paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/shatterseek_chromothripsis/results_ATLEAST2CALLERS/",SS_file[a],sep=""))
      SS_file_2$chrom <- unlist(lapply(SS_file_2$chrom, function(x) paste('chr', x,sep=""))) 
      SS_file_2 <- SS_file_2[,1:3]
      names(SS_file_2) <- c("chr", "start", "end")
      GR_coords_SS <- makeGRangesFromDataFrame(SS_file_2, keep.extra.columns=TRUE, ignore.strand=TRUE) #GRanges object to do the intersect
      list_GRcoords_SS <- c(list_GRcoords_SS, GR_coords_SS) #same lists but for AA
      list_CRtypes_SS <- c(list_CRtypes_SS, typeSS)
    }
  }
}


##### INTERSECT THE AMPLICONS DETECTED BY THE DIFFERENT ALGORITHMS #####

#If we only have results from 1 algorithm
if (is.null(list_GRcoords_SS) && !is.null(list_GRcoords_JBAA)){
  stop(paste(patient_id,": Shatterseek list is empty - NO CHANGES",sep="")) #break script in case there are no JABBA and AA results
} else if (!is.null(list_GRcoords_SS) && !is.null(list_GRcoords_JBAA)){
  indx_SS = NULL
  indx_JBAA = NULL
  for (z in 1:length(list_GRcoords_SS)){
    z.file=tempfile()
    write.table(data.frame(list_GRcoords_SS[[z]])[,1:3],file=z.file, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    for (z_1 in 1:length(list_GRcoords_JBAA)){
      z_1.file=tempfile()
      write.table(data.frame(list_GRcoords_JBAA[[z_1]])[,1:3],file=z_1.file, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
      INTERSECT_a <- unique(system(paste("intersectBed -a", z.file, "-b", z_1.file, "-wa", "-r"), intern=TRUE))
      INTERSECT_b <- unique(system(paste("intersectBed -a", z_1.file, "-b", z.file, "-wa", "-r"), intern=TRUE))
      if ((length(INTERSECT_a)/length(list_GRcoords_JBAA[[z]])>=0.75)&&(length(INTERSECT_b)/length(list_GRcoords_SS[[z_1]])>=0.75)&&(length(list_GRcoords_SS[[z_1]])/length(list_GRcoords_JBAA[[z]])>=0.25)&&(length(list_GRcoords_JBAA[[z]])/length(list_GRcoords_SS[[z_1]])>=0.25)){
        indx_SS <- c(indx_SS,z)
        indx_JBAA <- c(indx_JBAA,z_1)
      }
    }
  }
} else if (is.null(list_GRcoords_SS) && is.null(list_GRcoords_JBAA)){
  stop(paste(patient_id,": Both lists are empty - NO RESULTS",sep="")) #break script in case there are no JABBA and AA results
}

indx_SS <- unique(indx_SS)
indx_JBAA <- unique(indx_JBAA)

#We substract the CR that have been merged
if (!is.null(indx_SS) && !is.null(indx_JBAA)){
  system(paste("rm ", wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,"/",list_files_JBAA[indx_JBAA],sep=""))
  DF <- data.frame(list_GRcoords_SS[[indx_SS]])[,c(1:3)]
  names(DF) <- c("chr", "start", "end")
  DF$SS_type <- list_CRtypes_SS[[indx_SS]]
  DF$JBAA_type <- "none"
  LIST_ALL <- c(LIST_ALL, list(DF))
} else if (is.null(indx_SS) && is.null(indx_JBAA)){
  DF <- data.frame(list_GRcoords_SS[[1]])[,c(1:3)]
  names(DF) <- c("chr", "start", "end")
  DF$SS_type <- list_CRtypes_SS[[1]]
  DF$JBAA_type <- "none"
  LIST_ALL <- c(LIST_ALL, list(DF))
}

#Remove first entry on LIST_ALL which is empty
LIST_ALL <- LIST_ALL[-1]


##### WRITE RESULTS #####
#Iterate through the list
if (length(LIST_ALL)>0){
  for (l in 1:length(LIST_ALL)){
    DF <- as.data.frame(LIST_ALL[[l]])
    type_DF <- unique(c(DF$SS_type))
    write.table(DF[,c(1:3)], paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id, "/", patient_id, "_", type_DF, "_",(L+l),".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
  }
}

print(paste(patient_id, ": Completed", sep=""))
