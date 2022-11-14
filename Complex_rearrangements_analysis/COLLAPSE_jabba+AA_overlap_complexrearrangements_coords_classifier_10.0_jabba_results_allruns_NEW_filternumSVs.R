#############################################################
###### COLLAPSE JABBA + AA INFO COMPLEX REARRANGEMENTS ######
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
list_GRcoords_JABBA = NULL
list_CRtypes_JABBA = NULL
list_GRcoords_AA = NULL
list_CRtypes_AA = NULL
indx_JABBA = NULL
indx_AA = NULL
LIST_ALL <- list(vector(mode = "list", length = 1))

#working directory
wd <- "/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results"
dir.create(file.path(paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id,sep="")), showWarnings = FALSE)


##### OPEN INPUT FILES #####

#Number SVs patient
tmp <- system(paste("wc -l /data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/",patient_id,"/",patient_id,".ALLSVS.jabbainput.TIER1+2.bedpe | cut -d\" \" -f1",sep=""), intern=TRUE)
if (tmp==0){
  num_SVs <- 0
} else {
  num_SVs <- dim(read.table(paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/",patient_id,"/",patient_id,".ALLSVS.jabbainput.TIER1+2.bedpe",sep="")))[1]
}

#Open JABBA file for patient_id
J_F <- paste(wd,"/jabba_results_",typeofrun,AAversion,"/complex_rearrang_info/jabba_run_",typeofrun,"_all_COMPLEX_rearrangementsINFO_patient",patient_id,".txt",sep="")
if (file.exists(J_F)){
  JABBA_file <- read.delim(J_F)
} else {
  JABBA_file = NULL
}

#Open Amplicon Architect amplicons of the same patient
list_files_AA <- dir(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/amplicon_architect_classif/classification_",AAversion2,"/input_file_classif_amplicons_AA_classification_bed_files",sep=""), pattern= patient_id)
#list_files_AA_filt <- list_files_AA[!grepl('unknown|Linear amplification', list_files_AA)] #not interested in linear or unknown amplicons
list_files_AA_filt <- list_files_AA #since I'm interested in reclassifying all SVs, I don't filter unknown amplicons


##### GENERATE GRANGES LISTS FOR THE 2 ALGORITHMS #####

### JABBA
if (exists("JABBA_file")&&!is.null(JABBA_file)){
  if (dim(JABBA_file)[1]>0){
    JABBA_file_NEW <- JABBA_file[,c("type","footprint")] #we only need the type of rearrangement and the coordinates of each region
    JABBA_file_NEW$patient <- patient_id
    #Iterate for each rearrangement listing the coordinates
    for (x in 1:dim(JABBA_file_NEW)[1]){
      typeJABBA <- JABBA_file_NEW[x,1]
      coords <- data.frame(unlist(strsplit(as.character(JABBA_file_NEW[x,2]),";|,")))
      names(coords) <- "coords"
      split1 <- data.frame(do.call('rbind', strsplit(as.character(coords$coords),':',fixed=TRUE)))
      split2 <- data.frame(do.call('rbind', strsplit(as.character(split1$X2),'-',fixed=TRUE)))
      split3 <- data.frame(do.call('rbind', strsplit(as.character(split2$X2),'+',fixed=TRUE)))
      coords_NEW <- cbind (data.frame(split1[,1]), data.frame(split2[,1]), data.frame(split3[,1]))
      names(coords_NEW) <- c("chr", "start", "end")
      coords_NEW$chr <- paste("chr",coords_NEW$chr,sep="")
      GR_coords <- makeGRangesFromDataFrame(coords_NEW, keep.extra.columns=TRUE, ignore.strand=TRUE) #GRanges object to do the intersect
      list_GRcoords_JABBA <- c(list_GRcoords_JABBA, GR_coords) #list with all coordinates x rearrangement
      list_CRtypes_JABBA <- c(list_CRtypes_JABBA, typeJABBA) #list with all types x rearrangement
    }
  }
}

#### AMPLICON ARCHITECT
if (exists("list_files_AA_filt")){
  if (length(list_files_AA_filt)>0){
    #Iterate for each rearrangement file, listing the coordinates
    for (a in 1:length(list_files_AA_filt)){
      typeAA <- unlist(strsplit(unlist(strsplit(list_files_AA_filt[a],"amplicon[0-9]+_"))[2], "_[0-9]+_intervals"))[1]
      #Open AA info
      AA_file <- read.table(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/amplicon_architect_classif/classification_",AAversion2,"/input_file_classif_amplicons_AA_classification_bed_files/", list_files_AA_filt[a],sep=""))
      names(AA_file) <- c("chr", "start", "end")
      AA_file$chr <- paste("chr",AA_file$chr,sep="")
      if (dim(AA_file)[1]<=num_SVs*1.10){ #Remove all AA results that have more SVs than the total num of SVs found in our calls - I let a 10% margin
        GR_coords_AA <- makeGRangesFromDataFrame(AA_file, keep.extra.columns=TRUE, ignore.strand=TRUE) #GRanges object to do the intersect
        list_GRcoords_AA <- c(list_GRcoords_AA, GR_coords_AA) #same lists but for AA
        list_CRtypes_AA <- c(list_CRtypes_AA, typeAA)
      }
    }
  }
}


##### INTERSECT THE AMPLICONS DETECTED BY THE DIFFERENT ALGORITHMS #####

#If we have results from 2 algorithms
if (!is.null(list_GRcoords_JABBA) && !is.null(list_GRcoords_AA)){
  indx_JABBA = NULL
  indx_AA = NULL
  for (z in 1:length(list_GRcoords_JABBA)){
    z.file=tempfile()
    write.table(data.frame(list_GRcoords_JABBA[[z]])[,1:3],file=z.file, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
    for (z_1 in 1:length(list_GRcoords_AA)){
      z_1.file=tempfile()
      write.table(data.frame(list_GRcoords_AA[[z_1]])[,1:3],file=z_1.file, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
      INTERSECT_a <- unique(system(paste("intersectBed -a", z.file, "-b", z_1.file, "-wa", "-r"), intern=TRUE))
      INTERSECT_b <- unique(system(paste("intersectBed -a", z_1.file, "-b", z.file, "-wa", "-r"), intern=TRUE))
      if ((length(INTERSECT_a)/length(list_GRcoords_JABBA[[z]])>=0.75)&&(length(INTERSECT_b)/length(list_GRcoords_AA[[z_1]])>=0.75)&&(length(list_GRcoords_AA[[z_1]])/length(list_GRcoords_JABBA[[z]])>=0.25)&&(length(list_GRcoords_JABBA[[z]])/length(list_GRcoords_AA[[z_1]])>=0.25)){
        DF <- data.frame(list_GRcoords_AA[[z_1]])[,c(1:3)]
        names(DF) <- c("chr", "start", "end")
        DF$JABBA_type <- list_CRtypes_JABBA[[z]]
        DF$AA_type <- paste(list_CRtypes_AA[[z_1]], sep="_")
        LIST_ALL <- c(LIST_ALL, list(DF))
        indx_JABBA <- c(indx_JABBA,z)
        indx_AA <- c(indx_AA,z_1)
      }
    }
  }
} else if (is.null(list_GRcoords_JABBA) && is.null(list_GRcoords_AA)){
  stop(paste(patient_id,": Both lists are empty - NO RESULTS",sep="")) #break script in case there are no JABBA and AA results
}

indx_JABBA <- unique(indx_JABBA)
indx_AA <- unique(indx_AA)

#We substract the CR that have been merged
if (!is.null(indx_JABBA) && !is.null(indx_AA)){
  list_GRcoords_JABBA <- list_GRcoords_JABBA[-indx_JABBA]
  list_CRtypes_JABBA <- list_CRtypes_JABBA[-indx_JABBA]
  list_GRcoords_AA <- list_GRcoords_AA[-indx_AA]
  list_CRtypes_AA <- list_CRtypes_AA[-indx_AA]
}

#We fish the rest of the CR that are not merged
if (length(list_GRcoords_JABBA)>0 && length(list_GRcoords_AA)==0){
  for (z in 1:length(list_GRcoords_JABBA)){
    DF <- data.frame(list_GRcoords_JABBA[[z]])[,c(1:3)]
    names(DF) <- c("chr", "start", "end")
    DF$JABBA_type <- list_CRtypes_JABBA[[z]]
    DF$AA_type <- "none"
    LIST_ALL <- c(LIST_ALL, list(DF))
  }
} else if (length(list_GRcoords_JABBA)==0 && length(list_GRcoords_AA)>0){
  for (z in 1:length(list_GRcoords_AA)){
    DF <- data.frame(list_GRcoords_AA[[z]])[,c(1:3)]
    names(DF) <- c("chr", "start", "end")
    DF$JABBA_type <- "none"
    DF$AA_type <- paste(list_CRtypes_AA[[z]], sep="_")
    LIST_ALL <- c(LIST_ALL, list(DF))
  }
} else if (length(list_GRcoords_JABBA)>0 && length(list_GRcoords_AA)>0){
  for (z in 1:length(list_GRcoords_JABBA)){
    DF <- data.frame(list_GRcoords_JABBA[[z]])[,c(1:3)]
    names(DF) <- c("chr", "start", "end")
    DF$AA_type <- "none"
    DF$JABBA_type <- list_CRtypes_JABBA[[z]]
    LIST_ALL <- c(LIST_ALL, list(DF))
  }
  for (z_1 in 1:length(list_GRcoords_AA)){
    DF_AA <- data.frame(list_GRcoords_AA[[z_1]])[,c(1:3)]
    names(DF_AA) <- c("chr", "start", "end")
    DF_AA$AA_type <- paste(list_CRtypes_AA[[z_1]], sep="_")
    DF_AA$JABBA_type <- "none"
    LIST_ALL <- c(LIST_ALL, list(DF_AA))
  }
}

#Remove first entry on LIST_ALL which is empty
LIST_ALL <- LIST_ALL[-1]

##### WRITE RESULTS #####
#Iterate through the list
for (l in 1:length(LIST_ALL)){
  DF <- as.data.frame(LIST_ALL[[l]])
  type_DF <- unique(c(DF$JABBA_type, DF$AA_type))
  #select the CR type if multiple
  if (length(type_DF)==1 && type_DF=="unknown"){
    type_DF <- "Complex-non-cyclic"
  } else if (length(type_DF)>1){
    if (unique(DF$AA_type)=="ecDNA"||unique(DF$JABBA_type)=="dm"){
      type_DF <- "ecDNA"
    } else if (unique(DF$AA_type)=="BFB"){
      type_DF <- "bfb"
    } else if (unique(DF$JABBA_type)=="none" && unique(DF$AA_type)=="unknown"){
      type_DF <- "Complex-non-cyclic"
    } else if (unique(DF$JABBA_type)=="none"){
      type_DF <- gsub(" ","-",unique(DF$AA_type))
    } else if (unique(DF$AA_type)=="none"){
      type_DF <- gsub(" ","-",unique(DF$JABBA_type))
    } else {
      type_DF <- gsub(" ","-",unique(DF$JABBA_type))
    }
  }
  write.table(DF[,c(1:3)], paste(wd,"/jabba_results_",typeofrun,AAversion,"/COMPLX_REARRANG_bypatient_infocoords/", patient_id, "/", patient_id, "_", type_DF, "_",l,".txt", sep=""), sep="\t", row.names=FALSE, quote=FALSE, col.names=TRUE)
}

print(paste(patient_id, ": Completed", sep=""))
