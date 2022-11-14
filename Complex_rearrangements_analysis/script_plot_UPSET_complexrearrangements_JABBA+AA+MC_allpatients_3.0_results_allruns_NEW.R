###########################################################
#### SCRIPT TO PLOT UPSET PLOTS COMPLEX REARRANGEMENTS ####
###########################################################

#libraries
library(ComplexHeatmap)
library(ComplexUpset)
library(ggplot2)

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

all_muts <- unique(complex_svs)

#Exclude patients
all_muts <- all_muts[!grepl(all_muts$patient_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54"),]


############# CREATE MATRIX #############

table_muts <- table(all_muts[,c("patient_id","type")])
table_muts[table_muts>0] <- 1 #convert count frequencies to 1 (1 present; 0 absent)
dm_muts <- as.data.frame.matrix(table_muts) #convert table to data frame
df_muts <- data.frame(dm_muts)
df_muts$sample_id <- rownames(df_muts)
rownames(df_muts) <- NULL
df_muts <- df_muts[,c(dim(df_muts)[2],1:dim(df_muts)[2]-1)]
FINAL_matrix_muts <- merge(patients_type_file,df_muts,by="sample_id")


############# PLOT #############

intersect_cr <- rev(colnames(FINAL_matrix_muts[3:dim(FINAL_matrix_muts)[2]]))

PLOT1 <- upset(FINAL_matrix_muts,intersect_cr, base_annotations=list('Intersection size'=intersection_size(counts=FALSE,mapping=aes(fill=risk_group))+scale_y_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 2))+theme_classic()+ theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())),width_ratio=0.1,wrap=TRUE, set_sizes=FALSE,sort_sets=FALSE)

pdf(file=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/UPSETplot_jabba+AA+MCmergecomplexSVs_noqrp_nolinearamps_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 15)
PLOT1
dev.off()
