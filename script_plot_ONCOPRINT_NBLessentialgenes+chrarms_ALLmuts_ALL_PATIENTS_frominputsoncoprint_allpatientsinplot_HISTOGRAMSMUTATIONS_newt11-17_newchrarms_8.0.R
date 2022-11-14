########################################################
#### SCRIPT TO PLOT ONCOPRINT FINAL WITH HISTOGRAMS ####
########################################################

#libraries
library(ComplexHeatmap)
library(ggplot2)
library(reshape)


#######################
#### SVS HISTOGRAM ####
#######################

############# DECLARE LISTS #############

MUT_x_PATIENT_plot = NULL
MUT_x_all_PATIENTS_plot = NULL

############# OPEN AND PREPARE MUTATION FILES #############

#open file with patient_id, risk group and cohort
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
#Exclude patients due to abnormally high number of mutations
patients_type_file <- patients_type_file[patients_type_file$V1!="CB2044",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL47",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL53",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL54",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL49",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL50",]

#open results variant callers
#berlin+peifer
x_all<- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_ALLSVS_intraSVs+translocations_merge_svaba_delly2_novobreak_filter_ATLEAST2CALLERS_newinfostrands_precisioninfo_21.03.17.txt", sep="")

#Create factor with patient id
patients_list <- patients_type_file[,1]


############# GENERATE INPUT FOR PLOTS #############

###  SVS
#open each file and convert to format for num of mutations
for (i in 1:length(patients_list)){
    patient_id <- patients_list[i]
    cohort <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),2] #cojo la columna de type of patient
    type_patient <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),3]
    res_new2 <- x_all[x_all$V1==patient_id,]
    num_muts <- dim(res_new2)[1]  #numero total de mutaciones
    num_dels <- dim(res_new2[res_new2$V7=="deletion",])[1]
    num_dups <- dim(res_new2[res_new2$V7=="tandem-duplication",])[1]
    num_inv <- dim(res_new2[res_new2$V7=="inversion",])[1]
    num_ins <- dim(res_new2[res_new2$V7=="insertion",])[1]
    num_trans <- dim(res_new2[res_new2$V7=="translocation",])[1]
    type_mut <- "SVs"
    MUT_x_PATIENT_plot <- data.frame(patient_id,cohort, type_patient,num_muts,type_mut,num_dels,num_dups,num_inv,num_ins,num_trans)
    #data frame final
    MUT_x_all_PATIENTS_plot <- rbind(MUT_x_all_PATIENTS_plot, MUT_x_PATIENT_plot) #DATA FRAME TO PLOT NUM. OF MUTATIONS PER PATIENT
}

row.names(MUT_x_all_PATIENTS_plot) <- MUT_x_all_PATIENTS_plot[,1]
MUT_x_all_PATIENTS_plot_ONCOPRINT <- MUT_x_all_PATIENTS_plot[,c(6:10)]
MUT_x_all_PATIENTS_plot_ONCOPRINT_sv <- data.matrix(MUT_x_all_PATIENTS_plot_ONCOPRINT)


########################
#### CNVS HISTOGRAM ####
########################

############# DECLARE LISTS #############

MUT_x_PATIENT_plot = NULL
MUT_x_all_PATIENTS_plot = NULL


############# OPEN AND PREPARE MUTATION FILES #############

#open file with patient_id, risk group and cohort
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
#Exclude patients due to abnormally high number of mutations
patients_type_file <- patients_type_file[patients_type_file$V1!="CB2044",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL47",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL53",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL54",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL49",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL50",]

#Exclude patients without results from ASCAT
#patients_type_file <- patients_type_file[patients_type_file$V1!="NBL31",]
#patients_type_file <- patients_type_file[patients_type_file$V1!="NBL36",]
#patients_type_file <- patients_type_file[patients_type_file$V1!="NBL61",]

#open results variant callers
#berlin+peifer
x_all<- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_CNV_GAINS+LOSSES+AMP+HOMDELS_ASCAT_martin_results_cutoffploidy_ALLpatients_21.04.13.txt", sep="")

#Create factor with patient id
patients_list <- patients_type_file[,1]

############# GENERATE INPUT FOR PLOTS #############

###  CNVS
#open each file and convert to format for num of mutations
for (i in 1:length(patients_list)){
    patient_id <- patients_list[i]
    cohort <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),2] #cojo la columna de type of patient
    type_patient <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),3]
    res_new2 <- x_all[x_all$V1==patient_id,]
    num_muts <- dim(res_new2)[1]  #numero total de mutaciones
    num_gain <- dim(res_new2[res_new2$V7=="gain",])[1]
    num_loss <- dim(res_new2[res_new2$V7=="loss",])[1]
    num_amps <- dim(res_new2[res_new2$V7=="amplification",])[1]
    num_homdels <- dim(res_new2[res_new2$V7=="homozygous_loss",])[1]
    type_mut <- "CNVs"
    MUT_x_PATIENT_plot <- data.frame(patient_id,cohort,type_patient,num_muts,type_mut,num_gain,num_loss,num_amps,num_homdels)
    #data frame final
    MUT_x_all_PATIENTS_plot <- rbind(MUT_x_all_PATIENTS_plot, MUT_x_PATIENT_plot) #DATA FRAME TO PLOT NUM. OF MUTATIONS PER PATIENT
}

row.names(MUT_x_all_PATIENTS_plot) <- MUT_x_all_PATIENTS_plot[,1]
MUT_x_all_PATIENTS_plot_ONCOPRINT <- MUT_x_all_PATIENTS_plot[,c(6:9)]
MUT_x_all_PATIENTS_plot_ONCOPRINT_cnv <- data.matrix(MUT_x_all_PATIENTS_plot_ONCOPRINT)


########################
#### SNVS HISTOGRAM ####
########################

############# DECLARE LISTS #############

MUT_x_PATIENT_plot = NULL
MUT_x_all_PATIENTS_plot = NULL

############# OPEN AND PREPARE MUTATION FILES #############

#open file with patient_id, risk group and cohort
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
#Exclude patients due to abnormally high number of mutations
patients_type_file <- patients_type_file[patients_type_file$V1!="CB2044",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL47",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL53",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL54",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL49",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL50",]

#open results variant callers
#berlin+peifer
x_all_indels<- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_INDELS_mutect2_PASS_newinfostrands_21.01.26.txt", sep="")
x_all_snv<- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/results_21.01.26/RESULTS_all_patients_berlin+peifer_cohorts_SNVS_mutect_PASS_21.01.26.txt", sep="")

#transform type of mut. deletion/insertion to indel
x_all_indels$V7 <- "indel"

#grup lists indels+snv
x_all <- rbind(x_all_snv,x_all_indels)

#Create factor with patient id
patients_list <- patients_type_file[,1]


############# GENERATE INPUT FOR PLOTS #############

###SNVS + INDELS
for (i in 1:length(patients_list)){
    patient_id <- patients_list[i]
    cohort <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),2] #extract info of cohort
    type_patient <- patients_type_file[patients_type_file$V1==grep(patient_id, x = patients_type_file$V1,value=TRUE),3] #extract info of risk group
    res_new2 <- x_all[x_all$V1==patient_id,]
    num_muts <- dim(res_new2)[1]  #total num. of mutations
    num_snvs <- dim(res_new2[res_new2$V7=="snv",])[1] #total num. of snvs
    num_indels <- dim(res_new2[res_new2$V7=="indel",])[1] #total num. of indels
    type_mut <- "snvs+indels"
    MUT_x_PATIENT_plot <- data.frame(patient_id,cohort, type_patient,num_muts,num_snvs,num_indels)
    #data frame final
    MUT_x_all_PATIENTS_plot <- rbind(MUT_x_all_PATIENTS_plot,MUT_x_PATIENT_plot) #DATA FRAME TO PLOT NUM. OF MUTATIONS PER PATIENT
}


row.names(MUT_x_all_PATIENTS_plot) <- MUT_x_all_PATIENTS_plot[,1]
MUT_x_all_PATIENTS_plot_ONCOPRINT <- MUT_x_all_PATIENTS_plot[,c(5:6)]
MUT_x_all_PATIENTS_plot_ONCOPRINT_snv <- data.matrix(MUT_x_all_PATIENTS_plot_ONCOPRINT)



###################
#### ONCOPRINT ####
###################


############# DECLARE LISTS #############

FINAL_LIST_all_muts = NULL
PATIENTS_clinic = NULL
FINAL_LIST_GENES = NULL
FINAL_LIST_CHRARMS_G = NULL
FINAL_LIST_CHRARMS_L = NULL
FINAL_LIST_DUMMY = NULL
FINAL_LIST_TRANS = NULL

############# OPEN AND PREPARE INPUT FILES #############

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v3.0.txt",sep="")
#Change format of age
for (a in 1:dim(patients_type_file)[1]){
   if (patients_type_file[a,]$V7<365){
    patients_type_file[a,8] <- "<1year"
   } else if (patients_type_file[a,]$V7>=365 & patients_type_file[a,]$V7<=(5*365)){
    patients_type_file[a,8] <- "1-5years"
   } else if (patients_type_file[a,]$V7>(5*365)){
    patients_type_file[a,8] <- ">5years"
   }
}

#Format
patients_type_file <- patients_type_file[,c(1,2,3,5,6,8)]
names(patients_type_file) <- c("sample_id", "cohort","risk_group", "sex", "stage", "age")

#Exclude patients
patients_type_file <- patients_type_file[patients_type_file$sample_id!="CB2044",]
patients_type_file <- patients_type_file[patients_type_file$sample_id!="NBL47",]
patients_type_file <- patients_type_file[patients_type_file$sample_id!="NBL53",]
patients_type_file <- patients_type_file[patients_type_file$sample_id!="NBL54",]
patients_type_file <- patients_type_file[patients_type_file$sample_id!="NBL49",]
patients_type_file <- patients_type_file[patients_type_file$sample_id!="NBL50",]


############# OPEN AND PREPARE MUTATION FILES #############

#Generate a dummy list to keep all patients
dummy <- patients_type_file[,c(1,1,1,1)]
names(dummy) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")
dummy$Gene.Symbol <- "FLAG1"
dummy$mut_type <- "dummy"
dummy$conseq <- "."

#Inputs generated with the script script_generate_oncoprint_inputs_from_intersect_genesandchromosomearms_allmuttypes.R
#SNVS genes
snv_genes <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SNVS_oncoprint_NBLgenesessential.txt", sep="")
names(snv_genes) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#SNVS oncohistones
snv_oncohist <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SNVS_oncoprint_NBLoncohistones.txt", sep="")
names(snv_oncohist) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#INDELS genes
indels_genes <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/INDELS_oncoprint_NBLgenesessential.txt", sep="")
names(indels_genes) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#INDELS oncohistones
indels_oncohist <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/INDELS_oncoprint_NBLoncohistones.txt", sep="")
names(indels_oncohist) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#SVS genes
svs_genes <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_oncoprint_NBLgenesessential.txt", sep="")
names(svs_genes) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")
svs_genes$mut_type <- "sv"

#SVS genes
svs_genes_close <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_around20kbp_oncoprint_NBLgenesessential.txt", sep="")
names(svs_genes_close) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")
svs_genes_close$mut_type <- "structv_close"

#SVS oncohistones
#svs_oncohist <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/SVS_oncoprint_NBLoncohistones.txt", sep="")
#names(svs_oncohist) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")
#svs_oncohist$mut_type <- "sv"
#svs_oncohist = NULL #no results in sv in oncohistones

#CNVS genes
cnvs_genes <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/CNVS_oncoprint_NBLgenesessential.txt", sep="")
names(cnvs_genes) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#Genomic regions
cnvs_chrarms <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/CNVS_oncoprint_NBLchromarms_v2.txt", sep="")
names(cnvs_chrarms) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")

#Translocations chrarms
trans_chrarms <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/oncoprint/inputs_oncoprint_intersect_mutations_genesandregionsofinterest_NBL/Transloc11_17_oncoprint_NBLchromarms_v2.txt", sep="")
names(trans_chrarms) <- c("Gene.Symbol", "sample_id", "mut_type", "conseq")


############# FINAL MUTATION FILES #############

#Join all tables
#all_muts <- unique(rbind(snv_genes, snv_oncohist, indels_genes, indels_oncohist, svs_genes, svs_oncohist, cnvs_genes, cnvs_chrarms))

all_muts <- unique(rbind(dummy, snv_genes, snv_oncohist, indels_genes, indels_oncohist, svs_genes, svs_genes_close, cnvs_genes))
all_muts_chrarms <- unique(rbind(cnvs_chrarms))
all_muts_trans <- unique(rbind(trans_chrarms))

#Exclude patients
#all muts in genes
all_muts <- all_muts[all_muts$sample_id!="CB2044",]
all_muts <- all_muts[all_muts$sample_id!="NBL47",]
all_muts <- all_muts[all_muts$sample_id!="NBL53",]
all_muts <- all_muts[all_muts$sample_id!="NBL54",]
all_muts <- all_muts[all_muts$sample_id!="NBL49",]
all_muts <- all_muts[all_muts$sample_id!="NBL50",]
#all muts in chr. arms
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="CB2044",]
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="NBL47",]
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="NBL53",]
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="NBL54",]
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="NBL49",]
all_muts_chrarms <- all_muts_chrarms[all_muts_chrarms$sample_id!="NBL50",]
#all muts in trans
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="CB2044",]
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="NBL47",]
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="NBL53",]
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="NBL54",]
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="NBL49",]
all_muts_trans <- all_muts_trans[all_muts_trans$sample_id!="NBL50",]

#Prepare a DF for each selected consequence
Z_dummy <- all_muts[grep("dummy", all_muts$mut_type),]
#Mut genes
Z_missense <- all_muts[grep("^missense_variant", all_muts$conseq),]
Z_stop_gained <- all_muts[grep("^stop_gained", all_muts$conseq),]
Z_svs <- all_muts[grep("sv", all_muts$mut_type),]
Z_svs_close <- all_muts[grep("close", all_muts$mut_type),]
Z_amplification <- all_muts[grep("amplification", all_muts$mut_type),]
Z_homdels <- all_muts[grep("homozygous_loss", all_muts$mut_type),]
#Mut chr arms
Z_gainchr <- all_muts_chrarms[grep("gain", all_muts_chrarms$mut_type),]
Z_losschr <- all_muts_chrarms[grep("loss", all_muts_chrarms$mut_type),]
#Mut trans
Z_trans <- all_muts_trans[grep("translocation", all_muts_trans$mut_type),]

#Correct redundancy between SVS within gene and SVS close to gene
Z_svs$flag <- paste(Z_svs[,1],"_",Z_svs[,2], sep="")
Z_svs_close$flag <- paste(Z_svs_close[,1],"_",Z_svs_close[,2], sep="")

for (i in 1:dim(Z_svs_close)[1]){
    if (dim(Z_svs[Z_svs$flag==Z_svs_close$flag[i],])[1]>0){
       Z_svs_close <- Z_svs_close[-i,]
    }
}

#list DFs type of mutations (Z_snv, Z_indels)
tables_mutations <- objects(pattern="Z_")


############# CREATE MATRIX FOR EACH MUTATION TYPE #############

#List of matrix of mutated genes; 1 matrix per mut file
for (i in 1:length(tables_mutations)){
    DF_muts <- get(tables_mutations[i])
    table_muts <- table(DF_muts[,c("Gene.Symbol","sample_id")])
    table_muts[table_muts>0] <- 1 #convert count frequencies to 1 (1 present; 0 absent)
    df_muts <- as.data.frame.matrix(table_muts) #convert table to data frame
    FINAL_matrix_muts <- data.matrix(df_muts)
    FINAL_LIST_all_muts[[i]] <- FINAL_matrix_muts
}

#Unify all lists: lists with same patients and genes
FINAL_LIST_all_muts_unified <- unify_mat_list(FINAL_LIST_all_muts, default = 0) #convert freq=0 if not in list
names(FINAL_LIST_all_muts_unified) <- factor(tables_mutations)


#Make sure patients_type_file has the same patients that matrix (FINAL_LIST_all_muts)
patients_frommatrix <- sort(colnames(FINAL_LIST_all_muts_unified[[1]]))
for (p in 1:length(patients_frommatrix)){
    patients_line <- patients_type_file[patients_type_file$sample_id==patients_frommatrix[p],]
    PATIENTS_clinic <- rbind(PATIENTS_clinic,patients_line)
}

#Sort by risk group
PATIENTS_clinic_sort <- PATIENTS_clinic[order(PATIENTS_clinic[,3]),]

#Sort list of matrix the same way
for (l in 1:length(FINAL_LIST_all_muts_unified)){
    FINAL_LIST_all_muts_unified[[l]] <- FINAL_LIST_all_muts_unified[[l]][,PATIENTS_clinic_sort$sample_id]
}

names(FINAL_LIST_all_muts_unified) <- factor(tables_mutations)


########## SPLIT MATRIX ##########

#remove the dummy list
FINAL_LIST_all_muts_unified <- FINAL_LIST_all_muts_unified[-2]

#Split matrix in genes and chr arms
LIST_GENES_1 <- FINAL_LIST_all_muts_unified[-c(2,4,9)]
LIST_CHRARMS_G <- FINAL_LIST_all_muts_unified[2]
LIST_CHRARMS_L <- FINAL_LIST_all_muts_unified[4]
LIST_TRANS <- FINAL_LIST_all_muts_unified[9]

for (m in 1:length(LIST_GENES_1)){
    FINAL_LIST_GENES[[m]] <- LIST_GENES_1[[m]][grep("_|;",rownames(LIST_GENES_1[[m]]), invert=TRUE),]
}

for (m in 1:length(LIST_CHRARMS_G)){
    FINAL_LIST_CHRARMS_G[[m]] <- LIST_CHRARMS_G[[m]][grep("_gain",rownames(LIST_CHRARMS_G[[m]])),]
}

for (m in 1:length(LIST_CHRARMS_L)){
    FINAL_LIST_CHRARMS_L[[m]] <- LIST_CHRARMS_L[[m]][grep("_loss",rownames(LIST_CHRARMS_L[[m]])),]
}

for (m in 1:length(LIST_TRANS)){
    FINAL_LIST_TRANS[[m]] <- LIST_TRANS[[m]][grep("17",rownames(LIST_TRANS[[m]])),]
}

names(FINAL_LIST_GENES) <- names(LIST_GENES_1)
names(FINAL_LIST_CHRARMS_G) <- names(LIST_CHRARMS_G)
names(FINAL_LIST_CHRARMS_L) <- names(LIST_CHRARMS_L)
names(FINAL_LIST_TRANS) <- names(LIST_TRANS)


############# RETRIEVE CLINICAL ANNOTATION #############

risk_group <- PATIENTS_clinic_sort$risk_group
sex <- PATIENTS_clinic_sort$sex
stage <- PATIENTS_clinic_sort$stage
age <- PATIENTS_clinic_sort$age

#Reorder histograms for anotation
MUT_x_all_PATIENTS_plot_ONCOPRINT_snv_ord <- MUT_x_all_PATIENTS_plot_ONCOPRINT_snv[order(match(row.names(MUT_x_all_PATIENTS_plot_ONCOPRINT_snv),colnames(FINAL_LIST_all_muts_unified[[1]]))),]
MUT_x_all_PATIENTS_plot_ONCOPRINT_sv_ord <- MUT_x_all_PATIENTS_plot_ONCOPRINT_sv[order(match(row.names(MUT_x_all_PATIENTS_plot_ONCOPRINT_sv),colnames(FINAL_LIST_all_muts_unified[[1]]))),]
MUT_x_all_PATIENTS_plot_ONCOPRINT_cnv_ord <- MUT_x_all_PATIENTS_plot_ONCOPRINT_cnv[order(match(row.names(MUT_x_all_PATIENTS_plot_ONCOPRINT_cnv),colnames(FINAL_LIST_all_muts_unified[[1]]))),]

#AnnotatiON
annotation <- HeatmapAnnotation(snv=anno_barplot(MUT_x_all_PATIENTS_plot_ONCOPRINT_snv_ord, border=FALSE, bar_width=1, gp=gpar(fill=c("#FFCC99", "#663434"), col=0)), svs=anno_barplot(MUT_x_all_PATIENTS_plot_ONCOPRINT_sv_ord, border=FALSE, bar_width=1, gp=gpar(fill=c("#FFCC99", "#6699CC", "#CC9933", "#CCCC99", "#333366"), col=0)), cnvs=anno_barplot(MUT_x_all_PATIENTS_plot_ONCOPRINT_cnv_ord, border=FALSE, bar_width=1, gp=gpar(fill=c("#99CCFF", "#FDE496", "#CC3333", "#B3B3B3"), col=0)), risk_group = risk_group, sex = sex, stage = stage, age = age, col = list(risk_group = c("HR_MNA" = "#CF8484", "HR_non_MNA" = "#976E97", "non_HR" = "#9CB4CB"), sex = c("F"="#CAB097", "M"="#B9E0F3", "na" = "#EBEBEB"), stage = c("1"="#F2F2F2", "2"="#E1E1CF", "2A"="#E1E1BA", "2B"="#E1E19D", "3"="#898989", "4"="#FAD2FF", "4S"="#FFD2D2", "5"= "#333333" ), age=c("<1year"="#DBBF7F", "1-5years"="#AB8271", ">5years"="#61AB9A")), annotation_height = unit(c(30, 30, 30, 3, 3, 3, 3), "mm"), annotation_legend_param = list(risk_group = list(title = "Risk Group"), sex=list(title="Sex"), stage=list(title="Stage"), age=list (title = "Age")), gap=unit(7,"points"))


############# ONCOPRINT PLOT #############

#Colors per type of mutation
col = c(Z_missense="#E6C541", Z_stop_gained="#E6A541", Z_svs="#7960B4",  Z_svs_close="#609EB4", Z_amplification="#FFAAAA", Z_homdels="#CCCCCC", Z_gainchr="#5396D8", Z_losschr="#E09A55", Z_trans="#73C8AE")

#Dimension and colors for each mutation
alter_fun_GENES = list(
    background = alter_graphic("rect", fill = "#FFFEEE"),
    Z_homdels = alter_graphic("rect", fill = col["Z_homdels"]),
    Z_amplification = alter_graphic("rect", fill = col["Z_amplification"]),
    Z_svs_close = alter_graphic("rect", height = 0.50, fill = col["Z_svs_close"]),
    Z_svs = alter_graphic("rect", height = 0.50, fill = col["Z_svs"]),
    Z_missense = alter_graphic("rect", height = 0.20, fill = col["Z_missense"]),
    Z_stop_gained = alter_graphic("rect", height = 0.20, fill = col["Z_stop_gained"])
)

#Dimension and colors for each mutation
alter_fun_CHRARMS = list(
    background = alter_graphic("rect", fill = "#FFFFFF"),
    Z_gainchr = alter_graphic("rect", fill = col["Z_gainchr"]),
    Z_losschr = alter_graphic("rect", fill = col["Z_losschr"])
)

#Dimension and colors for each mutation
alter_fun_TRANS = list(
    background = alter_graphic("rect", fill = "#F0FFF4"),
    Z_trans = alter_graphic("rect", fill = col["Z_trans"])
)

column_title = "OncoPrint Mutated genes NBL (Peifer+Berlin)"
oncoprint_legend_param_genes = list(title = "Mutation type", at = c("Z_missense","Z_stop_gained", "Z_svs_close", "Z_svs", "Z_amplification", "Z_homdels"), labels = c("Missense","Stop gained", "SV close", "SV", "amplification", "homozygous_loss"))

oncoprint_legend_param_chrarms = list(title = "Mutation type", at = c("Z_gainchr", "Z_losschr"), labels = c("CN gain", "CN loss"))

oncoprint_legend_param_trans = list(title = "Mutation type", at = c("Z_trans"), labels = c("t(11;17)"))


#Plot (I split the plot per risk group)
PLOT1 <- oncoPrint(FINAL_LIST_CHRARMS_G, alter_fun = alter_fun_CHRARMS, col = col, column_title = column_title, heatmap_legend_param = oncoprint_legend_param_chrarms, column_km = 3, column_split = risk_group, top_annotation = annotation, remove_empty_rows = TRUE, show_column_names = FALSE, remove_empty_columns = FALSE, pct_side = "right", row_names_side = "left", right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(ylim =c(0,35)))) %v%
    oncoPrint(FINAL_LIST_CHRARMS_L, alter_fun = alter_fun_CHRARMS, col = col, column_title = column_title, heatmap_legend_param = oncoprint_legend_param_chrarms, column_km = 3, column_split = risk_group, top_annotation = NULL, remove_empty_rows = TRUE, show_column_names = FALSE, remove_empty_columns = FALSE, pct_side = "right", row_names_side = "left", right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(ylim =c(0,35)))) %v%
    oncoPrint(FINAL_LIST_TRANS, alter_fun = alter_fun_TRANS, col = col, column_title = column_title, heatmap_legend_param = oncoprint_legend_param_trans, column_km = 3, column_split = risk_group, top_annotation = NULL, remove_empty_rows = TRUE, show_column_names = FALSE, remove_empty_columns = FALSE, pct_side = "right", row_names_side = "left", right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(ylim =c(0,35)))) %v%
    oncoPrint(FINAL_LIST_GENES, alter_fun = alter_fun_GENES, col = col, column_title = column_title, heatmap_legend_param = oncoprint_legend_param_genes, top_annotation = NULL, remove_empty_rows = TRUE, show_column_names = FALSE, remove_empty_columns = FALSE, pct_side = "right", row_names_side = "left", right_annotation = rowAnnotation(row_barplot = anno_oncoprint_barplot(ylim =c(0,35))))


pdf(file="ONCOPRINT_NBLessentialgenes+chrarms+trans+HISTOGRAMS_ALLmuts_ALL_PATIENTS_PEIFER+BERLIN_filtered_PASS_newt11-17_newchrarms_21.12.8.pdf", width = 12, height = 11)
PLOT1
dev.off()
