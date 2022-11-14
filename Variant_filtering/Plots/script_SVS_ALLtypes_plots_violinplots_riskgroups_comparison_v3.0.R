#####################################################################
#### SCRIPT TO PLOT VIOLIN PLOTS NUMBER OF SVS PER MUTATION TYPE ####
#####################################################################


#libraries
library(ggplot2)
library(wesanderson)
library(reshape)
library(dplyr)
library(ggstatsplot)


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


############# PLOTS PER SV TYPE #############

###PLOT BOXPLOTS RISK GROUPS

#deletions
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_dels, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of deletions in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_deletions_distribution_ALL_PATIENTS_PEIFER+BERLIN_INTRA+INTER_SVs_filtered_atleast2callers_pre.ca0.05_21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_SVS)

#inversions
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_inv, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of inversions in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_inversions_distribution_ALL_PATIENTS_PEIFER+BERLIN_INTRA+INTER_SVs_filtered_atleast2callers_pre.ca0.05_21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_SVS)

#duplications
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_dups, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of duplications in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_duplications_distribution_ALL_PATIENTS_PEIFER+BERLIN_INTRA+INTER_SVs_filtered_atleast2callers_pre.ca0.05_21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_SVS)

#translocations
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_trans, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of translocations in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_translocations_distribution_ALL_PATIENTS_PEIFER+BERLIN_INTRA+INTER_SVs_filtered_atleast2callers_pre.ca0.05_21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_SVS)

#insertions
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_ins, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of insertions in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_insertions_distribution_ALL_PATIENTS_PEIFER+BERLIN_INTRA+INTER_SVs_filtered_atleast2callers_pre.ca0.05_21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_SVS)
 
