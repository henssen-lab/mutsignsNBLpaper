#####################################################
#### SCRIPT TO PLOT VIOLIN PLOTS PER TYPE OF CNV ####
#####################################################


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

#Exclude patients without results from ASCAT
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL31",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL36",]
patients_type_file <- patients_type_file[patients_type_file$V1!="NBL61",]

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


############# PLOTS PER TYPE OF CNV #############

###PLOT BOXPLOTS RISK GROUPS

#gains
PLOT_CNVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_gain, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of CN gains in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_CNgains_distribution_ALL_PATIENTS_PEIFER+BERLIN_ALLCNVS.21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_CNVS)

#losses
PLOT_CNVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_loss, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of CN losses in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_CNlosses_distribution_ALL_PATIENTS_PEIFER+BERLIN_ALLCNVS.21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_CNVS)

#amplifications
PLOT_CNVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_amps, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of amplifications in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_CNamplifications_distribution_ALL_PATIENTS_PEIFER+BERLIN_ALLCNVS.21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_CNVS)

#homozygous deletions
PLOT_CNVS <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient,y = num_homdels, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distribution num. of hom. dels in the different NBL risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_CNhomdels_distribution_ALL_PATIENTS_PEIFER+BERLIN_ALLCNVS.21.04.13_risk_groups.pdf", width = 10, height = 10, plot=PLOT_CNVS)


