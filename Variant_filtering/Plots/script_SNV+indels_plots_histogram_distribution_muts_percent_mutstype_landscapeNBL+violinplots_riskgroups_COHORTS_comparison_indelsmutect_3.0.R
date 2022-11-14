#############################################################
#### SCRIPT TO PLOT HISTOGRAMS NUMBER OF SNVS AND INDELS ####
#############################################################


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


############# PLOTS #############

###PLOT HISTOGRAM ABSOLUTE VALUES stacked from the different mutation types, split by risk group

#generate input list
input_plot_hist <- melt(MUT_x_all_PATIENTS_plot)

#Select the rows with num_muts total
input_all_muts <- input_plot_hist[input_plot_hist$variable!="num_muts",]

PLOT1 <- ggplot(input_all_muts, aes(x=reorder(patient_id, -value),y=value, fill=variable))
PLOT1_1 <- PLOT1 + geom_bar(stat="identity", position=position_stack(reverse = TRUE)) + facet_grid(~type_patient, scales = "free_x", space = "free_x") + scale_y_continuous(name="SNVs and Indels",limit=c(0,12000)) + scale_x_discrete(name="Sample ID") + labs(title = "Number of SNVs and Indels - Mutect2 PASS" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Patient type")) + theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + scale_fill_manual(values=c("#FFCC66", "#663333"))
ggsave(filename="Plot_SNVs+indels_distribution_absolutenummuts_ALL_PATIENTS_PEIFER+BERLIN_filtered_PASS_21.04.01.pdf", width = 20, height = 5, plot=PLOT1_1)


###PLOT BOXPLOTS RISK GROUPS

#all SNVs+indels
PLOT3 <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = type_patient, y = num_muts, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distrib. num. of SNVs and Indels in the diff. NBL risk groups", var.equal=FALSE, ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = type_patient, colour=type_patient), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
ggsave(filename="Plot_boxplot_SNVs+indels_distribution_ALL_PATIENTS_PEIFER+BERLIN_filtered_PASS_21.04.01_risk_groups.pdf", width = 10, height = 10, plot=PLOT3)


###PLOT BOXPLOTS COHORTS

#all SNVs+indels
PLOT4 <- ggstatsplot::ggbetweenstats(data = MUT_x_all_PATIENTS_plot, x = cohort,y = num_muts, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Distrib. num. of SNVs and Indels in the diff. NBL cohorts",var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np",  violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = cohort, colour=cohort), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#FFCC99", "#66CC99")) + ggplot2::scale_colour_manual(values=alpha(c("#FFCC99", "#66CC99"),0.2))
ggsave(filename="Plot_boxplot_SNVs+indels_distribution_ALL_PATIENTS_PEIFER+BERLIN_filtered_PASS_21.04.01_cohort.pdf", width = 10, height = 10, plot=PLOT4)

