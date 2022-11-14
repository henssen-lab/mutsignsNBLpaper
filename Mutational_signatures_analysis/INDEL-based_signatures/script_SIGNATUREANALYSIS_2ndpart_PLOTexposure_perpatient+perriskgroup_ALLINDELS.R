##############################################################
###                   SIGNATURE ANALYSIS                   ###
###    2ND SCRIPT: PLOT MUTATIONAL SIGNATURE EXPOSURES     ###
##############################################################


#libraries
library(stringr)
library(reshape2)
library(ggplot2)
library(wesanderson)


###### OPEN INPUT FILES ######

###Open exposure tables for COUNTS and FREQ exposures
exp_counts <- read.delim("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_COUNT_22.06.13.txt", sep="\t")
exp_freq <- read.delim("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt", sep="\t")


###### EXPOSURE PLOT COSMIC SIGNATURES PER PATIENT ######

#Create input for exposure plot ABSOLUTE
input_plot_exp2 <- melt(exp_counts)
signs_names <- names(exp_freq[c(-1,-length(names(exp_freq)))])
input_plot_exp2$variable <- ordered(input_plot_exp2$variable, levels=c(str_sort(signs_names, numeric=TRUE)))

# Plot exposures
PLOT1 <- ggplot(input_plot_exp2, aes(x=reorder(sample_id, -value),y=value, fill=variable))
PLOT2 <- PLOT1 + geom_col(position = position_stack(reverse = TRUE)) + facet_grid(~risk_group, scales = "free_x", space = "free_x") +scale_y_continuous(name="Count") + scale_x_discrete(name="Sample ID") + labs(title = "Signature's exposure" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Signature")) + theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/plots/Plot_ALLINDELS_COSMICsignatures_absexposure_COUNTSPERPATIENT_all114patients_22.06.13.pdf", width = 20, height = 5, plot=PLOT2)

#Create input for exposure plot FREQUENCY
input_plot_exp2 <- melt(exp_freq)
input_plot_exp2$variable <- ordered(input_plot_exp2$variable, levels=c(str_sort(signs_names, numeric=TRUE)))

# Plot exposures
PLOT1 <- ggplot(input_plot_exp2, aes(x=reorder(sample_id, -value),y=value, fill=variable))
PLOT2 <- PLOT1 + geom_col(position = position_fill(reverse = TRUE)) + facet_grid(~risk_group, scales = "free_x", space = "free_x") +scale_y_continuous(name="Freq") + scale_x_discrete(name="Sample ID") + labs(title = "Signature's exposure" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Signature")) + theme_classic() + theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())
ggsave(filename="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/plots/Plot_ALLINDELS_COSMICsignatures_absexposure_FREQPERPATIENT_all114patients_22.06.13.pdf", width = 20, height = 5, plot=PLOT2)


###### EXPOSURE PLOT COSMIC SIGNATURES PER RISK GROUP - COUNTS######

#Create input for exposure plot by risk group
total.exp.perrisk <- exp_counts[1:3,c(1:(dim(exp_counts)[2]-1))] #I generate the data frame rows: risk groups, cols: signatures
names(total.exp.perrisk)[1] <- "risk_group"
total.exp.perrisk[1,1] <- "HR_MNA"
total.exp.perrisk[2,1] <- "HR_non_MNA"
total.exp.perrisk[3,1] <- "non_HR"

#I correct by number of samples per risk group
for (i in 2:(dim(exp_counts)[2]-1)){
  exp_MNA <- sum(exp_counts[exp_counts$risk_group=="HR_MNA",i])/dim(exp_counts[exp_counts$risk_group=="HR_MNA",])[1]
  exp_HR <- sum(exp_counts[exp_counts$risk_group=="HR_non_MNA",i])/dim(exp_counts[exp_counts$risk_group=="HR_non_MNA",])[1]
  exp_no <- sum(exp_counts[exp_counts$risk_group=="non_HR",i])/dim(exp_counts[exp_counts$risk_group=="non_HR",])[1]
  total.exp.perrisk[1,i] <- exp_MNA
  total.exp.perrisk[2,i] <- exp_HR
  total.exp.perrisk[3,i] <- exp_no
}

input_plot_exp2 <- melt(total.exp.perrisk)
input_plot_exp2$variable <- ordered(input_plot_exp2$variable, levels=c(str_sort(signs_names, numeric=TRUE)))

#plot with percentage
PLOT2 <- ggplot(input_plot_exp2, aes(x=risk_group,y=value, fill=variable))
PLOT3 <- PLOT2 + geom_col(position = position_stack(reverse = TRUE)) +scale_y_continuous(name="# of INDELS caused by the signature / num of samples") + scale_x_discrete(name="Risk Group") + labs(title = "Counts of COSMIC signatures by risk group (corrected by number of samples)" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Signature")) + theme_classic()
ggsave(filename="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/plots/Plot_ALLINDELS_absexp_INDELS_COSMIC_signatures_COUNTS_RISKGROUP_22.06.13.pdf", width = 10, height = 7, plot=PLOT3)



###### EXPOSURE PLOT COSMIC SIGNATURES PER RISK GROUP - FREQUENCY ######


#Create input for exposure plot by risk group
total.exp.perrisk <- exp_freq[1:3,c(1:(dim(exp_freq)[2]-1))] #I generate the data frame rows: risk groups, cols: signatures
names(total.exp.perrisk)[1] <- "risk_group"
total.exp.perrisk[1,1] <- "HR_MNA"
total.exp.perrisk[2,1] <- "HR_non_MNA"
total.exp.perrisk[3,1] <- "non_HR"

for (i in 2:(dim(exp_freq)[2]-1)){
  exp_MNA <- sum(exp_freq[exp_freq$risk_group=="HR_MNA",i])
  exp_HR <- sum(exp_freq[exp_freq$risk_group=="HR_non_MNA",i])
  exp_no <- sum(exp_freq[exp_freq$risk_group=="non_HR",i])
  total.exp.perrisk[1,i] <- exp_MNA*100/dim(exp_freq[exp_freq$risk_group=="HR_MNA",])[1]
  total.exp.perrisk[2,i] <- exp_HR*100/dim(exp_freq[exp_freq$risk_group=="HR_non_MNA",])[1]
  total.exp.perrisk[3,i] <- exp_no*100/dim(exp_freq[exp_freq$risk_group=="non_HR",])[1]
}

#Create input for plot
input_plot_exp2 <- melt(total.exp.perrisk)
input_plot_exp2$variable <- ordered(input_plot_exp2$variable, levels=c(str_sort(signs_names, numeric=TRUE)))

#plot with percentage
PLOT2 <- ggplot(input_plot_exp2, aes(x=risk_group,y=value, fill=variable))
PLOT3 <- PLOT2 + geom_col(position = position_fill(reverse = TRUE)) + scale_y_continuous(name="% of INDELS caused by the signature") + scale_x_discrete(name="Risk Group") + labs(title = "Percentage of COSMIC signatures by risk group" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Signature")) + theme_classic()
ggsave(filename="/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/plots/Plot_ALLINDELS_absexp_INDELS_COSMIC_signatures_FREQ_RISKGROUP_22.06.13.pdf", width = 10, height = 7, plot=PLOT3)
