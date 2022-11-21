#########################################################################
#### SCRIPT TO PLOT COUNTS AND FREQ PER GROUP COMPLEX REARRANGEMENTS ####
#########################################################################


#libraries
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(wesanderson)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE INPUT FILES #############

#counts
MATRIX_complex_svs_M <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_COUNTS_perpatient_114patients_RISKGROUPinfo_",typeofrun,AAversion,".txt",sep=""),sep="\t")
#freq
MATRIX_complex_svs_M_FREQ <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_FREQ_perpatient_114patients_RISKGROUPinfo_",typeofrun,AAversion,".txt", sep=""),sep="\t")


###### EXPOSURE PLOT COSMIC SIGNATURES PER RISK GROUP - COUNTS CORRECTED NUM PATIENTS ######

#Create input for exposure plot by risk group
total.exp.perrisk <- MATRIX_complex_svs_M[1:3,c(1:(dim(MATRIX_complex_svs_M)[2]-1))] #I generate the data frame rows: risk groups, cols: signatures
names(total.exp.perrisk)[1] <- "risk_group"
total.exp.perrisk[1,1] <- "HR_MNA"
total.exp.perrisk[2,1] <- "HR_non_MNA"
total.exp.perrisk[3,1] <- "non_HR"

#I correct by the number of samples in each group
for (i in 2:(dim(MATRIX_complex_svs_M)[2]-1)){
  exp_MNA <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_MNA",i])/dim(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_MNA",])[1]
  exp_HR <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_non_MNA",i])/dim(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_non_MNA",])[1]
  exp_no <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="non_HR",i])/dim(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="non_HR",])[1]
  total.exp.perrisk[1,i] <- exp_MNA
  total.exp.perrisk[2,i] <- exp_HR
  total.exp.perrisk[3,i] <- exp_no
}

input_plot_exp2 <- melt(total.exp.perrisk)

#plot with percentage
PLOT2 <- ggplot(input_plot_exp2, aes(x=risk_group,y=value, fill=variable))
PLOT3 <- PLOT2 + geom_col(position = position_stack(reverse = TRUE)) +scale_y_continuous(name="# Complex rearrangement / num of patients") + scale_x_discrete(name="Risk Group") + labs(title = "Counts of Complex Rearrangements (JABBA+AA+MC) by risk group (corrected by num of samples in each group)" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Complex Rearrangements")) + theme_classic()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_complexrearrangements_types_COUNTS_correctednumpatients114_RISKGROUP_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 7, plot=PLOT3)


###### EXPOSURE PLOT COSMIC SIGNATURES PER RISK GROUP - COUNTS  ######

#Create input for exposure plot by risk group
total.exp.perrisk <- MATRIX_complex_svs_M[1:3,c(1:(dim(MATRIX_complex_svs_M)[2]-1))] #I generate the data frame rows: risk groups, cols: signatures
names(total.exp.perrisk)[1] <- "risk_group"
total.exp.perrisk[1,1] <- "HR_MNA"
total.exp.perrisk[2,1] <- "HR_non_MNA"
total.exp.perrisk[3,1] <- "non_HR"

#I correct by the number of samples in each group
for (i in 2:(dim(MATRIX_complex_svs_M)[2]-1)){
  exp_MNA <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_MNA",i])
  exp_HR <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="HR_non_MNA",i])
  exp_no <- sum(MATRIX_complex_svs_M[MATRIX_complex_svs_M$risk_group=="non_HR",i])
  total.exp.perrisk[1,i] <- exp_MNA
  total.exp.perrisk[2,i] <- exp_HR
  total.exp.perrisk[3,i] <- exp_no
}

input_plot_exp2 <- melt(total.exp.perrisk)

#plot with percentage
PLOT2 <- ggplot(input_plot_exp2, aes(x=risk_group,y=value, fill=variable))
PLOT3 <- PLOT2 + geom_col(position = position_stack(reverse = TRUE)) +scale_y_continuous(name="# Complex rearrangement") + scale_x_discrete(name="Risk Group") + labs(title = "Counts of Complex Rearrangements (JABBA+AA+MC) by risk group (absolute)" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Complex Rearrangements")) + theme_classic()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_complexrearrangements_types_COUNTS_notcorrectedpatients114_RISKGROUP_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 7, plot=PLOT3)


###### EXPOSURE PLOT COSMIC SIGNATURES PER RISK GROUP - FREQUENCY ######


#Create input for exposure plot by risk group
total.exp.perrisk <- MATRIX_complex_svs_M_FREQ[1:3,c(1:(dim(MATRIX_complex_svs_M_FREQ)[2]-1))] #I generate the data frame rows: risk groups, cols: signatures
names(total.exp.perrisk)[1] <- "risk_group"
total.exp.perrisk[1,1] <- "HR_MNA"
total.exp.perrisk[2,1] <- "HR_non_MNA"
total.exp.perrisk[3,1] <- "non_HR"

HR_MNA_tot <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="HR_MNA",(3:length(names(MATRIX_complex_svs_M_FREQ))-1)])
HR_tot <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="HR_non_MNA",(3:length(names(MATRIX_complex_svs_M_FREQ))-1)])
no_tot <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="non_HR",(3:length(names(MATRIX_complex_svs_M_FREQ))-1)])

for (i in 2:(dim(MATRIX_complex_svs_M_FREQ)[2]-1)){
  exp_MNA <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="HR_MNA",i])
  exp_HR <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="HR_non_MNA",i])
  exp_no <- sum(MATRIX_complex_svs_M_FREQ[MATRIX_complex_svs_M_FREQ$risk_group=="non_HR",i])
  total.exp.perrisk[1,i] <- exp_MNA*100/HR_MNA_tot
  total.exp.perrisk[2,i] <- exp_HR*100/HR_tot
  total.exp.perrisk[3,i] <- exp_no*100/no_tot
}

#Create input for plot
input_plot_exp2 <- melt(total.exp.perrisk)

#plot with percentage
PLOT2 <- ggplot(input_plot_exp2, aes(x=risk_group,y=value, fill=variable))
PLOT3 <- PLOT2 + geom_col(position = position_fill(reverse = TRUE)) + scale_y_continuous(name="% Complex rearrangement") + scale_x_discrete(name="Risk Group") + labs(title = "Percentage of Complex Rearrangements (JABBA+AA+MC) by risk group" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Complex Rearrangements")) + theme_classic()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_complexrearrangements_types_patients114_FREQ_RISKGROUP_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 7, plot=PLOT3)
