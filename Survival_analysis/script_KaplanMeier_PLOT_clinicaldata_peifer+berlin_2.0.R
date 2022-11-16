####### KAPLAN MEIER PLOT GENERATOR #######

#libraries
library(dplyr)
library(survival)
library(survminer)

#arguments
args <- commandArgs(TRUE)

cr_type <- args[1] #"chromothripsis"
run_type <- "pass"

#open clinical data
cli_data <- read.csv(file=paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/peifer+berlin_clinicaldata_excludedpatientsKM_",cr_type,".csv",sep=""), header=TRUE, sep=",")

#ALL PATIENTS
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_data$overall_survival, cli_data$status_DOA)
#function to produce the Kaplan-Meier estimates of the probability of survival over time PER RISK GROUP
fit_probsurvival_RG <- survfit(survival_object ~ cli_data[[run_type]], data = cli_data)
pvalue_ALL <- surv_pvalue(fit_probsurvival_RG)$pval
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/KM_plot_overall_survival_ALLpatients_",cr_type,"_",run_type,".pdf",sep=""), width = 10, height = 7)
ggsurvplot(fit_probsurvival_RG, pval=TRUE, data = cli_data, xlab = "Time (d)", ylab = "Overall survival", title = paste("Overall survival ALL patients with/without ",cr_type, " (", run_type,")",sep=""))
dev.off()

#HRMNA PATIENTS
cli_data_HRMNA <- cli_data[cli_data$risk_group=="HR_MNA",]
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_data_HRMNA$overall_survival, cli_data_HRMNA$status_DOA)
#function to produce the Kaplan-Meier estimates of the probability of survival over time PER RISK GROUP
fit_probsurvival_RG <- survfit(survival_object ~ cli_data_HRMNA[[run_type]], data = cli_data_HRMNA)
pvalue_HRMNA <- surv_pvalue(fit_probsurvival_RG)$pval
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/KM_plot_overall_survival_HRMNApatients_",cr_type,"_",run_type,".pdf",sep=""), width = 10, height = 7)
ggsurvplot(fit_probsurvival_RG, pval=TRUE, data = cli_data_HRMNA, xlab = "Time (d)", ylab = "Overall survival", title = paste("Overall survival HRMNA patients with/without ",cr_type, " (", run_type,")",sep=""))
dev.off()

#HRnonMNA PATIENTS
cli_data_HRnonMNA <- cli_data[cli_data$risk_group=="HR_non_MNA",]
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_data_HRnonMNA$overall_survival, cli_data_HRnonMNA$status_DOA)
#function to produce the Kaplan-Meier estimates of the probability of survival over time PER RISK GROUP
fit_probsurvival_RG <- survfit(survival_object ~ cli_data_HRnonMNA[[run_type]], data = cli_data_HRnonMNA)
pvalue_HRnonMNA <- surv_pvalue(fit_probsurvival_RG)$pval
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/KM_plot_overall_survival_HRnonMNApatients_",cr_type,"_",run_type,".pdf",sep=""), width = 10, height = 7)
ggsurvplot(fit_probsurvival_RG, pval=TRUE, data = cli_data_HRnonMNA, xlab = "Time (d)", ylab = "Overall survival", title = paste("Overall survival HRnonMNA patients with/without ",cr_type, " (", run_type,")",sep=""))
dev.off()

#nonHR PATIENTS
cli_data_nonHR <- cli_data[cli_data$risk_group=="non_HR",]
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_data_nonHR$overall_survival, cli_data_nonHR$status_DOA)
#function to produce the Kaplan-Meier estimates of the probability of survival over time PER RISK GROUP
fit_probsurvival_RG <- survfit(survival_object ~ cli_data_nonHR[[run_type]], data = cli_data_nonHR)
pvalue_nonHR <- surv_pvalue(fit_probsurvival_RG)$pval
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/KM_plot_overall_survival_nonHRpatients_",cr_type,"_",run_type,".pdf",sep=""), width = 10, height = 7)
ggsurvplot(fit_probsurvival_RG, pval=TRUE, data = cli_data_nonHR, xlab = "Time (d)", ylab = "Overall survival", title = paste("Overall survival nonHR patients with/without ",cr_type, " (", run_type,")",sep=""))
dev.off()

#save all the info in one line
DF_line <- data.frame(cr_type,pvalue_ALL,pvalue_HRMNA,pvalue_HRnonMNA,pvalue_nonHR)
write.table(DF_line, paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/pvalues_complexrearrang/KMsurvival_info_pvalues_",cr_type,".txt",sep=""),sep="\t",row.names=FALSE,quote=FALSE)



