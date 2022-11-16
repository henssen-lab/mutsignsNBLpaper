####### KAPLAN MEIER PLOT GENERATOR #######

#libraries
library(dplyr)
library(survival)
library(survminer)

#arguments
args <- commandArgs(TRUE)

cr_type <- "scenariosNEW_COMPHEATMAPK3"
run_type <- "pass"

#open clinical data
cli_data <- read.csv(file=paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/peifer+berlin_clinicaldata_excludedpatientsKM_",cr_type,".csv",sep=""), header=TRUE, sep=",")
cli_data <- cli_data[,c(1:7)]

cli_data$risk_group <- gsub('non_HR', 'A_non_HR', cli_data$risk_group)
cli_data$pass <- gsub('Cluster1', 'B_non_HR_clust1', cli_data$pass)
cli_data$risk_group <- gsub('HR_non_MNA', 'C_HR_non_MNA', cli_data$risk_group)
cli_data$pass <- gsub('Cluster2', 'D_HR_non_MNA_clust2', cli_data$pass)
cli_data$risk_group <- gsub('HR_MNA', 'E_HR_MNA', cli_data$risk_group)
cli_data$pass <- gsub('Cluster3', 'F_HR_MNA_clust3', cli_data$pass)

cli_data2 <- cli_data
cli_data2$risk_group <- cli_data$pass

cli_dataFINAL <- rbind(cli_data,cli_data2)


#ALL PATIENTS
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_dataFINAL$overall_survival, cli_dataFINAL$status_DOA)
#function to produce the Kaplan-Meier estimates of the probability of survival over time PER RISK GROUP
fit_probsurvival_RG <- survfit(survival_object ~ risk_group, data = cli_dataFINAL)
pvalue_ALL <- surv_pvalue(fit_probsurvival_RG)$pval
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/KM_plot_overall_survival_ALLpatients_",cr_type,"_",run_type,".COMPARISONRISKGROUPCLASSvsCLUSTERING.pdf",sep=""), width = 10, height = 7)
ggsurvplot(fit_probsurvival_RG, pval=TRUE, data = cli_dataFINAL, xlab = "Time (d)", ylab = "Overall survival", title = paste("Overall survival ALL patients with/without ",cr_type, " (", run_type,")",sep=""))
dev.off()

