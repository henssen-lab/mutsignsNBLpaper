###########################################
####### HAZARD RATIO PLOT GENERATOR #######
###########################################


#libraries
library(dplyr)
library(survival)
library(survminer)
library(foreign)

#arguments
cr_type <- "scenariosNEW_COMPHEATMAPK3_allinfoPAPERJOHANNES"
run_type <- "pass"


### OPEN INPUT FILES

#open clinical data
cli_data <- read.csv(file=paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/peifer+berlin_clinicaldata_excludedpatientsKM_",cr_type,".csv",sep=""), header=TRUE, sep=",")

### FORMAT FILES

#Create groups to compare
cli_data$risk_group <- gsub('non_HR', 'A_non_HR', cli_data$risk_group)
cli_data$pass <- gsub('Cluster1', 'B_non_HR_clust1', cli_data$pass)
cli_data$risk_group <- gsub('HR_non_MNA', 'C_HR_non_MNA', cli_data$risk_group)
cli_data$pass <- gsub('Cluster2', 'D_HR_non_MNA_clust2', cli_data$pass)
cli_data$risk_group <- gsub('HR_MNA', 'E_HR_MNA', cli_data$risk_group)
cli_data$pass <- gsub('Cluster3', 'F_HR_MNA_clust3', cli_data$pass)

cli_data2 <- cli_data
cli_dataFINAL <- cli_data2

### RISKGROUPS + ALL COVARIATES

#Response variable and cox model
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
Resp_var <- Surv(cli_dataFINAL$overall_survival, cli_dataFINAL$status_DOA)

#Backward elimination of variables (Backward elimination is possibly the best of the above methods for identifying the important variables, and it allows one to examine the full model, which is the only fit providing accurate standard errors and P-values https://doi.org/10.1038/sj.bjc.6601117)
model1 <- coxph(Resp_var~risk_group, data=cli_dataFINAL)

#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/HAZARDRATIO_plot_ALLpatients_RISKGROUPunivariate.pdf",sep=""), width = 10, height = 2)
ggforest(model1, data = cli_dataFINAL)
dev.off()

### SCENARIOS + ALL COVARIATES

#Response variable and cox model
#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
Resp_var <- Surv(cli_dataFINAL$overall_survival, cli_dataFINAL$status_DOA)
model1 <- coxph(Resp_var~pass, data=cli_dataFINAL)

#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/HAZARDRATIO_plot_ALLpatients_CLUSTERSunivariate.pdf",sep=""), width = 10, height = 2)
ggforest(model1, data = cli_dataFINAL)
dev.off()