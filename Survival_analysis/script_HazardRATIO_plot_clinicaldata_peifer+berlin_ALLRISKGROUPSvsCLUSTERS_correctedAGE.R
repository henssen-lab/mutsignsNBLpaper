###########################################
####### HAZARD RATIO PLOT GENERATOR #######
###########################################


#libraries
library(dplyr)
library(survival)
library(survminer)

#arguments
args <- commandArgs(TRUE)

cr_type <- args[1] #"scenariosNEW_COMPHEATMAPK3_allinfo"
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
cli_data2$risk_group <- cli_data$pass

cli_dataFINAL <- rbind(cli_data,cli_data2)


### ALL RISKGROUPS vs ALL CLUSTERS

#survival object with overall survival data and status alive=0 (alive in remission + alive in relapse) and dead=1
survival_object <- Surv(cli_dataFINAL$overall_survival, cli_dataFINAL$status_DOA)
#function to produce the HAZARD RATIO plot
fit.coxph <- coxph(survival_object ~ risk_group + age  , data = cli_dataFINAL)
p_val <- summary(fit.coxph)$coefficients[,5]
#Plot
pdf(file = paste("/Users/Elias/Documents/Landscape_neuroblastoma_project.CHARITE/Clinical_data_analysis/Plots_allcomplexrearrang/HAZARDRATIO_plot_ALLpatients_compareALLRISKGROUPSvsCLUSTERS_CORRECTBYAGE.pdf",sep=""), width = 10, height = 4)
ggforest(fit.coxph, data = cli_dataFINAL)
dev.off()



