####################################################################
#####               SCRIPT TO PLOT CORRELATION                 #####
#####    OF SBS40 AND SBS5 CLOCK LIKE SIGNATURES WITH AGE      #####
####################################################################


#libraries
library("ggpubr")

## Prepare data frame
ALL_tests = NULL

#Patients clinical data
P_data <- read.delim("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_age_os_nummutations_forcorrelation.txt")
P_data <- P_data[,c(1,3)]

#SNV signatures
SNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/exposure_SNV_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt")
SNVsign <- SNVsign[,-dim(SNVsign)[2]]


#Are the data from each of the 2 variables (x, y) follow a normal distribution
#use Shapiro test to evaluate this
#If normal, we can use pearson method, otherwise spearman.

#DF with signatures and age
DT_ALL_PATIENTS <- merge(SNVsign, P_data, by.x="sample_id",by.y="patient_id")

#Normality test
shapiro.test(DT_ALL_PATIENTS[,3])$p.value #SBS40
shapiro.test(DT_ALL_PATIENTS[,4])$p.value #SBS5
shapiro.test(DT_ALL_PATIENTS[,6])$p.value #age


### PLOT CORRELATION WITH AGE

#SBS5
pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/correlation_signatures_withclinical/results_postREVISION_22.10.04/plots/SBS5correlationwithAGE_pearson_22.10.04.pdf", width = 10, height = 7)
ggscatter(DT_ALL_PATIENTS, x = "age", y = "SBS5", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "age (days)", ylab = "Rel. exposure SBS5", size=1, color="#7E5109", ggtheme=theme_classic(), add.params =list(color = "#FAD7A0"))
dev.off()

#SBS40
pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/correlation_signatures_withclinical/results_postREVISION_22.10.04/plots/SBS40correlationwithAGE_pearson_22.10.04.pdf", width = 10, height = 7)
ggscatter(DT_ALL_PATIENTS, x = "age", y = "SBS40", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "age (days)", ylab = "Rel. exposure SBS40", size=1, color="#1B4F72", ggtheme=theme_classic(), add.params =list(color = "#85C1E9"))
dev.off()

