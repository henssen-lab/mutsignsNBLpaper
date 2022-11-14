###################################################################################
#### SCRIPT TO PLOT PIECHARTS COMPLEX REARRANGEMENTS DISTRIBUTION WHOLE COHORT ####
###################################################################################

#libraries
library(scales)
library(ggplot2)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE MUTATION FILES #############

#complex rearrangements counts
complex_svs <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_COUNTS_perpatient_114patients_RISKGROUPinfo_FILTEREDATLEAST2CALLERS+AAnormal.txt",sep=""), sep="\t")

all_muts <- unique(complex_svs)
all_muts <- all_muts[,-c(1,dim(all_muts)[2])]


############# CREATE DF #############

#COUNTS OF CR IN COHORT - we count the total number of complex rearrangements in the cohort
all_muts_COUNTS <- data.frame(colSums(all_muts))
all_muts_COUNTS$type <- rownames(all_muts_COUNTS)
names(all_muts_COUNTS) <- c("freq", "type")

#COUNTS OF PATIENTS WITH CR IN COHORT - we count the total number of patients having complex rearrangements in the cohort
all_muts_p <- all_muts
all_muts_p[all_muts_p>0] <- 1
all_muts_PAT <- data.frame(colSums(all_muts_p))
all_muts_PAT$type <- rownames(all_muts_PAT)
names(all_muts_PAT) <- c("freq", "type")

############# PLOT #############

pdf(file=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/PIECHART_distribution_CR_across_WHOLECOHORT_countsCR+countspatients_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 15)

#COUNTS OF CR IN COHORT - we count the total number of complex rearrangements in the cohort
pie_labels <- paste0(all_muts_COUNTS$type, " : ", round(100 * all_muts_COUNTS$freq/sum(all_muts_COUNTS$freq), 2), "%")
plot1 <- pie(x=all_muts_COUNTS$freq,labels=pie_labels,main="Distribution of complex rearrangements across the cohort\n(114 patients) - count of CR")

#COUNTS OF PATIENTS WITH CR IN COHORT - we count the total number of patients having complex rearrangements in the cohort
pie_labels <- paste0(all_muts_PAT$type, " : ", round(100 * all_muts_PAT$freq/sum(all_muts_PAT$freq), 2), "%")
plot1 <- pie(x=all_muts_PAT$freq,labels=pie_labels,main="Distribution of complex rearrangements across the cohort\n(114 patients) - count of patients having CR")

dev.off()

