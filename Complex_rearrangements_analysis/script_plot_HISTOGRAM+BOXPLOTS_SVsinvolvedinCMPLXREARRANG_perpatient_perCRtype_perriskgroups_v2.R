####################################################################
###### PLOT PROPORTION SVS INVOLVED IN COMPLEX REARRANGEMENTS ######
####################################################################

#Libraries
library(ggplot2)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE CLINICAL FILES #############

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v3.0.txt",sep="")

#Format
patients_type_file <- patients_type_file[,c(1,3)]
names(patients_type_file) <- c("sample_id", "risk_group")

#Exclude patients
patients_type_file <- patients_type_file[grep(patients_type_file$sample_id,pattern="CB2044|NBL47|NBL53|NBL54|NBL49|NBL50",invert=TRUE),]


############# OPEN AND PREPARE INPUT FILES #############

### Open files
#input plot 2
SV_counts1 <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/INTERSECT_SVS_complexrearrangements/COUNTS_SVS_implicated_complexrearrangements_PERCOMPLXREARRANGTYPE_WHOLECOHORT.txt",sep=""),sep="\t")

#input plot 1
SV_counts2 <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/INTERSECT_SVS_complexrearrangements/COUNTS_SVS_implicated_complexrearrangements_PERPATIENT_PER_RISK_GROUP.txt",sep=""),sep="\t")
#get counts for all patients (not only the ones having CR)
SV_counts2_M <- merge(patients_type_file,SV_counts2,by.x="sample_id",by.y="patient_id",all.x=TRUE)
SV_counts2_M <- SV_counts2_M[,c(1,2,4)]
SV_counts2_M[is.na(SV_counts2_M$proportion_SVs_patient),3] <- 0


############# PLOTS #############

###PLOT BOXPLOTS RISK GROUPS

#all SVs
PLOT_SVS <- ggstatsplot::ggbetweenstats(data = SV_counts2_M, x = risk_group,y = proportion_SVs_patient, pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = "Proportion SVs in complex rearrangements per risk groups", var.equal=FALSE,ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = risk_group, colour=risk_group), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2)) + geom_hline(yintercept = 10)
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_BOXPLOT_proportionSVsincomplexrearrang_114patients_COUNTPERPATIENT_risk_groups.pdf",sep=""), width = 10, height = 10, plot=PLOT_SVS)


###PLOT HISTOGRAM %OF SVS IN EACH COMP REARRANG TYPE

PLOTa <- ggplot(SV_counts1, aes(x=type_CR, y=percent_SVs_CR))
PLOTb <- PLOTa + geom_bar(position="stack", stat="identity") + scale_y_continuous(name="% of total SVs") + scale_x_discrete(name="Type of complex rearrangement") + labs(title = "Percentage of SVs for each type of Complex Rearrangement - Delly2+Svaba+Novobreak ATLEAST2CALLERS" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "Patient type")) + theme_classic() + theme()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_percentageSVs_ineachcomplexrearrangTYPE_WHOLECOHORT.pdf",sep=""), width = 7, height = 7, plot=PLOTb)
