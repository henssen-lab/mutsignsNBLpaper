########################################################################
#### SCRIPT TO PLOT VIOLIN PLOTS COMPLEX REARRANGEMENTS PER CR TYPE ####
########################################################################


#libraries
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(wesanderson)
library(ggstatsplot)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE INPUT FILES #############

#counts
MATRIX_complex_svs_M <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/MATRIX_complexrearrangements_types_COUNTS_perpatient_114patients_RISKGROUPinfo_",typeofrun,AAversion,".txt",sep=""),sep="\t")

MATRIX_complex_svs_M$all_CR <- rowSums(MATRIX_complex_svs_M[,c(2:10)])
MATRIX_complex_svs_M <- MATRIX_complex_svs_M[,c(1:10,12,11)]

###### VIOLIN PLOTS DISTRIBUTION COUNTS OF CR PER RISK GROUP  ######

### ALL CR TYPES

for (c in 2:(length(names(MATRIX_complex_svs_M))-1)){
	crtype <- names(MATRIX_complex_svs_M)[c]
	PLOT_CR <- ggstatsplot::ggbetweenstats(data = MATRIX_complex_svs_M, x = risk_group,y = !!ensym(crtype), pairwise.comparisons = TRUE, pairwise.display = "sig", p.adjust.method = "fdr", mean.plotting=FALSE, title = paste("Distribution num. of ",crtype," in the different NBL risk groups",sep=""), var.equal=FALSE, ggtheme = ggplot2::theme_classic(), type="np", violin.args = list(width = 0.5, alpha = 0, colour=NA), point.args = list(position = ggplot2::position_jitterdodge(dodge.width = 0.60),alpha = 0.4, size = 5, stroke = 0)) + geom_violin(aes(fill = risk_group, colour=risk_group), width = 0.5, alpha = 0.2) + ggplot2::scale_fill_manual(values=c("#CC3333", "#993399", "#6699CC")) + ggplot2::scale_colour_manual(values=alpha(c("#CC3333", "#993399", "#6699CC"),0.2))
	ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/VIOLINPlots_complexrearrangements_type",crtype,"_COUNTS_compare_RISKGROUP_",typeofrun,AAversion,".pdf",sep=""), width = 10, height = 7, plot=PLOT_CR)
}

