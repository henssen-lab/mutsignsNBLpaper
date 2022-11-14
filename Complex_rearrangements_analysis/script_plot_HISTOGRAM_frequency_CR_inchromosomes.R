########################################################################################
#### SCRIPT TO PLOT HISTOGRAM COMPLEX REARRANGEMENTS ACROSS CHROMOSOMES PER CR TYPE ####
########################################################################################


#libraries
library(stringr)
library(reshape2)
library(ggplot2)
library(wesanderson)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE INPUT FILES #############

DF_coords_rearrang <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/table_ALLregions_ALLcomplexrearrangements_114patients_infopatient_inforiskgroup.txt",sep=""),sep="\t")
DF_coords_rearrang <- unique(DF_coords_rearrang[,c(1,4,6)])

#frequency table counting per patient
TABLE_coords_rearrang <- data.frame(table(DF_coords_rearrang$chr, DF_coords_rearrang$CR_type))
names(TABLE_coords_rearrang) <- c("chr", "CR_type", "Freq")


###### HISTOGRAM PLOT DISTRIBUTION OF COMPLEX REARRANGEMENTS ACROSS THE CHROMOSOMES  ######

TABLE_coords_rearrang$chr <- factor(TABLE_coords_rearrang$chr, levels= str_sort(unique(TABLE_coords_rearrang$chr),numeric=TRUE))

#plot
PLOTa <- ggplot(TABLE_coords_rearrang, aes(x=chr, y=Freq, fill=CR_type))
PLOTb <- PLOTa + geom_bar(position="stack", stat="identity") + scale_y_continuous(name="Freq of CR per chrom in patients") + scale_x_discrete(name="Chromosome") + labs(title = "Frequency of CR types across all the chromosomes (patient based)" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "CR type")) + theme_classic() + theme()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_HISTOGRAM_CRfrequencyinCHROMOSOMES.pdf",sep=""), width = 7, height = 9, plot=PLOTb)
