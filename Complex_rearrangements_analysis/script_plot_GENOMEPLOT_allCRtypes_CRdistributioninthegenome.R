########################################################################
#### SCRIPT TO PLOT VIOLIN PLOTS COMPLEX REARRANGEMENTS PER CR TYPE ####
########################################################################


#libraries
library(stringr)
library(reshape2)
library(ggplot2)
library(viridis)
library(karyoploteR)
library(GenomicRanges)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal


############# OPEN AND PREPARE INPUT FILES #############

DF_coords_rearrang <- read.delim(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/table_ALLregions_ALLcomplexrearrangements_114patients_infopatient_inforiskgroup.txt",sep=""),sep="\t")

types_CR_list <- sort(unique(DF_coords_rearrang$CR_type))


###### GENOME PLOT DISTRIBUTION OF COMPLEX REARRANGEMENTS ACROSS THE GENOME  ######


pdf(file=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/GENOME_PLOT_densityofregionsaffectedbyCR_acrossgenome_allCR_perCRtype_114patients.pdf",sep=""), width = 15, height = 6)

#parameters
pp <- getDefaultPlotParams(plot.type = 4)
pp$leftmargin <- 0.15
pp$topmargin <- 15
pp$bottommargin <- 15
pp$ideogramheight <- 5
pp$data1inmargin <- 10
#plot chromosomes and bands
kp <- plotKaryotype(genome="hg19", plot.type=4, ideogram.plotter = NULL, labels.plotter = NULL, plot.params = pp, main="Complex Rearrangements Density (1Mb window)")
kpAddCytobandsAsLine(kp,color.schema="biovizbase",lwd=7)
kpAddChromosomeNames(kp,chr.names=c(1:22,"X","Y"))
for (i in 1:length(types_CR_list)){
	at <- autotrack(current.track = i, total.tracks = length(types_CR_list))
	kpDataBackground(kp, r0=at$r0, r1=at$r1, color = NULL)
	GR_coords_rearrang <- makeGRangesFromDataFrame(DF_coords_rearrang[DF_coords_rearrang$CR_type==types_CR_list[i],])
	kp1 <- kpPlotDensity(kp, GR_coords_rearrang, window.size = 1e6, r0=at$r0, r1=at$r1, col=rainbow(length(types_CR_list))[i], border=transparent(rainbow(length(types_CR_list))[i],amount=0.90))
	computed.ymax <- ceiling(kp1$latest.plot$computed.values$max.density)
	kpAxis(kp1, ymin=0, ymax=computed.ymax, numticks = 2, r0=at$r0, r1=at$r1)	
	kpAddLabels(kp1, labels = types_CR_list[i], r0=at$r0, r1=at$r1)
	#kpAxis(kp)	
}
kp1
dev.off()
