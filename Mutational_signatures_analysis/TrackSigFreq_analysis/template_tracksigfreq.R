###### R SESSION ######
R
#TEMPLATE

#libraries
library(TrackSig)
library(ggplot2)
library(optparse)

option_list = list(
  make_option(c("-v", "--vcffile"), type="character", default=NULL, help="VCF file", metavar="character"),
  make_option(c("-c", "--cnafile"), type="character", default=NULL, help="CNA file", metavar="character"),
  make_option(c("-p", "--purity"), type="numeric", default=NULL, help="purity tumor", metavar="character"),
  make_option(c("--file1"), type="character", default=NULL, help="name plot 1", metavar="character"),
  make_option(c("--file2"), type="character", default=NULL, help="name plot 2", metavar="character"),
  make_option(c("--title"), type="character", default=NULL, help="tittle plot", metavar="character"),
  make_option(c("--sampleid"), type="character", default=NULL, help="Sample ID", metavar="character")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

VCFFILE = opt$vcffile
CNAFILE = opt$cnafile
PURITY = opt$purity
PLOT1 = opt$file1
PLOT2 = opt$file2
TITLE = opt$title
SAMPLEID = opt$sampleid

###### OPEN INPUT FILES ######
### Open files vcf, cna for each patient
vcfFile = file.path(VCFFILE)
cnaFile = file.path(CNAFILE)

##put the purity number (default 1)
purity = PURITY
# SBS3, SBS5, SBS18, SBS40
elias <- TrackSig:::alex_merged[,c(2,4,13,35), drop=FALSE]

## First, restrict the list of signatures to fit exposure for (threshold) signatures with exposure under this across all timepoints will not be fit
detectedSigs <- detectActiveSignatures(vcfFile = vcfFile, cnaFile = cnaFile,
                                       purity = purity, threshold = 0, referenceSignatures = elias) #0 cause we want to see all the signatures

## Compute the trajectory for all timepoints.
set.seed(1224)

traj <- TrackSig(sampleID = "$i", activeInSample = detectedSigs,
                 vcfFile = vcfFile, cnaFile = cnaFile, purity = purity, referenceSignatures = elias)

# function TrackSig has three available methods for segmentation, controlled by the parameter scoreMethod
	#Signature (described in the TrackSig paper)
	#SigFreq (described in the TrackSigFreq paper)
	#Frequency (not explicitly described, but corresponds to the frequency likelihood in the TrackSigFreq paper).


###### PLOT TRAJECTORY ######
## trajectory with linear x-axis
pdf(file=PLOT1, width = 10, height = 4)
plotTrajectory(traj, linearX = T) + labs(title = TITLE)

dev.off()

## trajectory with non-linear x-axis
pdf(file=PLOT2, width = 10, height = 4)
nonLinPlot <- plotTrajectory(traj, linearX = F, anmac = T) + labs(title = TITLE)
addPhiHist(traj, nonLinPlot)
dev.off()


















