#########################################################################
###                      SIGNATURE ANALYSIS                           ###
### 1ST SCRIPT: GENERATE FILE WITH TRINUCLEOTIDE CONTEXT FOR ALL SNVS ###
#########################################################################


#libraries
library(BSgenome.Hsapiens.UCSC.hg19)
library(mutSignatures)
library(stringr)
library(readxl)
library(reshape2)
library(ggplot2)
library(wesanderson)

#Reference genome hg19
hg19 <- BSgenome.Hsapiens.UCSC.hg19


###### OPEN INPUT FILES ######

### Open file with patient_id type and coverage
patients_type_file <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/clinical_data/patients_risk_groups_simplified_Peifer+Berlin_infodataavail_v2.2.txt",sep="")
names(patients_type_file) <- c("sample_id", "cohort", "risk_group")
patients_type_file <- patients_type_file[,c(1:3)]

### Open files from intersect with genes (mutation info + gene mutated)
#SNVs
x_snvs <- read.table("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_21.04.12_ALLSNVS/RESULTS_all_patients_berlin+peifer_cohorts_SNVS_mutect_PASS_21.01.26.txt.bedformat", sep="")
x_snvs <- x_snvs[,c(1,2,11,13,14,11,11,9,11,11,4,10)]
names(x_snvs) <- c("chr", "pos_bkp1", "id", "ref", "alt", "qual", "filter", "info", "format", "xtr1", "sample_id", "mut_type")

#Exclude patients
x_snvs <- x_snvs[grep(x_snvs$sample_id,pattern="CB2044|NBL47|NBL49|NBL50|NBL53|NBL54",invert=TRUE),]


###### SIGNATURE ANALYSIS ######

### Prepare data for signature analysis
#Filter non SNV
x_snvs_filt <- filterSNV(dataSet = x_snvs,seq_colNames = c("ref","alt"))
#Attach nucleotide context: 1 nucleotide before, 1 after
source("/fast/users/rodrigue_c/work/scripts/function_attachContextelias_mutsignatures.R")
x_snvs_filt_CONTEXT <- attachContextelias(mutData = x_snvs_filt, chr_colName = "chr", start_colName = "pos_bkp1", end_colName = "pos_bkp1", nucl_contextN = 3, BSGenomeDb = hg19, context_colName="context")
#Remove mismatches
x_snvs_filt_CONTEXT_2 <- removeMismatchMut(mutData = x_snvs_filt_CONTEXT, refMut_colName = "ref", context_colName = "context", refMut_format = "N")

### Assess mutation type
#Compute mutType (trinucleotides for the signature ex.A[C]T)
x_snvs_muttype <- attachMutType(mutData = x_snvs_filt_CONTEXT_2, ref_colName = "ref", var_colName = "alt", context_colName = "context")
#Mutation Counts
x_snvs_muttype_counts <- countMutTypes(mutTable = x_snvs_muttype, mutType_colName = "mutType", sample_colName = "sample_id")
save(x_snvs_muttype_counts, file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/ALLSNVS_mutation_type_counts_TRInucleotidecontext_22.06.13.RData")

