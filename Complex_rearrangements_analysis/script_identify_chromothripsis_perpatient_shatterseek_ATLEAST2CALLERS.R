### USING SHATTERSEEK TO DETECT CHROMOTHRIPSIS PER PATIENT
### FROM SV AND CNV DATA

#Rscript script_identify_chromothripsis_perpatient_shatterseek.R CB2001

#Libraries
library(ShatterSeek)
library(cowplot)
library(gridExtra)

#Arguments
args <- commandArgs(TRUE)
patient_id <- args[1]

#Initialize data.frame
FINAL_chromo = NULL

#Open cnv
cnv <-  read.delim("/fast/users/rodrigue_c/work/neuroblastoma_variant_calling_results_elias/ascat_results_hg19_martin/ALL.cnseg.tsv",sep="")
cnv_patient <- cnv[cnv$sample==patient_id,]
cnv_patient$total_cn <- cnv_patient$CNmajor + cnv_patient$CNminor
cnv_patient <- cnv_patient[,c(2,3,4,8)]
names(cnv_patient) <- c("chr", "start", "end", "total_cn")

#Collapse cn segments, script from github following their recommendation
dd <- cnv_patient
dd$total_cn[dd$total_cn == 0] <- 150000
dd$total_cn[is.na(dd$total_cn)] <- 0
library(GenomicRanges)
dd <- as(dd,"GRanges")
cov <- coverage(dd,weight = dd$total_cn)
dd1 <- as(cov,"GRanges")
dd1 <- as.data.frame(dd1)
dd1 <- dd1[dd1$score !=0,]
dd1 = dd1[,c(1,2,3,6)]
names(dd1) <- names(cnv_patient)[1:4]
dd1$total_cn[dd1$total_cn == 150000] <- 0
cnv_patient= dd1; rm(dd)

#Open sv
sv <- read.table(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_inputs/",patient_id,"/",patient_id,".ALLSVS.jabbainput.TIER1+2.bedpe",sep=""),sep="\t")
sv <- sv[sv$V12!="unknown",]
sv <- sv[sv$V1!="Y",]
sv <- sv[sv$V4!="Y",]

#change sv_type
for (s in 1:dim(sv)[1]){
    if (sv$V12[s]=="deletion"){
        sv$type[s] <- "DEL"
    } else if (sv$V12[s]=="tandem-duplication"){
        sv$type[s] <- "DUP"
    } else if (sv$V12[s]=="translocation"){
        sv$type[s] <- "TRA"
    } else if (sv$V12[s]=="inversion"){
        sv$type[s] <- "t2tINV"
    }
}

#Transform data
SV_data <- SVs(chrom1=as.character(sv$V1),
    pos1=as.numeric(sv$V2),
    chrom2=as.character(sv$V4),
    pos2=as.numeric(sv$V5),
    SVtype=as.character(sv$type),
    strand1=as.character(sv$V9),
    strand2=as.character(sv$V10))


CN_data <- CNVsegs(chrom=as.character(cnv_patient$chr),
    start=cnv_patient$start,
    end=cnv_patient$end,
    total_cn=cnv_patient$total_cn)


#Look for chromothripsis
chromothripsis <- shatterseek(SV.sample=SV_data, seg.sample=CN_data)

#Filter chromothripsis following high confident criteria
#High confidence: at least 6 interleaved intrachromosomal SVs, 7 contiguous segments oscillating between 2 CN states, the fragment joins test, and either the chromosomal enrichment or the exponential distribution of breakpoints test
chromo_filt <- chromothripsis@chromSummary
for (c in 1:dim(chromo_filt)[1]){
    if (!is.na(chromo_filt$number_DEL[c]) & !is.na(chromo_filt$number_DUP[c]) & !is.na(chromo_filt$number_h2hINV[c]) & !is.na(chromo_filt$number_t2tINV[c]) & !is.na(chromo_filt$pval_fragment_joins[c]) & (!is.na(chromo_filt$chr_breakpoint_enrichment[c]) | !is.na(chromo_filt$pval_exp_cluster[c])) & !is.na(chromo_filt$max_number_oscillating_CN_segments_2_states_chr[c])){
        if ((chromo_filt$number_DEL[c]+chromo_filt$number_DUP[c]+chromo_filt$number_h2hINV[c]+chromo_filt$number_t2tINV[c])>=6 & chromo_filt$pval_fragment_joins[c]<0.05 & (chromo_filt$chr_breakpoint_enrichment[c]<0.05 | chromo_filt$pval_exp_cluster[c]<0.05) & chromo_filt$max_number_oscillating_CN_segments_2_states_chr[c]>=7){
            line_PASS <- chromo_filt[c,]
            FINAL_chromo <- rbind(FINAL_chromo, line_PASS)
        }
    }
}

#save results if there are chromothripsis in our data
if (length(FINAL_chromo)>0){
    write.table(FINAL_chromo, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/shatterseek_chromothripsis/results_ATLEAST2CALLERS/",patient_id,"_shatterseek_chromothripsis.txt",sep=""),sep="\t",row.names=FALSE, quote=FALSE)
}


