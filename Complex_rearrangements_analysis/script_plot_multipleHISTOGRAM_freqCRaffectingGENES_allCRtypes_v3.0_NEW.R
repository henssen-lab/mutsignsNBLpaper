############################################################################################
#### SCRIPT TO PLOT HISTOGRAM SHOWING THE FREQUENCY OF CR AFFECTING GENES IN OUR COHORT ####
############################################################################################


#libraries
library(reshape2)
library(ggplot2)
library(wesanderson)
library(forcats)

### Arguments
args <- commandArgs(trailingOnly=TRUE)

typeofrun <- args[1] # select the run between FILTEREDATLEAST2CALLERS, PASS, or nofiltered
AAversion <- paste("+AA", args[2], sep="") # select between forced, or normal

INPUT = NULL


############# OPEN AND PREPARE INPUT FILES #############

DF_genesCR <- read.table(paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/GENES_analysis_CR/table_GENESoverlapingcomplexrearrangements_INTERSECT_114patients_infopatient_inforiskgroup_proteincoding.txt",sep=""),sep="\t")
names(DF_genesCR) <- c("gene", "CR_type", "risk_group", "patient_id")
DF_genesCR <- unique(DF_genesCR[grep(DF_genesCR$gene,pattern="HLA-",invert=TRUE),])

freq_table_CR_genes <- data.frame(table(DF_genesCR$gene, DF_genesCR$CR_type))
names(freq_table_CR_genes) <- c("gene", "CR_type", "Freq")

#gene refs
#genes1 <- read.csv("/data/gpfs-1/users/rodrigue_c/work/refs/cancer_gene_census.csv")
#genes1 <- data.frame(unique(genes1[,1]))
#names(genes1) <- "gene"
genes2 <- read.table("/data/gpfs-1/users/rodrigue_c/work/refs/KNOWNgenes_DNArepair.txt")
names(genes2) <- "gene"
genes3 <- read.table("/data/gpfs-1/users/rodrigue_c/work/refs/KNOWNgenes_NBL_ALL.txt")
names(genes3) <- "gene"

#genes_all <- unique(rbind(genes1, genes2, genes3))
genes_all <- unique(rbind(genes2, genes3))


###### HISTOGRAM TOP GENES AFFECTED PER CR TYPE - TOP 50  ######
types_CR <- unique(DF_genesCR$CR_type)

#GENERATE INPUT PLOT
for (i in 1:length(types_CR)){
  DF_genesCR_select <- DF_genesCR[DF_genesCR$CR_type==types_CR[i],]
  freq_table_only_genes_CR <- data.frame(table(DF_genesCR_select$gene))
  freq_table_only_genes_M <- merge(genes_all,freq_table_only_genes_CR,by.x="gene", by.y="Var1")
  freq_table_only_genes_O <-freq_table_only_genes_M[order(freq_table_only_genes_M$Freq, decreasing=TRUE),]
  freq_table_only_genes_TOP <- freq_table_only_genes_O
  input1 <- merge(freq_table_only_genes_TOP,DF_genesCR_select,by="gene")
  INPUT <- rbind(INPUT, input1)
}
INPUT$Freq.y <- 1
INPUT$Freq.p <- round(INPUT$Freq*100/114,digits=2) #calculate the % of patients


#PLOT num patients
#PLOTa <- ggplot(INPUT, aes(x=fct_reorder2(gene,gene,gene,.desc = TRUE), fill=CR_type)) + geom_bar() + coord_flip() + facet_grid(.~CR_type) + scale_x_discrete(name=" NBL + repair Genes") + scale_y_continuous(name="Frequency of CR (# patients)") + labs(title = "Frequency of CR types affecting genes" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "CR type")) + theme_linedraw() + theme()
#save
#ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_HISTOGRAM_frequencyofpatients_GENESaffectedbyCR_114patients_NBL+repairgenes_numpatients.pdf",sep=""), width = 10, height = 7, plot=PLOTa)


#PLOT % patients
INPUT2 <- unique(INPUT[,c(1,3,7)])
PLOTb <- ggplot(INPUT2, aes(x=fct_reorder2(gene,gene,gene,.desc = TRUE), fill=CR_type)) + geom_bar(aes(y = Freq.p),stat="identity") + coord_flip() + facet_grid(.~CR_type) + scale_x_discrete(name=" NBL + repair Genes") + scale_y_continuous(name="Frequency of CR (% patients)") + labs(title = "Frequency of CR types affecting genes" ,subtitle = "Peifer+Berlin cohorts" ) + guides(fill = guide_legend(title = "CR type")) + theme_linedraw() + theme()
ggsave(filename=paste("/fast/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_",typeofrun,AAversion,"/CR_allplots/Plot_HISTOGRAM_frequencyofpatients_GENESaffectedbyCR_114patients_NBL+repairgenes_percentpatients.pdf",sep=""), width = 10, height = 7, plot=PLOTb)


