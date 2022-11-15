
#libraries
library(factoextra)
library(cluster)
library(clValid)
library(fpc)
library(NbClust)
library(dplyr)
library(ClustImpute)

#arguments
args <- commandArgs(TRUE)

k_clust <- as.numeric(args[1]) #3 number of clusters for the analysis


#### OPEN INPUT FILES ####

#SNV signatures
SNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/mut_signatures_analysis/analysis_22.06.13_ALLSNVS_COSMIC+SIGNAL/exposure_SNV_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt", sep="\t", stringsAsFactors=FALSE)
risk_DF <- SNVsign[,c(1,dim(SNVsign)[2])]
SNVsign <- SNVsign[,-dim(SNVsign)[2]]

#INDELS signatures
INDsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/indels_signatures_analysis/analysis_22.06.13_ALLINDELS_perpatient_perriskgroup/exposure_INDELS_COSMIC_signature_absexp_per_patient_FREQ_22.06.13.txt", sep="\t", stringsAsFactors=FALSE)
INDsign <- INDsign[,-dim(INDsign)[2]]

#SV signatures
SVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis/sv_Palimpsest_signatures_analysis/SV_signatures_ALLpatients_2022.06.15/exposure_SV_denovo+SIGNAL_signature_absexp_per_patient_FREQ_22.06.13.txt", sep="\t", stringsAsFactors=FALSE)
SVsign <- SVsign[,-dim(SVsign)[2]]

#CNV signatures
CNVsign <- read.delim("/fast/work/users/rodrigue_c/landscape_NBL_project/wgs_analysis//cnv_signatures_analysis/cnv_signature_analysis_NEWfilters_22.06.15/exposure_CNA_denovo_signature_freqexp_per_patient_FREQ_22.06.13_REFITFREQ.txt", sep="\t", stringsAsFactors=FALSE)
CNVsign <- CNVsign[,-dim(CNVsign)[2]]

#Complex rearrangements
comprearrang <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/jabba_results/jabba_results_FILTEREDATLEAST2CALLERS+AAnormal/MATRIX_complexrearrangements_types_FREQ_perpatient_114patients_RISKGROUPinfo_FILTEREDATLEAST2CALLERS+AAnormal.txt", sep="\t", stringsAsFactors=FALSE)
comprearrang <- comprearrang[,-dim(comprearrang)[2]]


#### FORMAT INPUT FILES ####

#Merge everything by sample_id
m1 <- merge(SNVsign, INDsign, by="sample_id", all=TRUE)
m2 <- merge(m1, SVsign, by="sample_id", all=TRUE)
m3 <- merge(m2, CNVsign, by="sample_id", all=TRUE)
merge_final <- merge(m3, comprearrang, by="sample_id", all=TRUE)

#Format (only numeric values)
rownames(merge_final) <- merge_final$sample_id
df_cluster_input <- merge_final[,-1]
#df_cluster_input[is.na(df_cluster_input)] <-0 #transform NAs from CNVsign into 0s for the patients that don't have results
df_cluster_input_N <- ClustImpute(df_cluster_input,nr_cluster=k_clust, nr_iter=25, c_steps=25, n_end=10) #In case we want to impute the missing values
df_cluster_input <- df_cluster_input_N$complete_data

#### PROCESS INPUT DATA ####

#Filter out non-informative signatures+complex rearrangements
#The ones not exhibiting significant differences between risk groups (ID8, SV2, SV4, CX2, CX15, BFB, chromoplexy, chromothripsis, pyrgo, rigma)
#df_cluster_input_FILT <- df_cluster_input[,-which(names(df_cluster_input) %in% c("ID8","SV_denovo_2","SV_denovo_4","CX2","CX15","chromothripsis","Clustered.rearrangement","Complex.non.cyclic","pyrgo","rigma"))]
df_cluster_input_FILT <- df_cluster_input
#data standarization
df_cluster_input_FILT <- data.frame(scale(df_cluster_input_FILT))


#### HOW MANY CLUSTERS? ####

# Open pdf file
pdf(file= paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/plots/Assess_BESTnumberofclusters_runk",k_clust,"_22.08.17.pdf",sep=""))
par(mfrow= c(2,2))

#We use different methods
#gap stats
fviz_nbclust(df_cluster_input_FILT, hkmeans, method = c("gap_stat"),nboot=500)+ theme_classic()
#wss
fviz_nbclust(df_cluster_input_FILT, hkmeans, method = c("wss"),nboot=500) + geom_vline(xintercept = 4, linetype = 2, colour="dodgerblue4")+ theme_classic()
fviz_nbclust(df_cluster_input_FILT, hkmeans, method = c("silhouette"),nboot=500) + theme_classic()

#All methods in index including silhouette
clust_DF <- data.frame(NULL)

index <- c("kl","ch", "hartigan", "cindex", "db", "silhouette","duda", "pseudot2", "ratkowsky", "ball","ptbiserial", "frey", "mcclain", "dunn", "sdindex", "sdbw")

for (x in 1:length(index)){
  num_clust <- NbClust(data = df_cluster_input_FILT, diss = NULL, distance = "euclidean", min.nc = 2, max.nc = 10, method = "ward.D2",index =index[x])$Best.nc[1]
  res_clust <- c(index[x], num_clust)
  clust_DF <- rbind(clust_DF,res_clust)
}

names(clust_DF) <- c("method_eval", "num_clust")
#plot
ggplot(clust_DF,aes(x=num_clust)) + geom_histogram(stat="count",fill="dodgerblue4")+ theme_classic()

dev.off()


#### CLUSTERING ANALYSIS ####

#The best number of clusters is k=3 (corresponds to risk group classification)
result_clustering.hk <-hkmeans(df_cluster_input_FILT, k_clust, hc.metric="euclidean", hc.method="ward.D2",iter.max =25)

#Generate table with counts to get coeffs
clusters_patients <- data.frame(result_clustering.hk$cluster)
names(clusters_patients) <- "clusters"
clusters_patients$sample_id <- rownames(clusters_patients)
clusters_patients_RISK <- merge(clusters_patients, risk_DF, by="sample_id")
table_count_clust <- table(clusters_patients_RISK$clusters, clusters_patients_RISK$risk_group)
#Coeffs of correctness
table_exito_coeff <- t(data.frame(table_count_clust[1,]))
rownames(table_exito_coeff) <- NULL
for (t in 1:dim(table_count_clust)[2]){
  exito_coeff <- max(table_count_clust[,t])*100/sum(table_count_clust[,t])
  table_exito_coeff[1,t] <- exito_coeff
}

#Table with classification cluster vs real risk groups
table_count_clust_PRE <- table_count_clust
table_count_clust_PRE
#Save results
write.table(clusters_patients_RISK, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUSTER_DATAFRAME_classification_patients_riskgroups_k",k_clust,"_hkmeansmethod_22.08.17.txt",sep=""),sep="\t", row.names=FALSE, quote=FALSE)
write.table(table_count_clust_PRE, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUSTER_TABLE_classification_patients_riskgroups_k",k_clust,"_hkmeansmethod_22.08.17.txt",sep=""),sep="\t", quote=FALSE)


# Compute cluster stats to compare with risk group classification
#The corrected Rand index provides a measure for assessing the similarity between two partitions
#Its range is -1 (no agreement) to 1 (perfect agreement)
riskgroups <- as.numeric(as.factor(clusters_patients_RISK$risk_group))
clust_stats <- cluster.stats(d = dist(df_cluster_input_FILT), riskgroups, result_clustering.hk$cluster)
#corrected rand index
exito_coeff_PRE <- clust_stats$corrected.rand

# descriptive statistics at the cluster level
summary_stats <- df_cluster_input_FILT %>% mutate(Cluster = result_clustering.hk$cluster) %>% group_by(Cluster) %>% summarise_all("mean")
summary_statsDF <- data.frame(summary_stats)
#generate table for forest plot
RESULTS_df=NULL
all_resultsDF <- data.frame(df_cluster_input_FILT %>% mutate(Cluster = result_clustering.hk$cluster) %>% group_by(Cluster))
for (n in 1:max(all_resultsDF$Cluster)){
  DF <- t(all_resultsDF[all_resultsDF$Cluster==n,-dim(all_resultsDF)[2]])
  for (l in 1:dim(DF)[1]){
    DF2 <- data.frame(DF[l,1])
    DF2$feature <- rownames(DF)[l]
    DF2$mean <- mean(DF[l,])
    DF2$lower <- min(DF[l,])
    DF2$upper <- max(DF[l,])
    DF2$group <- n
    RESULTS_df <- rbind(RESULTS_df, DF2)
    DF2 = NULL
  }
}

RESULTS_df <- RESULTS_df[,-1]

#get clustering distances per patient
clust_dist_pat <- data.frame(result_clustering.hk$data)
clust_dist_pat$sample_id <- rownames(clust_dist_pat)
clust_dist_pat_info <- merge(clust_dist_pat,clusters_patients_RISK,by="sample_id")

#Save results
write.table(summary_statsDF, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/SUMMARY_STATISTICS_clustering_k",k_clust,"_hkmeansmethod_22.08.17.txt",sep=""),sep="\t", row.names=FALSE, quote=FALSE)
write.table(RESULTS_df, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/ALLRESULTS_STATISTICS_meanupperlower_clustering_k",k_clust,"_hkmeansmethod_22.08.17.txt",sep=""),sep="\t", row.names=FALSE, quote=FALSE)
write.table(clust_dist_pat_info, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUST_DISTANCES_perPATIENT_clustering_k",k_clust,"_hkmeansmethod_22.08.17.txt",sep=""),sep="\t", row.names=TRUE, quote=FALSE)


#### PLOT CLUSTERING ####

#Silhouette of clusters
sil_plot_inp <- silhouette(result_clustering.hk$cluster, dist(df_cluster_input_FILT, method="euclidean"))
sil.data <- data.frame(cluster = factor(sil_plot_inp[, 1]), sil_width = sil_plot_inp[, 3])
sil.data$sample_id <- clusters_patients_RISK$sample_id
sil.data$risk_group <- clusters_patients_RISK$risk_group
misplaced_samples <- sil.data[sil.data$sil_width<0,]

#Plot
PLOT_sil <- ggplot(sil.data, aes(x = row.names(sil.data), y = sil_width, fill = cluster, col = cluster)) + geom_bar(stat = "identity", width = 0.5) + coord_flip() + labs(x = "") + scale_x_discrete(limits = row.names(sil.data[order(sil.data$cluster,sil.data$sil_width), ]),labels = sil.data[order(sil.data$cluster,sil.data$sil_width),3]) + theme_classic()+ theme(text = element_text(size=7))
ggsave(paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/plots/PLOT_silhouette_clustering_hkmeans_FILTEREDvariables_k",k_clust,"_22.08.17.pdf",sep=""), plot=PLOT_sil, width = 5, height = 5)


#CLUSTERS PLOT
input_plot<-df_cluster_input_FILT
input_plot$risk_group <- clusters_patients_RISK$risk_group
#dendogram - no me da los clusters bien
#fviz_dend(result_clustering.hk, k = k_clust, k_colors = c("#2E9FDF", "#00AFBB", "#E7B800", "#FC4E07"),color_labels_by_k = TRUE, cex = 0.6, palette = "jco", rect = TRUE, rect_border = "jco", rect_fill = TRUE)
#cluster
pdf(file= paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/plots/PLOT_clusters_hkmeansclustering_FILTEREDvariables_k",k_clust,"_22.08.17.pdf",sep=""))
if (k_clust==3){
  fviz_cluster(result_clustering.hk, data = input_plot, stand =TRUE, palette = "jco", repel = TRUE, ggtheme = theme_classic(), ellipse.type = "convex", shape=19, geom = "point",ellipse.alpha=0.1, show.clust.cent=FALSE) + geom_point(aes(color = input_plot$risk_group)) + scale_color_manual(values=c("#0073C2FF","#EFC000FF","#868686FF","#0073C2FF","#868686FF","#EFC000FF"))
} else {
  fviz_cluster(result_clustering.hk, data = input_plot, stand =TRUE, palette = "jco", repel = TRUE, ggtheme = theme_classic(), ellipse.type = "convex", shape=19, geom = "point",ellipse.alpha=0.1, show.clust.cent=FALSE) + geom_point(aes(color = input_plot$risk_group)) + scale_color_manual(values=c("#0073C2FF","#EFC000FF","#CD534CFF","#868686FF","#0073C2FF","#868686FF","#EFC000FF"))
}
dev.off()
