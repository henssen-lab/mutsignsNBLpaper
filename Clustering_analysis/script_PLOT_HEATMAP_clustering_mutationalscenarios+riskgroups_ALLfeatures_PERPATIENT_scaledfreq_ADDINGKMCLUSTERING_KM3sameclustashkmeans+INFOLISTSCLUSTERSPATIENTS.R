#####################################################################
#### SCRIPT HEATMAP TO SHOW CLUSTERING OF PATIENTS - SCALED FREQ ####
#####################################################################

##### SCRIPT TO PLOT HEATMAP WITH ALL THE SCALED FREQUENCIES WE USE TO DO THE CLUSTERING, PER PATIENT, PER SCENARIO/RISK GROUP

#libraries
library(ComplexHeatmap)
library(circlize)
library(stringr)

#open list with all the KM pvalues
list_distances <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUST_DISTANCES_perPATIENT_clustering_k3_hkmeansmethod_22.08.17.txt")
rownames(list_distances) <- list_distances[,1]
list_distances <- list_distances[,-1]

### PLOT PER CLUSTER
col_fun = colorRamp2(c(-1,0,1), c("gray75","white","darkred"))
risk_group <- list_distances$risk_group
annotation <- HeatmapAnnotation(risk_group = risk_group, col = list(risk_group = c("HR_MNA" = "#CF8484", "HR_non_MNA" = "#976E97", "non_HR" = "#9CB4CB")), annotation_height = unit(c(3), "mm"), annotation_legend_param = list(risk_group = list(title = "Risk Group")), gap=unit(7,"points"))

#Clustering
CL1 <- list_distances[,-c(33,34)]
MAT_dist_CL1 <- t(CL1)
#each time you run it, the clustering is different
#I first run the hkmeans clustering (with the recommended permutations) and then I run this until I get the correct clustering
HT1 <- Heatmap(MAT_dist_CL1, row_order =c(4,3,1,2,5,6,7,8,9,10,11,12,13,14,15,16,18,17,19,21,23,22,20,24,25,26,27,28,29,30,31,32),name = "scaled freq", col = col_fun, cluster_rows=FALSE, cluster_columns=TRUE, row_names_side = "left", column_names_side = "bottom", heatmap_legend_param = list(at=c(-1,0,1),legend_height = unit(4, "cm")),clustering_method_columns="ward.D2",clustering_distance_columns="euclidean",column_names_gp=gpar(fontsize = 8),top_annotation=annotation,column_km=3)
pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/plots/PLOT_HEATMAP_3mutscenarios_scaledFREQ_ALLfeatures_ALLpatients_sign+CR_dendogram_HC+KM3clust.pdf", width = 10, height = 7)
HT1
dev.off()

#Look at table distribution of risk groups in clusters
ord <- column_order(HT1) #VERY IMPORTANT - each time we run column_order is like re-doing the clustering, we obtain a different clustering - MAKE SURE IT CORRESPONDS TO THE PLOT
cnames <- data.frame(colnames(MAT_dist_CL1))
cnames$cluster <- "no"
names(cnames) <- c("sample_id", "cluster")

for (i in 1:length(ord)){
  for (x in 1:dim(as.data.frame(ord[i]))[1]){
    coeff <- as.data.frame(ord[i])[x,1]
    clust_name <- unlist(strsplit(names(as.data.frame(ord[i])),"X"))[2]
    cnames$cluster[coeff] <- paste("Cluster",clust_name,sep="")
  }
}
#Generate the table
risk_df <- data.frame(list_distances$risk_group)
risk_df$sample_id <- (rownames(list_distances))
names(risk_df) <- c("risk_group", "sample_id")
clusters_patients_RISK <- merge(cnames, risk_df, by="sample_id")
table_count_clust <- table(clusters_patients_RISK$cluster, clusters_patients_RISK$risk_group)
table_count_clust
#Save results
write.table(table_count_clust, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUSTER_TABLE_classification_patients_riskgroups_k3_hc+kmeansCOMPLEXHEATMAPmethodRESULTS.txt",sep=""),sep="\t", quote=FALSE)
write.table(clusters_patients_RISK, paste("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/CLUSTER_RESULTSPATIENTS_classification_patients_riskgroups_k3_hc+kmeansCOMPLEXHEATMAPmethodRESULTS.txt",sep=""),sep="\t", quote=FALSE,row.names=FALSE)
