##########################################################
#### SCRIPT HEATMAP CLUSTERING - MUTATIONAL SCENARIOS ####
##########################################################

##### SCRIPT TO PLOT HEATMAP FOR ALL THE CLUSTERING FEATURES WITH THE CLUSTERING MEAN DISTANCE

#libraries
library(ComplexHeatmap)
library(circlize)
library(stringr)

#open list with all the KM pvalues
list_distances <- read.delim("/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/SUMMARY_STATISTICS_clustering_k3_hkmeansmethod_22.08.17.txt")
list_distances <- list_distances[,-1]
MAT_dist <- t(list_distances)
colnames(MAT_dist) <- c("Scenario_1","Scenario_3","Scenario_2")

### PLOT HEATMAP FOR ADJ PVALUES
col_fun = colorRamp2(c(-1,0,1), c("gray75","white","darkred"))
pdf(file="/data/gpfs-1/users/rodrigue_c/work/landscape_NBL_project/wgs_analysis/clustering_analysis/newresults_CMPLXREARRNG_ATLEAST2CALLERS+normal/plots/PLOT_HEATMAP_3mutscenarios_meandistance_ALLfeatures_sign+CR.pdf", width = 4, height = 10)
Heatmap(MAT_dist, row_order =c(4,3,1,2,5,6,7,8,9,10,11,12,13,14,15,16,18,17,19,21,23,22,20,24,25,26,27,28,29,30,31,32), column_order=c(1,3,2),name = "mean\nclust distance", col = col_fun, cluster_rows=FALSE, cluster_columns=FALSE, row_names_side = "left", column_names_side = "top", heatmap_legend_param = list(at=c(-1,0,1),legend_height = unit(4, "cm")))
dev.off()
