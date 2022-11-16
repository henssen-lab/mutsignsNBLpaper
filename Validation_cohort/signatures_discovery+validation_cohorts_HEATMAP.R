###########################################################################################################################
# 			HEATMAP
###########################################################################################################################
load(mis_Objetos.RData)

heat_meli <- read.table("figure3.csv",sep=";")
names(heat_meli) <- c("Y", "X", "Z")
heat_meli$Y <- factor(heat_meli$Y,levels = c("SBS3","SBS5","SBS18","SBS40","ID1","ID2","ID4","ID5","ID6","ID8","ID9","R2","R6a","R6b","CX1","CX2","CX3","CX5","CX7","CX11","CX14","CX15"))
all$sig <- factor(all$sig, levels = c("SBS3", "SBS5", "SBS18", "SBS40"))

pdf(file="correlation_signatures_cohorts_peifer_berlin_pedion_2.pdf", width = 4, height = 14)
ggplot(heat_meli, aes(X, Y, fill= Z)) + geom_tile(color = "white", lwd = 1.2,linetype = 1) +scale_colour_Publication()+ theme_Publication(base_size=8) + scale_fill_distiller(palette = "RdPu") + theme(legend.position = "none")
dev.off()