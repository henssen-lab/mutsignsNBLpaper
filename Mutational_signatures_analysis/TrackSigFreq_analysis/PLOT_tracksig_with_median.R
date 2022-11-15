###########################################################################################################################
#			TRACKSIG WITH MEDIAN
###########################################################################################################################
load(mis_Objetos.RData)

pdf(file="tracksig_W_medianhand_riskgroup.pdf", width = 11, height = 8)

Line <- ggplot(all2, aes(x=mixture, y=exposure, colour=sig, shape = sample, group=interaction(sig, sample))) + geom_line(size=0.5, alpha=0.30) + facet_wrap(~group, nrow =3) + labs(title="PeDiOn SBS") 
#grid.arrange(Line,(Line +scale_colour_Publication()+ theme_Publication()),nrow=1)
Line1 <- Line + geom_line(aes(size = sample)) + scale_size_manual(values = c( median = 0.9)) +scale_colour_Publication()+ theme_Publication() + xlim(1.0, 0) + labs(x = "Cancer cell fraction" , y = "Signature exposure(%)")  + scale_color_manual(values=c('#73C8AE','#FBE596','#333467','#F4794E'))# SI
#Line1 <- Line +scale_colour_Publication()+ theme_Publication() + xlim(1.25, 0) + labs(x = "Cancer cell fraction" , y = "Signature exposure(%)")  + scale_color_manual(values=c('#73C8AE','#FBE596','#333467','#F4794E'))
#Line2 <- Line1 + stat_summary(aes(group=sig), geom="line", fun = "mean", size=0.8, alpha=0.50)
Line1
dev.off()
c

data <- read.table("median_tracksig.csv", sep="\t", header=TRUE, dec=",")
colnames(data)[1] <- "exposure"
all2 <- rbind(all,data)


Line <- ggplot(all, aes(x=mixture, y=exposure, colour=sig, shape = sample, group=interaction(sig, sample))) + geom_line(size=0.5, alpha=0.40) + facet_wrap(~group, nrow =3) + labs(title="PeDiOn SBS") 
#grid.arrange(Line,(Line +scale_colour_Publication()+ theme_Publication()),nrow=1)
Line1 <- Line +scale_colour_Publication()+ theme_Publication() + xlim(1.25, 0) + labs(x = "Cancer cell fraction" , y = "Signature exposure(%)")  + scale_color_manual(values=c('#73C8AE','#FBE596','#333467','#F4794E'))
Line2 <- Line1 + stat_summary(aes(group=sig), geom="line", fun = "mean", size=0.8, alpha=0.50)
