###########################################################################################################################
#			CORRELATION cohorts
###########################################################################################################################
data <- read.table("figure4.csv", sep="\t", header=TRUE, dec=',')
colnames(data)[1] <- "signature"

pvaule <- data$p_value
round(p.adjust(pvaule, "BH"), 3)
data$pvalueadj <- round(p.adjust(pvaule, "BH"), 3)

pdf(file="mimi2.pdf", width = 8, height = 8)
ggplot(data, aes(x=PeDiOn, y=Peifer.Berlin)) + geom_point(aes(col = risk_group, size=pvalueadj)) + geom_smooth(method=lm , color="red", fill="#69b3a2", se=TRUE) + scale_size_continuous(range = c(2, 7), limits = c(.001, .9), breaks = c(.01, .05, .1), labels = c("0.01", "0.05", "0.1")) + geom_text(label=data$signature, nudge_x = 1.55, nudge_y = 2.2, check_overlap = FALSE) + scale_x_continuous(limits = c(0, 70)) + scale_y_continuous(limits = c(0, 70)) + theme_Publication()
dev.off()


