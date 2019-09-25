#Require needed libraries
library(plyr)
library(tidyr)
library(dplyr)
library(ggplot2)
library(scales)
library(Hmisc)
theme_set(theme_minimal())
#Load and clean data

paths <- dir("data", pattern = "\\.csv$", full.names = TRUE)
mydata <- ldply(paths, read.csv)
tidyData <- separate(mydata, Condition, into = c("Condition", "Concentration"), sep = "_")

#summarize results by % of cells with autophagy (>=6 dots)
tidyData$Autophagy <- tidyData$Count >= 6
all_grouping <- group_by(tidyData, Experiment, Condition, Concentration)
summary <- summarise(all_grouping, percentage_autophagy = sum(Autophagy)/length(Autophagy)*100)
final_summary <- ddply(summary, c("Condition", "Concentration"), summarise, MeanAutophagy = mean(percentage_autophagy), SD = sd(percentage_autophagy))

#statistical analysis (Fisher's exact test)
#Counts
summary_counts <- summarise(all_grouping, Proportions = sum(Autophagy), TotalCells = length(Autophagy))
data_fisher <- ddply(summary_counts, c("Condition", "Concentration"), summarise, Autophagy = sum(Proportions), noAutophagy = sum(TotalCells-Proportions))

A10uM <- data_fisher$Autophagy[data_fisher$Concentration =="10uM"]
nA10uM <- data_fisher$noAutophagy[data_fisher$Concentration == "10uM"]
A50uM <- data_fisher$Autophagy[data_fisher$Concentration == "50uM"]
nA50uM <- data_fisher$noAutophagy[data_fisher$Concentration == "50uM"]

fisher10uM <- matrix( c(A10uM, nA10uM), byrow = TRUE, nr = 2)
fisher50uM <- matrix( c(A50uM, nA50uM), byrow = TRUE, nr = 2)

pvalue10uM <- fisher.test(fisher10uM)
pvalue50uM <- fisher.test(fisher50uM)


#visualization of data
tidyPlot <- ggplot(tidyData, aes(x = Condition, y = Count, fill = Concentration)) +
				geom_boxplot() + geom_point() + theme_bw() +
				scale_fill_grey(start = 0, end = .9)


		
dotplot_mean <- ggplot(summary, aes(x = Concentration, y = percentage_autophagy, shape = Condition, col = Condition, width = 0.8)) +
		geom_point(size = 3, position = position_dodge(width=0.8)) + 
		stat_summary(fun.data = mean_sdl, fun.args = list(mult = 0), color = "black", geom = "errorbar", size = 1, position = position_dodge(width = 0.8)) +
		scale_color_grey(start = 0, end = 0.5, name ="Condition") + scale_shape(name = "Condition") +
		xlab("") + ylab ("Autophagic Cells(%)")
		
facetplot_mean <- ggplot(summary, aes(x = Concentration, y = percentage_autophagy, shape = Condition, col = Condition, width = 0.8)) +
		geom_point(size = 5, position = position_dodge(width=0.8)) + 
		stat_summary(fun.data = mean_sdl, fun.args = list(mult = 0), color = "black", geom = "errorbar", size = 1, position = position_dodge(width = 0.8)) +
		scale_color_grey(start = 0, end = 0.5, name ="Condition") + scale_shape(name = "Condition") +
		xlab("") + ylab ("Autophagic Cells(%)")
ggsave(tidyPlot, file = 'Hek_autophagy.pdf', width=4, height=4)
ggsave(dotplot_mean, file = 'dotplotmean.svg', width=4, height=4)
