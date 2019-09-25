data<-read.table('peptide_toxicity.txt', header=T)
library(plyr)
library(ggplot2)

data['Total'] <- data$Alive + data$Dead
data['Percentage.alive'] <- data$Alive * 100/data$Total
data['Percentage.dead'] <- 100.00 -data$Percentage.alive

data_summary <- ddply(data, "Condition", summarise, mean.alive = mean(Percentage.alive),
        mean.dead = mean(Percentage.dead))

surv_boxplot <- boxplot(data$'Percentage.alive'~data$'Condition', data, main='Survival Percentage', xlab='Condition',
        ylab='Percentage', col=c('lightskyblue1', 'lightskyblue2', 'lightskyblue', 'lightskyblue3', 'steelblue'))

surv_plot <- ggplot(data, aes(x = Condition, y = Percentage.alive, fill = Condition)) +
				geom_boxplot() + geom_point() + theme_bw() +
				ylim(0, 100) + scale_fill_grey(start = 0, end = .9)
				

ggsave(surv_plot, file = 'MEFs_survival.pdf')

aov_cond <- aov(data$Percentage.alive ~ data$Condition)
print(summary(aov_cond))

tuk <- TukeyHSD(aov_cond)
print(tuk)
