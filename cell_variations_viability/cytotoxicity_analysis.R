datos<-read.table('cytotoxicity_data.txt', header=T)
library(plyr)
library(ggplot2)

datos['Total'] <- datos$Alive + datos$Dead
datos['Percentage.alive'] <- datos$Alive * 100/datos$Total
datos['Percentage.dead'] <- 100.00 -datos$Percentage.alive

data_summary <- ddply(datos, "Condition", summarise, mean.alive = mean(Percentage.alive),
        mean.dead = mean(Percentage.dead))

variations_plot <- ggplot(datos, aes(x = Condition, y = Percentage.alive, fill = Condition)) +
				geom_boxplot() + geom_point(aes(y = Percentage.alive, group = Condition), position = position_dodge(width = 0.75)) +
				facet_wrap(~ seed, scales = 'free') + theme_bw() + ylim(50, 100) + scale_fill_grey(start = 0, end = .9)
				


ggsave(variations_plot, file = 'variations_cytotoxicity.pdf')

# Anova over two variables (condition and number of seeded cells)

two_aov <- aov(datos$Percentage.alive ~ datos$Condition + datos$seed)

# Post analysis of the effects over all combinations of variables

two_tuk <- TukeyHSD(two_aov)

print(two_tuk)

# Anova over one variable (two variables combined in one string)

datos$combined <- paste(datos$Condition, datos$seed, sep="_")

aov_comb <- aov(datos$Percentage.alive ~ datos$combined)

tuk_comb <- TukeyHSD(aov_comb)
print(tuk_comb)
