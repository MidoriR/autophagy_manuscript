library(grid)
library(lattice)
library(gridExtra)
library(ggplot2)

#Load tables
cytotoxicity=read.table('citotoxicidadIP1-Macrofagos4R_sinT0.csv', head=TRUE, sep=',')
autophagy=read.table('datos_cytoid_jacqueline_boxplot4R.csv', head=TRUE, sep=',')

#Generate plots
p1 <- ggplot(cytotoxicity, aes(x = cytotoxicity$Condition, y = cytotoxicity$Value)) + geom_boxplot() + theme_bw() + 
  scale_x_discrete(name = 'Condition') + scale_y_continuous(name = 'Cell survival (%)') + theme(axis.text.x=element_text(angle=90,hjust=1))

p2 <- ggplot(autophagy, aes(x = autophagy$Condition, y = autophagy$Value)) + geom_boxplot() + theme_bw() +
  scale_x_discrete(name = 'Condition') + scale_y_continuous(name = 'Autophagy (%)') + theme(axis.text.x=element_text(angle=90,hjust=1))


#Save composed image
pdf('fig13.pdf', width = 8, height = 12)
grid.arrange(p1,p2,ncol=1,nrow=2)
dev.off()

#Perform anova test
aov_cytotoxicity <- aov(lm(cytotoxicity$Value ~ cytotoxicity$Condition))
aov_autophagy <- aov(lm(autophagy$Value ~ autophagy$Condition))

#Calculate p values from aov in a Post-hoc analysis with Tukey test

tukey_cytotoxicity <- TukeyHSD(aov_cytotoxicity)
tukey_autophagy <- TukeyHSD(aov_autophagy)

print('P values for cytotoxicity:')
print(tukey_cytotoxicity)

print('P values for autophagy:')
print(tukey_autophagy)

#Save output to a txt file
capture.output(tukey_cytotoxicity, file='Cytotoxicity.txt')
capture.output(tukey_autophagy, file='Autophagy.txt')