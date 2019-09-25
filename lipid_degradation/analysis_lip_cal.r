library(tidyr)
library(ggplot2)

#Import data

datos<-read.csv('FreeFattyAcids.csv', header=TRUE)
datos$X <- NULL

#reshape data for plotting (wide to long df)

datos_long<- gather(datos, FattyAcid, Concentration, X12.0:Total)

#Create pdf figure

pdf(file='fa.pdf')
sp<- ggplot(datos_long, aes(x=Sample, y=Concentration))
sp+geom_point(aes(color=Sample, shape=Sample))+facet_wrap(~FattyAcid, ncol=7)+theme_bw()+
scale_color_brewer()

dev.off()

#Calculate p values for comparisons within fatty acid chain length of all groups
p_vals <- datos_long %>% group_by(FattyAcid) %>% summarize(pval=anova(lm(Concentration ~ Sample))[,5][1])

#Adjust p values for multiple comparisons
p_vals$adjusted_pval <- (p_vals$pval, 'fdr')

#Write result to table
write.table(p_vals, 'p_values.csv', sep=',')

print(p_vals)