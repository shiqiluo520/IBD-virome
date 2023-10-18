#lme. Input: csv file with Shannon diversity and richness of virome and bacteriome and their metadata for each sample. Output: LM plot showing Shannon diversity or richness between virome and bacteriome.
library(nlme)
data <- read.csv('virome_bacteriome_alpahdiv_mata.csv',row.names = 1,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
B.Chao1_V.Chao1 <- lme(Bacteriome.Chao1 ~ Virome.Chao1, random =~1|Facility/Animal.strain/Genotype, data = data)
LMPlot <- ggplot(B.Chao1_V.Chao1, aes(x = Virome.Chao1, y = Bacteriome.Chao1)) +
  stat_smooth(method="lm", se=TRUE, span=0.95, level=0.95,alpha = 0.1, colour = "black") +
  geom_point(color = category,shape = 16, size = 12 ) +
  theme_bw()+
  labs(x = paste('Fixed effect: ', 'Bacteriome.Chao1'), y = paste('Response variable: ', 'Virome.Chao1'))

ggsave('B.Chao1_V.Chao1.png', LMPlot, width = 11, height = 10.5)


#corrplot.  Input: csv file with Shannon diversity and richness of virome and bacteriome for each sample. Output: cor plot showing Shannon diversity or richness between virome and bacteriome.
library(corrplot)
data <- read.csv('virome_bacteriome_alpahdiv.csv',row.names = 1,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
cor_correlation <- cor(data,method="spearman")
res1 <- cor.mtest(data, conf.level = .95)
p.mat = res1$p
p <- corrplot(cor_correlation, tl.col = "black", tl.srt = 45, bg = "White",tl.cex = 1.2,
         title = "corrplot",mar=c(0,0,1,0),
         type = "lower",p.mat = res1$p, sig.level = c(.001, .01, .05),outline="white",
         insig = "label_sig",pch.cex = 1.2, pch.col = "white")
ggsave('corrplot.png', p, width = 11, height = 10.5)
