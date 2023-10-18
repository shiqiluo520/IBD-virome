#alpha-div. Input: tsv.file with abundance of vOTUs or bASVs in each sample. Output: csv.file with alpha diversity for each sample.
abundance <- read.csv('abundance.tsv',row.names = NULL,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
alpha <- function(x, tree = NULL, base = exp(1)) {
  est <- estimateR(x)
  Richness <- est[1, ]
  Chao1 <- est[2, ]
  ACE <- est[4, ]
  Shannon <- diversity(x, index = 'shannon', base = base)
  Simpson <- diversity(x, index = 'simpson')    #Gini-Simpson ??????
  Pielou <- Shannon / log(Richness, base)
  goods_coverage <- 1 - rowSums(x == 1) / rowSums(x)
  result <- data.frame(Richness, Shannon, Simpson, Pielou, Chao1, ACE, goods_coverage)
}
alpha <- alpha(abundance, base = 2)
write.csv(alpha, 'alpha.diversity.csv', quote = FALSE)

#beta-div. Input: tsv.file with abundance of vOTUs or bASVs in each sample and csv file with metadata for each sample. Output: PCoA plot using bray-curtis dissimilarity of vOTUs or bASVs in each sample.
library(vegan)
abundance <- read.csv('abundance.tsv',row.names = NULL,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
convCounts2Abun<-function(otucounttable){
  colsum<-apply(otucounttable,2,sum)
  colmax<-matrix(rep(colsum,each=nrow(otucounttable)),nrow=nrow(otucounttable))
  otu.abundance.table<-otucounttable*100/colmax
}
abundance_per <-convCounts2Abun(abundance)
abundance_per <- t(abundance_per)
group <- read.csv('metadata.csv', sep = ',', stringsAsFactors = FALSE)
distance <- vegdist(abundance_per, method = 'bray')
ID = as.data.frame(row.names(abundance_per))
colnames(ID)[1] <-"ID"
ID_group<- inner_join(ID, group, by = 'ID')
set.seed(1234)
adonis_result <- adonis2(distance~Pathological.state, ID_group, permutations = 999)
pcoa <- cmdscale(distance, k = (nrow(abundance_per) - 1), eig = TRUE)
pcoa_eig <- (pcoa$eig)[1:3] / sum(pcoa$eig)
sample_site <- merge(sample_site, group, by = 'ID', all.x = TRUE)
pcoa_plot <- ggplot(sample_site, aes(PCoA1, PCoA2)) +
  geom_point(aes(color=Strain.genotype), size = 5,stroke =1.5 ,alpha=0.9) +
  scale_shape_manual(values = c(16))+scale_color_manual(values=c( "#606060", "#D0B1B1","#D4D4D4","#91AACF","#B0A87B"))+
  labs(x = paste('PCoA 1: ', round(100 * pcoa_eig[1], 2), '%'), y = paste('PCoA 2: ', round(100 * pcoa_eig[3], 2), '%')) +
  stat_ellipse(aes(linetype= Pathological.state),type = "norm",level=0.85)+
  annotate("text", x=Inf, y=Inf,
           label=paste0("R2 = ", round(adonis_result$R2[1], 3), "\n", "p = ", round(adonis_result$"Pr(>F)"[1], 3)),
           hjust=1, vjust=1, size=6)

ggsave('beta-diversity.png', pcoa_plot, width = 10, height = 6)
