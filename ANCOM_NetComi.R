#ANCOM-BC. Input: pseq. file with vOTUs, taxonomy and metadata. Output: volcano plot showing the differential abundant vOTUs in inflamed and non-inflamed mice.
library(ANCOMBC)
#change chater to factor in a phyloseq
sample_data(pseq)$Inflammation <- as.factor(sample_data(pseq)$Inflammation)
set.seed(123)
DA = ancombc2(data = pseq,assay_name = "counts", tax_level = NULL,fix_formula = "Inflammation", rand_formula = NULL,
                 p_adj_method = "BH", prv_cut = 0, lib_cut = 0,
                 group = NULL, struc_zero = FALSE, neg_lb = FALSE,alpha = 0.05, n_cl = 2, verbose = TRUE,
                 global = FALSE, pairwise = FALSE,
                 dunnet = FALSE, trend = FALSE,
                 iter_control = list(tol = 1e-5, max_iter = 20,
                                     verbose = FALSE),
                 em_control = list(tol = 1e-5, max_iter = 100),
                 lme_control = NULL, mdfdr_control = NULL,
                 trend_control = NULL)
ancom_res_df <- data.frame(
  Contig = unlist(DA$res$taxon),
  logFoldChange = unlist(DA$res$lfc_Inflammation),
  q_val = unlist(DA$res$q_Inflammation))
ancom_res_df$level = as.factor(ifelse(ancom_res_df_BH$logFoldChange >2 & ancom_res_df_BH$q_val < 0.05, "enriched",ifelse(ancom_res_df_BH$logFoldChange < -2 & ancom_res_df_BH$q_val < 0.05, "depleted", "nosig")))
p <- ggplot(ancom_res_df, aes(x=logFoldChange, y=-log2(q_val))) + 
  geom_point(aes(color=level))  +
  scale_colour_manual(values=c("#FFCC32","#6500CC","#AEAEAE"))
ggsave(file=paste("DA_vOTUs.png", sep=""), p, width = 24, height = 15)

#Network. Input: csv.file with abundance of differential abundant vOTUs and bASVs in inflamed or non-inflamed mice. Output: Network plot showing the differential abundant vOTUs in inflamed and non-inflamed mice with bacterial genera.
library(NetCoMi)
library(corrplot)
vOTU_bASV <- read.csv('vOTUs_bASV.csv',row.names = NULL,sep = ',', stringsAsFactors = FALSE, check.names = FALSE)
set.seed(123456)
assoMat_vOTU_bASV <- cor(vOTU_bASV,method="spearman")
net_vir_16S <- netConstruct(data = assoMat_vOTU_bASV,
                                     dataType = "condDependence", 
                                     sparsMethod = "none")
netprops_vir_16S<- netAnalyze(net_vir_16S, clustMethod = "cluster_fast_greedy",hubPar = c("degree"),hubQuant = 0.85)
pdf('Output/network.pdf', width = 25, height = 25)
plot(netprops_vir_16S,
     layout = "spring",
     repulsion = 1,
     labelScale = TRUE,
     rmSingles = "all",
     nodeSize = "degree",
     nodeSizeSpread = 2,
     labelFont = 4,
     nodeColor = "feature", 
     colorVec = colorsetting,
     borderWidth = 1,
     hubBorderCol = "darkgray",
     highlightHubs = TRUE,
     hubBorderWidth = 1,
     cexHubs = 1,
     cexNodes = 1,
     cexLabels = 0,
     cexHubLabels	= 5, #hub node label
     nodeTransp = 0,
     shortenLabels = "none",
     hubLabelFont =2,
     edgeWidth = 1)
dev.off()
