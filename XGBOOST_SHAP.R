#XGBOOST Input: csv.file with abundance and metadata each sample. Output: Plot of SHAP value of each factors to inflammation.
library(SHAPforxgboost)
data <- read.csv("abundance_metadata_vOTUs.csv", row.names = NULL,sep=",")
data <- mutate(data, Inflammation = as.factor(Inflammation))
y <- "Inflammation"
x <-c("Virome.Chao1","Virome.Shannon","Bacteriome.Shannon","Bacteriome.Chao1","Facility","Animal.strain","Genotype")
set.seed(9999)
IL10_sample  <- sample(nrow(IL10), 0.8 * nrow(IL10))
IL10_train  <- xgb.DMatrix(data.matrix(IL10[IL10_sample, x]),
                           label = IL10[IL10_sample, y])
IL10_valid  <- xgb.DMatrix(data.matrix(IL10[-IL10_sample, x]),
                           label = IL10[-IL10_sample, y])
params <- list(
  objective = "reg:squarederror",
  learning_rate = 0.05,
  subsample = 0.9,
  colsample_bynode = 1,
  reg_lambda = 2,
  max_depth = 5
)
fit_xgb <- xgb.train(
  params,
  data = IL10_train,
  watchlist = list(valid = IL10_valid),
  early_stopping_rounds = 20,
  print_every_n = 100,
  nrounds = 10000
)

#SHAP
library(xgboost)
#Compact SHAP analysis
X <- data.matrix(IL10[sample(nrow(IL10), 30), x])
shap <- shap.prep(fit_xgb, X_train = X)
p <- shap.plot.summary(shap)
png("SHAP.png", width = 6, height = 2.2, units = "in", res = 300) #Figure 1B
p
dev.off() 

p1 <- shap.plot.dependence(size0 = 4, alpha=0.5,data_long = shap, x = 'Virome.Chao1', y = 'Virome.Chao1', color_feature = 'Inflammation') + ggtitle("Virome.Chao1 to inflammation")+theme(plot.title = element_text(size = 12, face = "bold"))
png("SHAP-virome.Chao1.png", width = 6, height = 2.2, units = "in", res = 300) #Fig.S4A
p1
dev.off()

