rm(list = ls())

library(ggplot2)
library(hrbrthemes)
library(viridis)
library(ggridges)
library(ggsci)
library(tidyverse)
library(wesanderson)

## Load model comparisons: mechanistic vs random genes
load("./Objs/KTSP/ModelCompare_KTSP_DiffNoFeat.rda")
load("./Objs/RF/ModelCompare_RF_DiffNoFeat_new.rda")
load("./Objs/SVM/ModelCompare_SVM_DiffNoFeat_new.rda")
load("./Objs/XGB/ModelCompare_XGB_DiffNoFeat_new.rda")

#########
ModelCompare_KTSP_DiffNoFeat$algorithm <- "K-TSPs"
ModelCompare_RF_DiffNoFeat$algorithm <- "RF"
ModelCompare_SVM_DiffNoFeat$algorithm <- "SVM"
ModelCompare_XGB_DiffNoFeat$algorithm <- "XGB"

## Bind the 4 together in one data frame
AllModelCompare_Bladder_DiffNoFeat <- rbind(ModelCompare_KTSP_DiffNoFeat, ModelCompare_RF_DiffNoFeat, ModelCompare_SVM_DiffNoFeat, ModelCompare_XGB_DiffNoFeat)

AllModelCompare_Bladder_DiffNoFeat$data_type <- factor(AllModelCompare_Bladder_DiffNoFeat$data_type, levels = c("Training", "Testing"))

AllModelCompare_Bladder_DiffNoFeat$modelType <- factor(AllModelCompare_Bladder_DiffNoFeat$modelType, levels = c("Agnostic_DEGs", "Agnostic_Pairs", "Mechanistic"))

levels(AllModelCompare_Bladder_DiffNoFeat$modelType) <- c("Agnostic DEGs", "Agnostic pairs", "Mechanistic")

############################################################################
## Remove the training data
table(AllModelCompare_Bladder_DiffNoFeat$data_type)

AllModelCompare_Bladder_DiffNoFeat <- AllModelCompare_Bladder_DiffNoFeat %>%
  filter(data_type == "Testing")

AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn <- factor(AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn, levels =  c("37_Pairs", "50_Pairs", "100_Pairs", "250_Pairs", "74_Genes", "100_Genes", "200_Genes", "500_Genes"))
levels(AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn) <- c("37 Pairs", "50 Pairs", "100 Pairs", "250 Pairs", "74 Genes", "100 Genes", "200 Genes", "500 Genes")

sel <- which(AllModelCompare_Bladder_DiffNoFeat$algorithm == "K-TSPs" & AllModelCompare_Bladder_DiffNoFeat$modelType == "Mechanistic" & AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn != "37 Pairs")
AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn[sel] <- "37 Pairs"

sel2 <- which(AllModelCompare_Bladder_DiffNoFeat$modelType == "Mechanistic" & AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn != "37 Pairs")
AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn[sel2] <- "37 Pairs"

table(AllModelCompare_Bladder_DiffNoFeat$modelType, AllModelCompare_Bladder_DiffNoFeat$NofFeatAgn)

AllModelCompare_Bladder_DiffNoFeat$FeatureType <- ifelse(AllModelCompare_Bladder_DiffNoFeat$modelType == "Agnostic DEGs", "Genes", "Pairs")

###########################
############################################################################
## Plot
My_Theme = theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.title.y = element_text(size = 14),
  axis.text.y = element_text(size = 12),
  strip.text.x = element_text(size = 14),
  plot.title = element_text(size=15, face = "bold", hjust = 0.5)
)


## Density plot
png(filename = "./Figs/Bladder_BS_AllModels_Density_DiffNoFeatAndPairs_new.png", width = 3000, height = 1500, res = 200)
BS_AUC_ModelCompare <- ggplot(AllModelCompare_Bladder_DiffNoFeat, aes(x = AUC, y = modelType, fill = NofFeatAgn, height = ..ndensity..)) + 
  geom_density_ridges(stat = "density", alpha = 0.9, bw = 0.8, adjust= 0.01, scale=1.2) +
  scale_x_continuous(limits = c(0.5, 0.9), breaks = seq(0.5, 0.9, by = 0.1), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title="Mechanistic (37 FFLs) vs agnostic (different N of features/pairs) performance at predicting bladder cancer progression in the testing data") + 
  ylab("Model Type") +
  My_Theme +
  #scale_fill_manual(values = c("#21908CFF", "cadetblue1", "#3B9AB2", "cyan4", "#440154FF", "darkcyan", "#21908CFF", "#21908CFF")) +
  #scale_fill_jco() +
  scale_alpha_manual(values = c(0.4, 1))+
  scale_fill_viridis(discrete = TRUE)+
  coord_cartesian(clip = "off") +
  facet_wrap(~algorithm, dir = "v", ncol = 1, scales = "free_y") + theme(axis.text = element_text(face = "bold"), 
                                            panel.background = element_rect(fill = "gray93"), 
                                            plot.background = element_rect(fill = "white")) +labs(fill = "N of features/pairs", alpha = "data type")

BS_AUC_ModelCompare
dev.off()

