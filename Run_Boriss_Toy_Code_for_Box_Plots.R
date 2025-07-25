


#Set your working directory
setwd("/Users/alipka/Library/CloudStorage/Box-Box/Sabbatical_Roslin_Institute/R_workspace/omnigenic_sim")
home.dir <- getwd()

#####Read in all of the packages that are necessary
#Read in prerequiste libaries for GAPIT
library('MASS')
library(gplots)
library(sommer)
library(reshape)
library(ggplot2)
library(cowplot)



#########

correlation.core.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.core.QTNs.RDS") 
correlation.core.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.core.SNPs.RDS")
correlation.peripheral.QTNs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.peripheral.QTNs.RDS") 
correlation.peripheral.SNPs <- readRDS("/Users/alipka/Library/CloudStorage/Box-Box/IR-281/results-analysis-20250421/master.these.spearman.rank.correlation.between.GWAS.peripheral.SNPs.RDS")


#Format the results for trait values, genetic values, and breeding values so that they can be
# used to make graphs
Rep <- c(1:3)
AddC_AddP <- c(1,2,4)
EpiCC_AddC <- c(0.05, 0.5, 1, 2)
EpiPP_AddP <- c(0.05, 0.5, 1, 2)
EpiCP_AddP <- c(0.05, 0.5, 1, 2)
Selection <- c("Directional", "Stabilizing", "Disruptive")
Intensity <- c("10%", "20%")
Genes <- c("Core", "Peripheral")

grid <- expand.grid(Rep, AddC_AddP, EpiCC_AddC, EpiPP_AddP, EpiCP_AddP, 
                    Selection, Intensity, Selection, Intensity, Genes)
set.seed(1000)
# core_corr <- c(seq(-1, -0.4, 0.1), seq(0.4, 1, 0.1))
# core_corr <- c(
results <- data.frame(Rep = grid$Var1,
                      AddC_AddP = grid$Var2, 
                      EpiCC_AddC = grid$Var3, 
                      EpiPP_AddP = grid$Var4, 
                      EpiCP_AddP = grid$Var5,
                      Target_population = paste(grid$Var6, grid$Var7, sep = "_"),
                      Comp_population = paste(grid$Var8, grid$Var9, sep = "_"), 
                      Genes = grid$Var10,
                      Correlation = sample(seq(-1, 1, 0.1), nrow(grid), replace = T))

results <- within(results, {
  Target_population <- factor(Target_population, levels = c("Directional_10%", "Directional_20%", 
                                                            "Stabilizing_10%", "Stabilizing_20%", 
                                                            "Disruptive_10%", "Disruptive_20%"))
  Comp_population <- factor(Comp_population, levels = c("Directional_10%", "Directional_20%", 
                                                        "Stabilizing_10%", "Stabilizing_20%", 
                                                        "Disruptive_10%", "Disruptive_20%"))
  Genes <- factor(Genes, levels = c("Core", "Peripheral"))
})

library(tidyverse)

p1 <- ggplot(results, aes(x = Comp_population, y = Correlation)) +
  geom_boxplot(aes(fill = Target_population, linetype = Genes)) +
  labs(x = "Comparison population", y = "Spearman correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p1

ggsave("Results_summary1.jpg", p1, width = 10, height = 6, dpi = 600)

p2 <- ggplot(results, aes(x = Comp_population, y = Correlation)) +
  geom_boxplot(aes(fill = Target_population)) +
  labs(x = "Comparison population", y = "Spearman correlation") +
  facet_wrap(~Genes, nrow = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2
ggsave("Results_summary2.jpg", p2, width = 10, height = 6, dpi = 600)

p3 <- ggplot(results, aes(x = Comp_population, y = Correlation)) +
  geom_boxplot(aes(fill = Genes)) +
  facet_wrap(~Target_population, nrow = 2) +
  labs(x = "Comparison population", y = "Spearman correlation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p3
ggsave("Results_summary3.jpg", p3, width = 10, height = 6, dpi = 600)




