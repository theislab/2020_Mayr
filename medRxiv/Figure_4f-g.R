## panel 4f random forest lung function prediction of single cell cohorts##

## files needed:
## BALF_proteome_annot_imput.txt
## categorical annotations_Gauting2.txt
## Numerical annotations_Gauting2.txt
## pseudocounts_schiller.csv
## pseudocounts_reyfman.csv
## pseudocounts_banovich.csv
## metacelltypes_health-status_mastertable.csv

# Load R libs ####
library(Matrix)
library(randomForest)
library(preprocessCore)
library(ggplot2)
library(pheatmap)
cm_colors=colorRampPalette(c("#253494","white","red"))(100)
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files")

# Load BALF data ####
balf <- read.delim("BALF_proteome_annot_imput.txt")
balf_expr <- balf[39:nrow(balf), 1:124]
balf_expr <- data.matrix(balf_expr)
balf_genes <- as.character(balf[39:nrow(balf), 'Gene.names'])

# Quantile normalize data ####
tmp <- normalize.quantiles(balf_expr)
colnames(tmp) <- colnames(balf_expr)
balf_quant_normalized <- tmp

# Load annotations ####
#setwd('/Users/lsimon/OneDrive/Miko/Helmholtz/Schiller/data/Clinical study/')
cat <- read.table('categorical annotations_Gauting2.txt', header = T, sep = '\t', row.names = 1)
num <- read.table('Numerical annotations_Gauting2.txt', header = T, sep = '\t', row.names = 1)
cat <- cat[colnames(balf_expr),]
num <- num[colnames(balf_expr),]

# Correlate protein expression with lung function ####
col_nom <- "Dlco..SB.....Soll."
lung_function <- num[, col_nom]
correls <- apply(balf_expr, 1, function(x) cor(method = "spearman", x, use = "complete.obs", lung_function))

# Load pseudobulks ####
sums_schiller <- data.matrix(read.csv("pseudocounts_schiller.csv", row.names = 1))
sums_reyfman <- data.matrix(read.csv("pseudocounts_reyfman.csv", row.names = 1))
sums_banovich <- data.matrix(read.csv("pseudocounts_banovich.csv", row.names = 1))
treat_schiller <- scan(what = character(), sep = "\n", "treat_ild_schiller.txt")
treat_reyfman <- scan(what = character(), sep = "\n", "treat_ild_reyfman.txt")
treat_banovich <- scan(what = character(), sep = "\n", "treat_ild_banovich.txt")

# Define functions ####
normalize_pseudobulks <- function(sums){
  insilico <- normalize.quantiles(sums)
  insilico <- t(apply(insilico, 1, function(x) (x - mean(x))/sd(x)))
  rownames(insilico) <- rownames(sums)
  colnames(insilico) <- colnames(sums)
  insilico <- t(apply(insilico, 1, function(x) (x - mean(x))/sd(x)))
  insilico
}

# Define set of genes present in proteome and all pseudobulks ####
present_in_all <- intersect(rownames(sums_schiller), intersect(rownames(sums_banovich), rownames(sums_reyfman)))
ok <- intersect(present_in_all, balf_genes)
correls_ok <- correls[match(ok, balf_genes)]
balf_ok <- balf_quant_normalized[match(ok, balf_genes),]
names(correls_ok) <- rownames(balf_ok) <- ok

# Train random forest on protein data to predict lung function (Dlco) ####
good <- names(which(abs(correls_ok) > 0.2))
tmp <- balf_ok[good,]
tmp <- tmp[which(apply(tmp, 1, var) > 0),]
tmp <- t(apply(tmp, 1, function(x) (x - mean(x))/sd(x)))
aframe <- data.frame(lung_function, t(tmp))
aframe <- aframe[-which(is.na(aframe$lung_function)),]
rf <- randomForest(lung_function ~., aframe)
plot(aframe$lung_function, predict(rf), xlab = "Observed", ylab = "predicted", main = "Lung function")

# Test on pseudobulks ####
predict_pseudobulks <- function(insilico){
  tmp <- t(insilico[good,])
  tmp <- data.frame(tmp)
  predicted_lung_function <- predict(rf, tmp)
  predicted_lung_function
}
preds_schiller <- predict_pseudobulks(insilico = normalize_pseudobulks(sums_schiller))
preds_reyfman <- predict_pseudobulks(insilico = normalize_pseudobulks(sums_reyfman))
preds_banovich <- predict_pseudobulks(insilico = normalize_pseudobulks(sums_banovich))

# Combine results into one table ####
aframe <- data.frame(prediction = c(preds_schiller, preds_reyfman, preds_banovich),
                     treat = c(treat_schiller, treat_reyfman, treat_banovich),
                     study = c(rep("schiller", length(preds_schiller)), rep("reyfman", length(preds_reyfman)), rep("banovich", length(preds_banovich))))
aframe$treat[which(aframe$treat == "TRUE")] <- "IPF"
aframe$treat[which(aframe$treat == "FALSE")] <- "HC"
aframe$treat[which(aframe$treat == "1")] <- "HC"
aframe$treat[which(aframe$treat == "2")] <- "IPF"
aframe$treat <- as.character(aframe$treat)
aframe$group <- paste(aframe$study, aframe$treat)

# Run statistics ####
summary(aov(prediction ~ study + treat, aframe))

# Visualize predictions ####
ggplot(aframe, aes(group, prediction, fill = treat)) +
  geom_boxplot(outlier.size = NULL) + geom_jitter() +
  ylab("Predicted lung function [Dlco]")




## panel 4g show example genes and their gene expression that drive random forest lung function prediction##

# Load foldchanges from ILD scRNAseq differential expression analysis ####
megatable <- read.csv('metacelltypes_health-status_mastertable.csv', row.names = 1)
megatable$Neu.t.value <- NULL
megatable$Neu.p.value <- NULL
megatable$Neu.avg.detection <- NULL
megatable$NA.t.value <- NULL
megatable$NA.p.value <- NULL
megatable$NA.avg.detection <- NULL
megatable$DC.EREG.t.value <- NULL
megatable$DC.EREG.p.value <- NULL
megatable$DC.EREG.avg.detection <- NULL

# Get important genes from random forest ####
tmp <- rf$importance
tmp <- tmp[sort.list(-tmp[,1]),]
good <- names(head(tmp,30))

# Show cell-type specific foldchanges of imporant genes ####
ok <- intersect(good, rownames(megatable))
fcs <- grep('t.value', fixed = T, colnames(megatable))
fcs <- data.matrix(megatable[, fcs])
colnames(fcs) <- gsub('.t.value', '', fixed = T, colnames(fcs))
tmp <- fcs[ok,]
tmp[which(is.na(tmp))] <- 0
pheatmap(tmp, breaks = seq(-4, 4, length = 101), col = colorRampPalette(c("#253494","white","red"))(101), cluster_cols=T)
