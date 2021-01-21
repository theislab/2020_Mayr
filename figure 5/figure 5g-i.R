## Figure 5Protein signatures in BALF predict lung function decline and the corresponding cellular changes. (
## panel 5g-i random forest lung function prediction of single cell cohorts and bulk cohorts##

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


##
## RF on published bulk RNAseq of IPF samples from different histopathological stages (GEO GSE124685)
##
# Load early stage IPF bulk samples ####
fpkm <- read.delim(gzfile("GSE124685_Kaminski_rnaseq/GSE124685_farida_a_FPKM_Matix_84_samples_txt.txt.gz"))
sra <- read.delim("GSE124685_Kaminski_rnaseq/GSE124685_SraRunTable.txt")
match <- read.delim("GSE124685_Kaminski_rnaseq/GSE124685_GSE_sample_matching.txt", header = F)
sample_info <- read.xlsx("GSE124685_Kaminski_rnaseq/Data for HS.xlsx", sheetIndex = 1)
merged <- merge(sra, match, by.x = "Sample_Name", by.y = "V1")
treat <- sample_info$Disease.Severity

# Restrict to unique gene symbols ####
genes <- fpkm$gene_id
duplicates <- names(which(table(genes) > 1))
fpkm <- fpkm[-which(genes %in% duplicates | is.na(genes)),]

# Merge samples ####
expr <- data.matrix(fpkm[, -c(1:9)])
expr_ok <- expr
rownames(expr_ok) <- fpkm$tracking_id
tmp <- normalize.quantiles(expr_ok)
colnames(tmp) <- colnames(expr_ok)
rownames(tmp) <- rownames(expr_ok)
expr_ok <- t(apply(tmp, 1, function(x) (x - mean(x))/sd(x)))

# Define set of genes present in proteome and all pseudobulks ####
ok <- intersect(rownames(expr_ok), balf_genes)
correls_ok <- correls[match(ok, balf_genes)]
balf_ok <- balf_quant_normalized[match(ok, balf_genes),]
names(correls_ok) <- rownames(balf_ok) <- ok

# Train random forest on protein data to predict lung function (Dlco) ####
good <- names(which(abs(correls_ok) > 0.2))
tmp <- balf_ok[good,]
tmp <- tmp[which(apply(tmp, 1, var) > 0),]
tmp <- t(apply(tmp, 1, function(x) (x - mean(x))/sd(x)))
lung_function <- num[, col_nom]
aframe <- data.frame(lung_function, t(tmp))
aframe <- aframe[-which(is.na(aframe$lung_function)),]
rf <- randomForest(lung_function ~., aframe)

test <- data.frame(t(expr_ok[intersect(colnames(aframe), rownames(expr_ok)),]))
predicted_scores <- predict(rf, test)

boxplot(split(predicted_scores, treat),
        main = 'Kaminski RNAseq', ylab = 'Predicted lung function')

aframe <- data.frame(prediction = predicted_scores, treat)
ggplot(aframe, aes(treat, prediction, fill = treat)) +
  geom_boxplot(outlier.size = NULL) + geom_jitter() +
  ylab("Predicted lung function [Dlco]")

