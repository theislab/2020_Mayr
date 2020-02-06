## Figure 4: Protein signatures in ELF predict lung function decline and the corresponding cellular changes

## panel 4a: Diagnosis specific cell type contribution

# files needed:
# BALF_proteome_z-score_diagnosis-groups.xlsx
# SchiReyBan_final_allMarkers_cell_type.xlsx

# This code will perform deconvolution

# Load R libs ####
library(readxl)
library(pheatmap)
library(BiocManager)
library(preprocessCore)
library(dplyr)
library(plotrix)
cm_colors=colorRampPalette(c("#253494","white","red"))(100)
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files")

# Load clinical metainfo correlations ####
prot <- read_excel("BALF_proteome_z-score_diagnosis-groups.xlsx", sheet = 1)
prot_correlations <- prot[,grep("z-score", colnames(prot), fixed = T)]

# Load FindAllMarkers table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)
markers$log_FC <- as.numeric(markers$log_FC)

# Define cell type signatures ####
tmp <- markers[which(markers$p_val_adj < 0.1 & markers$log_FC > 3),]
signatures <- split(tmp$gene, as.character(tmp$cluster))

# Run deconvolution ####
ks_results <- lapply(names(signatures), function(sig){
  print(sig)
  sig <- signatures[[sig]]
  ok <- intersect(unlist(signatures), prot$`Gene names`)
  correl <- data.matrix(prot_correlations[match(ok, prot$`Gene names`),])
  rownames(correl) <- ok
  res <- do.call(rbind, lapply(colnames(correl), function(col){
    col <- correl[, col]
    ok <- intersect(sig, rownames(correl))
    if(length(ok) <=5) return(c(NA, NA))
    set <- col[match(ok, rownames(correl))]
    rest <- col[-match(ok, rownames(correl))]
    ks <- ks.test(set, rest)
    c(ks$p.value, mean(set) - mean(rest))
  }))
})

signed_pvals <- do.call(cbind, lapply(ks_results, function(x){
  pvals <- (-log10(x[,1]))
  pvals[which(x[,2] < 0)] <- pvals[which(x[,2] < 0)]*(-1)
  pvals
}))
colnames(signed_pvals) <- names(signatures)
rownames(signed_pvals) <- gsub("z-score", "", fixed = T, colnames(prot_correlations))

# Generate heatmap ####
tmp <- signed_pvals
tmp[which(is.na(tmp))] <- 0
tmp <- tmp[,which(apply(tmp, 2, var) > 0)]
tmp[which(tmp > 4)] <- 4
tmp[which(tmp < (-2))] <- (-2)
paletteLength=100
mybreaks <- c(seq(min(tmp), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(tmp)/paletteLength, max(tmp), length.out=floor(paletteLength/2)))
pheatmap(tmp, breaks = mybreaks, color= cm_colors)




## panel 2b: Clinical paramter specific cell type contribution

# files needed:
# Table S2.xlsx
# SchiReyBan_final_allMarkers_cell_type.xlsx


# This code will perform deconvolution
cm_colors=colorRampPalette(c("#253494","white","red"))(80)
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files")

# Load clinical metainfo correlations ####
prot <- read_excel("Table S2.xlsx", sheet = 1)
prot_correlations <- prot[,grep("Pearson correlation Corr.Coef.", colnames(prot), fixed = T)]

# Load FindAllMarkers table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)
markers$log_FC <- as.numeric(markers$log_FC)

# Define cell type signatures ####CM: Take top100 per celly type###
tmp2 <- markers[order(-markers$log_FC),]
tmp2 <-group_by(tmp2,cluster)
tmp <-top_n(tmp2, 175, log_FC)
signatures <- split(tmp$gene, as.character(tmp$cluster))

# Define cell type signatures ####
#tmp <- markers[which(markers$p_val_adj < 0.1 & markers$log_FC > 2),]
#signatures <- split(tmp$gene, as.character(tmp$cluster))

# Run deconvolution ####
ks_results <- lapply(names(signatures), function(sig){
  print(sig)
  sig <- signatures[[sig]]
  ok <- intersect(unlist(signatures), prot$`Gene names`)
  correl <- data.matrix(prot_correlations[match(ok, prot$`Gene names`),])
  rownames(correl) <- ok
  res <- do.call(rbind, lapply(colnames(correl), function(col){
    col <- correl[, col]
    ok <- intersect(sig, rownames(correl))
    if(length(ok) <=5) return(c(NA, NA))
    set <- col[match(ok, rownames(correl))]
    rest <- col[-match(ok, rownames(correl))]
    ks <- ks.test(set, rest)
    c(ks$p.value, mean(set) - mean(rest))
  }))
})

signed_pvals <- do.call(cbind, lapply(ks_results, function(x){
  pvals <- (-log10(x[,1]))
  pvals[which(x[,2] < 0)] <- pvals[which(x[,2] < 0)]*(-1)
  pvals
}))
colnames(signed_pvals) <- names(signatures)
rownames(signed_pvals) <- gsub("Pearson correlation Corr.Coef._", "", fixed = T, colnames(prot_correlations))

# Generate heatmap ####
tmp <- signed_pvals
tmp[which(is.na(tmp))] <- 0
tmp <- tmp[,which(apply(tmp, 2, var) > 0)]
tmp[which(tmp > 4)] <- 4
tmp[which(tmp < (-4))] <- (-4)
pheatmap(tmp, color= cm_colors)




## panel 2c,d: Example plots for correlations with cell types and the genes driving the correlation
# files needed:
# Table S2.xlsx
# SchiReyBan_final_allMarkers_cell_type.xlsx

# Load R libs ####
library(readxl)
library(pheatmap)
library(BiocManager)
library(preprocessCore)
library(dplyr)
library(plotrix)
cm_colors=colorRampPalette(c("#253494","white","red"))(100)
# Define function to generate ECDF plots ####
# Load clinical metainfo correlations ####
prot <- read_excel("Table S2.xlsx", sheet = 1)
prot_correlations <- prot[,grep("Pearson correlation Corr.Coef.", colnames(prot), fixed = T)]

# Load FindAllMarkers table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)
markers$log_FC <- as.numeric(markers$log_FC)

# Define cell type signatures ####CM: Take top100 per celly type###
tmp2 <- markers[order(-markers$log_FC),]
tmp2 <-group_by(tmp2,cluster)
tmp <-top_n(tmp2, 175, log_FC)
signatures <- split(tmp$gene, as.character(tmp$cluster))
par(mfrow = c(1,1))
genECDFplot <- function(celltype, param){
  genes <- intersect(signatures[[celltype]], as.character(prot$`Gene names`))
  print(genes)
  background <- setdiff(prot$`Gene names`, genes)
  plot(ecdf(unlist(prot[match(background, prot$`Gene names`), param])),
       xlab = param, main = celltype, ylab = 'Percentile')
  lines(ecdf(unlist(prot[match(genes, prot$`Gene names`), param])), col = ' red')
  abline(v = 0, lty = 2)
  legend('bottomright', c(celltype, 'Background'), col = c('red', 'black'), bty = 'n', pch = 19)
}
genECDFplot(celltype = "Myofibroblasts", param = "Pearson correlation Corr.Coef._DLCO (SB) [Hb corrected % Soll]")
genECDFplot(celltype = "Plasma cells", param = "Pearson correlation Corr.Coef._Zytospin: BAL Alveolarmacrophages")




