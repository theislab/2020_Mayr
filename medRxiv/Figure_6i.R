## Figure 6: Plasma
## panel 6i: Plasma patients sorted in lung function categories were deconvoluted into cell types##

# files needed:
# plasma_groﬂhadern+HA-healthy_z-score_FVCquantiles_2
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
prot <- read_excel("plasma_groﬂhadern+HA-healthy_z-score_FVCquantiles_2.xlsx", sheet = 1)
prot_correlations <- prot[,grep("quantile", colnames(prot), fixed = T)]


# Load FindAllMarkers table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)
markers$log_FC <- as.numeric(markers$log_FC)

# Define cell type signatures ####
tmp <- markers[which(markers$p_val_adj < 0.1 & markers$log_FC > 2),]
signatures <- split(tmp$gene, as.character(tmp$cluster))

# Run deconvolution ####
ks_results <- lapply(names(signatures), function(sig){
  print(sig)
  sig <- signatures[[sig]]
  ok <- intersect(unlist(signatures), prot$PG.Genes)
  correl <- data.matrix(prot_correlations[match(ok, prot$PG.Genes),])
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
rownames(signed_pvals) <- gsub("quantile", "FVC ", fixed = T, colnames(prot_correlations))

# Generate heatmap ####
tmp <- signed_pvals
tmp[which(is.na(tmp))] <- 0
tmp <- tmp[,which(apply(tmp, 2, var) > 0)]
tmp[which(tmp > 1.5)] <- 1.5
tmp[which(tmp < (-1))] <- (-1.5)
paletteLength=100
mybreaks <- c(seq(min(tmp), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(tmp)/paletteLength, max(tmp), length.out=floor(paletteLength/2)))
pheatmap(tmp, breaks = mybreaks, color= cm_colors)


