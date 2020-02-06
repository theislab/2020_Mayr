## Figure 5: Smoking
## panel 5l: Cell type contribution to active/previous, active/never, previous/never smokers fold changes #

# files needed:
# BALF_smoking.xlsx
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
prot <- read_excel("BALF_smoking.xlsx",sheet =1)
prot_correlations <- prot[,grep("Student's T-test Difference", colnames(prot), fixed = T)]

# Load DGE table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)
markers$log_FC <- as.numeric(markers$log_FC)

# Define cell type signatures ####
tmp <- markers[which(markers$p_val_adj < 0.1 & markers$log_FC > 4),]
signatures <- split(tmp$gene, as.character(tmp$cluster))

# Run deconvolution ####
ks_results <- lapply(names(signatures), function(sig){
  sig <- signatures[[sig]]
  ok <- intersect(unlist(signatures), prot$`Gene names`)
  print(ok)
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
rownames(signed_pvals) <- gsub("Student's T-test Difference", "", fixed = T, colnames(prot_correlations))

# Generate heatmap ####
tmp <- signed_pvals
tmp[which(is.na(tmp))] <- (0)
tmp <- tmp[,which(apply(tmp, 2, var) > 0)]
tmp[which(tmp > 4)] <- (4)
tmp[which(tmp < (-4))] <- (-4)
paletteLength = 100
mybreaks <- c(seq(min(tmp), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(tmp)/paletteLength, max(tmp), length.out=floor(paletteLength/2)))
pheatmap(tmp, color=cm_colors,breaks = mybreaks)




## panel 5m: Top features that drive cell type correaltions of active/never smokers #

# files needed:
# metacelltypes_smoking_mastertable.csv
# BALF_smoking_yes-no.xlsx

# Define cell type signatures ####
smoking <- read_excel("BALF_smoking_yes-no.xlsx",sheet =1)
good <- smoking[which(smoking$`Student's T-test Difference y_n` >= 1.5 | smoking$`Student's T-test Difference y_n` <=-1.5 ),]
good <- good$`Gene names`


# Load foldchanges from ILD scRNAseq differential expression analysis ####
megatable <- read.csv('metacelltypes_smoking_mastertable.csv', row.names = 1)
megatable$Neu.t.value <- NULL
megatable$Neu.p.value <- NULL
megatable$Neu.avg.detection <- NULL
megatable$NA.t.value <- NULL
megatable$NA.p.value <- NULL
megatable$NA.avg.detection <- NULL
megatable$DC.EREG.t.value <- NULL
megatable$DC.EREG.p.value <- NULL
megatable$DC.EREG.avg.detection <- NULL

# Show cell-type specific foldchanges of imporant genes ####
ok <- intersect(good, rownames(megatable))
fcs <- grep('t.value', fixed = T, colnames(megatable))
fcs <- data.matrix(megatable[, fcs])
colnames(fcs) <- gsub('.t.value', '', fixed = T, colnames(fcs))
tmp <- fcs[ok,]
tmp[which(is.na(tmp))] <- 0
pheatmap(tmp, breaks = seq(-4, 4, length = 101), col = colorRampPalette(c("#253494","white","red"))(101), cluster_cols=T)
