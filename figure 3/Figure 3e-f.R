## Figure 3: Multi-cohort single cell data reveals transcriptional changes in >40 cell types and altered cell type frequencies in diseases progression

## panel 3e The heatmap shows changes of our cell type signatures in published bulk RNA-seq data (GEO
## GSE124685) across different histopathological stages that represent increasing extent of fibrosis from stage 1-3, as
## determined by quantitative micro-CT imaging and tissue histology. 

# This code will perform deconvolution of bulk RNA-seq data (GEOGSE124685)
# files needed:
# GSE124685_IPF1_vs_Ctrl.xlsx
# GSE124685_IPF2_vs_Ctrl.xlsx
# GSE124685_IPF3_vs_Ctrl.xlsx
# SchiReyBan_final_allMarkers_cell_type.xlsx

# Load R libs ####
library(readxl)
library(pheatmap)
library(ggplot2)
library(cowplot)
library(ggrepel)
library(gridExtra)
library(dplyr)
theme_set(theme_cowplot())
cm_colors=colorRampPalette(c("#253494","white","red"))(100)
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files")

# Load foldchanges from previously published data ####
res_ipf1 <- read_excel("GSE124685_IPF1_vs_Ctrl.xlsx")
res_ipf2 <- read_excel("GSE124685_IPF2_vs_Ctrl.xlsx")
res_ipf3 <- read_excel("GSE124685_IPF3_vs_Ctrl.xlsx")

colnames(res_ipf1)[1] <- colnames(res_ipf2)[1] <- colnames(res_ipf3)[1] <- "GENE_SYMBOL"

fcs_ipf1 <- res_ipf1$logFC
names(fcs_ipf1) <- res_ipf1$GENE_SYMBOL
fcs_ipf1[which(fcs_ipf1 > 10)] <- 10
fcs_ipf1[which(fcs_ipf1 < (-10))] <- (-10)

fcs_ipf2 <- res_ipf2$logFC
names(fcs_ipf2) <- res_ipf2$GENE_SYMBOL
fcs_ipf2[which(fcs_ipf2 > 10)] <- 10
fcs_ipf2[which(fcs_ipf2 < (-10))] <- (-10)

fcs_ipf3 <- res_ipf3$logFC
names(fcs_ipf3) <- res_ipf3$GENE_SYMBOL
fcs_ipf3[which(fcs_ipf3 > 10)] <- 10
fcs_ipf3[which(fcs_ipf3 < (-10))] <- (-10)

# Load FindAllMarkers table ####
markers <- read_excel("SchiReyBan_final_allMarkers_cell_type.xlsx", sheet = 1)
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)

# Define cell type signatures - play with parameters here ####
#tmp <- markers[which(markers$p_val_adj < 1E-100 & markers$log_FC > 2),]

#use top 40 per cell type ####
tmp2 <- markers[order(-markers$log_FC),]
tmp2 <-group_by(tmp2,cluster)
tmp <-top_n(tmp2, 40, log_FC)

signatures <- split(tmp$gene, as.character(tmp$cluster))
signatures <- lapply(signatures, function(x) x[1:100])



# Run deconvolution ####
runDeconvolve <- function(fcs){
  ks_results <- do.call(rbind, lapply(names(signatures), function(sig){
    print(sig)
    sig <- signatures[[sig]]
    ok <- intersect(unlist(signatures), names(fcs))
    fcs <- fcs[ok]
    
    ok <- intersect(sig, names(fcs))
    if(length(ok) <=5) return(c(NA, NA))
    set <- fcs[match(ok, names(fcs))]
    rest <- fcs[-match(ok, names(fcs))]
    ks <- ks.test(set, rest)
    c(ks$p.value, mean(set) - mean(rest))
  }))
  rownames(ks_results) <- names(signatures)
  colnames(ks_results) <- c("pval", "fc_diff")
  ks_results
}
deconvolve_ipf1 <- runDeconvolve(fcs_ipf1)
deconvolve_ipf2 <- runDeconvolve(fcs_ipf2)
deconvolve_ipf3 <- runDeconvolve(fcs_ipf3)


# Generate heatmap for cell type deconvolution of IPFs ####
signedPval <- function(pval, fcs){
  pval[which(pval < 1e-50)] <- 1e-50
  pval <- (-log10(pval))
  pval[which(fcs < 0)] <- pval[which(fcs < 0)]*(-1)
  pval
}
tmp <- cbind(signedPval(deconvolve_ipf1[,"pval"], deconvolve_ipf1[,"fc_diff"]),
             signedPval(deconvolve_ipf2[,"pval"], deconvolve_ipf2[,"fc_diff"]),
             signedPval(deconvolve_ipf3[,"pval"], deconvolve_ipf3[,"fc_diff"]))
colnames(tmp) <- c("IPF1", "IPF2", "IPF3")
tmp
tmp[which(is.na(tmp))]<-0
pheatmap(tmp, cluster_cols = F, breaks = seq(-20, 20, length = 101), color= cm_colors)

### Figure 3f ####
# Generate heatmaps for genes driving the specific cell type signatures ####
genHeatmap_mean <- function(signature, orderpval = F){
  nom <- signature
  signature <- signatures[[nom]]
  ok <- intersect(signature, rownames(zscore))
  anno_col <- data.frame(treat)
  rownames(anno_col) <- colnames(zscore)
  
  if(orderpval){
    tmp <- markers[which(markers$cluster == nom),]
    tmp <- tmp$log_FC[match(ok, tmp$gene)]
    ok <- ok[order(tmp)]
  }
  
  asplit <- split(colnames(zscore), treat)
  means <- do.call(cbind, lapply(asplit, function(x) rowMeans(zscore[,x])))
  
  if(orderpval){
    pheatmap(means[ok,], cluster_cols = F, cluster_rows = F, color= cm_colors)  
  }
  if(!orderpval){
    pheatmap(means[ok,], cluster_cols = F, color= cm_colors)  
  }
  
}

## The heatmaps show z-scores for the individual marker genes of the indicated cell types across IPF stages and controls.

genHeatmap_mean("Ciliated cells", orderpval = T)
genHeatmap_mean("Myofibroblasts", orderpval = T)
genHeatmap_mean("Plasma cells", orderpval = T)



