## Figure 5 Protein signatures in BALF predict lung function decline and the corresponding cellular changes
## Panel 5f: The heatmap illustrates gene expression changes associated with lung fibrosis across indicated cell types (columns) for selected BALF protein biomarkers (rows).
## The dotplot visualizes the frequency changes of the indicated cell types inferred from the deconvolution of bulk mRNA data of ILD samples compared to controls (GSE47460). Samples used in this study from the LTRC n= 254 ILD patients and n= 108 controls.

# files needed:
# SchiReyBan_final_allMarkers_meta_cell_type.csv
# Dataset EV4b.xlsx


# Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(xlsx)
library(ggrepel)
library(reshape)
library(gridExtra)

# load series and platform data from GEO
gset <- getGEO("GSE47460", GSEMatrix =TRUE, AnnotGPL=TRUE)[[1]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "undefined"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data
outcome <- gset@phenoData@data$`disease state:ch1`
outcome <- gsub("Chronic Obstructive Lung Disease", "COPD", fixed = T, outcome)
outcome <- gsub("Control", "Ctrl", fixed = T, outcome)
outcome <- gsub("Interstitial lung disease", "ILD", fixed = T, outcome)
gset$description <- outcome
design <- model.matrix(~ description + 0, gset)
colnames(design) <- gsub("description", "", colnames(design))

# Run differential expression (COPD vs ctrl & ILD vs ctrl)
fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(ILD-Ctrl, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
res_ild <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(gset))

# Load cell type markers ####
markers <- read.csv("SchiReyBan_final_allMarkers_meta_cell_type.csv")
markers$p_val_adj <- as.numeric(markers$p_val_adj)
markers$p_val <- as.numeric(markers$p_val)

# Define cell type signatures - play with parameters here ####
tmp <- markers[which(markers$p_val_adj < 0.1 & markers$log_FC > 1),]
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
fcs_ild <- res_ild$logFC
names(fcs_ild) <- res_ild$GENE_SYMBOL
deconvolve_ild <- runDeconvolve(fcs_ild)

# Plot volcano ####
genVolcano <- function(res, main){
  res <- data.frame(res)
  res$pval[which(res$pval < 1e-50)] <- 1e-50
  res$pval <- (-log10(res$pval))
  ggplot(res, aes(fc_diff, pval)) + geom_point() +
    geom_text_repel(data = res, aes(label = rownames(res))) +
    geom_vline(xintercept = 0) +
    ggtitle(main)
}
genVolcano(deconvolve_ild, "Microarray ILD vs ctrl")

# Load t-value tables ####
tmp <- read.xslx("Dataset EV4b.xlsx", row.names = 1)
tvals <- tmp[, grep('t.value', fixed = T, colnames(tmp))]
pvals <- tmp[, grep('p.value', fixed = T, colnames(tmp))]

# Load BALF regression results ####
balf_res <- read_excel('Protein_lung_function_correlation_regression.xlsx', sheet = 1)
res <- balf_res

sig <- which(res$adjusted_pval < 0.25)

aframe <- data.frame(coef = res$estimate, pval = -log10(res$p_value), gene = res$Gene.names)
ggplot(aframe, aes(x = coef, y = pval, color = pval)) + geom_point(alpha = 0.5) +
  xlab('Coefficient') + ylab('Significance [-log10 p-value]') +
  scale_colour_gradient2(low = "black", mid = 'grey', high = "red") + 
  geom_text_repel(data = aframe[sig,], aes(label = gene)) +
  theme_bw()

# Define set of genes to plot
#genes <- na.omit(unique(res$Gene.names[sig]))
#genes <- unlist(lapply(genes, function(x) strsplit(x, ';', fixed = T)[[1]][1]))
genes <- read_excel('heatmap+rf.xlsx')
genes <- genes$genes
genes <- intersect(genes, rownames(tvals))

# Create heatmap of t-values ####
load('BlueYellowColormaps_V1.RData')
tmp <- data.matrix(tvals[genes,])
tmp[which(is.na(pvals[genes,]))] <- 0
#tmp[which(pvals[genes,] > 0.1)] <- 0

colnames(tmp) <- gsub('.t.value', '', fixed = T, colnames(tmp))
nom <- colnames(data.frame(t(deconvolve_ild)))
nom <- intersect(colnames(tmp), nom)
tmp <- tmp[, nom]

p <- pheatmap(tmp,
              breaks = seq(-4, 4, length = length(blue2yellow) + 1), col = yellow2blue)
ord <- p$tree_col$order

# Define function for dotplot ####
DotPlot <- function(
  tval,
  pval,
  cex.use = 0.5,
  dot_max = 5,
  cols.use = NULL,
  thresh.col = 2.5,
  dot.min = 0.25,
  title)
{
  avg.exp <- data.matrix(tval)
  avg.exp[which(is.na(avg.exp))] <- 0
  avg.alpha <- data.matrix(pval)
  avg.alpha[which(is.na(avg.alpha))] <- 1
  avg.alpha <- -log10(avg.alpha)
  ii <- cut(avg.exp, breaks = seq(-4, 4, len = 100), include.lowest = TRUE)
  cols.use <- colorRampPalette(c("darkblue", "darkred"))(99)[ii]
  cols.use <- colorRampPalette(c('red', 'grey', 'blue'))(99)[ii]
  n.col <- length(cols.use)
  data.y <- rep(x = 1:ncol(x = avg.exp), nrow(x = avg.exp))
  data.x <- unlist(x = lapply(X = 1:nrow(x = avg.exp), FUN = rep, ncol(x = avg.exp)))
  data.avg <- unlist(x = lapply(
    X = 1:length(x = data.y),
    FUN = function(x) {
      return(avg.exp[data.x[x], data.y[x]])
    }
  ))
  
  tmp <- data.avg
  tmp[which(tmp > 2)] <- 2
  tmp[which(tmp < (-2))] <- (-2)
  exp.col <- color.scale(tmp, extremes = c('blue', 'red'))
  
  data.cex <- unlist(x = lapply(
    X = 1:length(x = data.y),
    FUN = function(x) {
      return(avg.alpha[data.x[x], data.y[x]])
    }
  )) * cex.use + dot.min
  
  data.cex[which(data.cex > dot_max)] <- dot_max
  
  plot(
    x = data.y,
    y = data.x,
    cex = data.cex,
    pch = 16,
    col = exp.col,
    xaxt = "n",
    xlab = "",
    ylab = "",
    yaxt = "n",
    main = title
  )
  axis(side = 2, at = 1:nrow(x = avg.alpha), rownames(avg.alpha), las = 1)
  axis(side = 1, at = 1:ncol(x = avg.alpha), colnames(avg.alpha), las = 2)
}

# Dotplot of cell type frequency changes ####
tmp <- deconvolve_ild[ord,]
tmp_pvals <- t(data.matrix(tmp[,1]))
tmp_fcs <- t(data.matrix(tmp[,2]))
par(mar = c(10, 2, 5, 2))
DotPlot(tmp_fcs, tmp_pvals,
        dot.min = 1, cex.use = 0.2, dot_max = 4,
        title = 'Frequency changes')