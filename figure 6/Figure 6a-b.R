## Figure 6 Cell type specific transcriptional ILD signatures translate into protein signatures associated with lung function in the BALF
##
## Panel 6a-b: Scatter plots stratify genes based on the protein lung function associations (y-axis) and cell type specific ILD associations (x-axis). 
## The size of the dots represents the detection level of each gene in the corresponding cell type. Colors highlight genes with marginal associations at the protein and RNA levels.

## files needed:
## Protein_lung_function_correlation_regression.xlsx
## Dataset EV4b.xlsx
## BALF_imputed.RData ## from Figure 4 ##

# Load R libs ####
library(ggplot2)
library(xlsx)
library(missMDA)
library(cowplot)
theme_set(theme_cowplot())
library(ggrepel)

# Load BALF data ####
balf_res <- read.xlsx(file = 'Protein_lung_function_correlation_regression.xlsx', sheetIndex = 1)
balf_res <- balf_res[sort.list(balf_res$adjusted_pval),]

tmp <- do.call(rbind, lapply(1:nrow(balf_res), function(x){
  genes <- strsplit(balf_res$Gene.names[x], ';', fixed = T)[[1]]
  if(length(genes) > 1){
    return(data.frame(gene = genes, pos = x))
  }
  if(length(genes) == 1){
    return(data.frame(gene = genes, pos = x))
  }
}))
balf <- balf_res[tmp$pos,]
balf$gene <- tmp$gene

# Load scRNAseq data ####
ild <- read.xslx('Dataset EV4b.xlsx', row.names = 1)
pvals <- ild[,grep('p.value', colnames(ild), fixed = T)]
tvals <- ild[,grep('t.value', colnames(ild), fixed = T)]
detection <- ild[,grep('detection', colnames(ild), fixed = T)]
colnames(detection) <- colnames(tvals) <- colnames(pvals) <- gsub('.p.value', '', fixed = T, colnames(pvals))

# Find common genes ###
ok <- intersect(rownames(tvals), balf$gene)
correls <- as.numeric(cor(tvals[ok,], as.numeric(balf$t_value[match(ok, balf$gene)]),
                          use = 'pairwise.complete'))
names(correls) <- colnames(pvals)

gene_counts <- unlist(lapply(colnames(tvals), function(x) length(which(!is.na(tvals[ok,x])))))
names(gene_counts) <- colnames(tvals)
good <- names(which(gene_counts > 500))

# Plot correlations as barplot ####
aframe <- data.frame(cor = correls, name = colnames(pvals))
aframe <- aframe[which(aframe$name %in% good),]
aframe$name <- factor(as.character(aframe$name), levels = rev(aframe$name[order(aframe$cor)]))
aframe$sig <- aframe$cor > 0
ggplot(data=aframe, aes(x=name, y=cor, fill = sig)) +
  geom_bar(stat="identity") + coord_flip() + xlab('Metacelltype') +
  scale_fill_manual(values=c("red", "blue")) +
  ylab(paste('Correlation coefficient',
             '[BALF lung function & scRNAseq ILD]', sep = '\n'))

# Plot single comparisons ####
gen_scatter_plot <- function(celltype, pval = 1e-3){
  aframe <- data.frame(scrnaseq_tval = tvals[ok, celltype],
                       scrnaseq_pval = pvals[ok, celltype],
                       scrnaseq_detection = detection[ok, celltype],
                       balf_tval = as.numeric(balf$t_value[match(ok, balf$gene)]),
                       balf_pval = as.numeric(balf$p_value[match(ok, balf$gene)]),
                       gene = ok)
  aframe$sig <- 'neither'
  aframe$sig[which(aframe$scrnaseq_pval < pval)] <- 'scRNAseq'
  aframe$sig[which(aframe$balf_pval < pval)] <- 'BALF'
  aframe$sig[which(aframe$balf_pval < pval & aframe$scrnaseq_pval < pval)] <- 'both'
  
  group.colors <- c(neither = "grey", scRNAseq = "red", BALF ="blue", both = "purple")
  
  intercept <- coefficients(aov(aframe$balf_tval ~ aframe$scrnaseq_tval))[1]
  slope <- coefficients(aov(aframe$balf_tval ~ aframe$scrnaseq_tval))[2]
  
  ggplot(aframe, aes(x = scrnaseq_tval, y = balf_tval, color = sig, size = scrnaseq_detection)) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values=group.colors) + 
    geom_vline(xintercept = 0) + geom_hline(yintercept = 0) +
    ggtitle(celltype) + labs(color = paste('Significant [<', pval, ']')) +
    xlab('t-value [scRNAseq - ILD status]') + ylab('t-value [BALF - meta lung function]') +
    #geom_abline(intercept = intercept, slope = slope) +
    geom_text_repel(data = aframe[which(aframe$sig != 'neither'),],
                    aes(label = gene, color = sig))
}

gen_scatter_plot('alveolar.epithelium')
gen_scatter_plot('Basal.cells')
gen_scatter_plot('Club.cells')

# Highlight some example genes ####
load('BALF_imputed.RData') ## from figure 4 ##

lung_function <- num[, colnames(num)[c(3, 5:8, 13)]]
lung_function_imputed <- imputePCA(lung_function)$completeObs
pca <- prcomp(lung_function_imputed, scale. = T)
meta_lung_function <- pca$x[,1]

examples <- c('TNC', 'SFN', 'MMP7', 'SUSD2')

gen_scatter <- function(gene){
  ok <- grep(gene, prot_info$Gene.names)
  plot(imputed[ok,], meta_lung_function,
       ylab = 'Meta lung function', xlab = 'Expression', pch = 16, main = gene)
  abline(coefficients(aov(meta_lung_function ~ imputed[ok,])))
  legend('topleft', paste(paste('Coefficient:', signif(balf$estimate[which(balf$gene == gene)], 2)),
                          paste('P-value:', signif(balf$adjusted_pval[which(balf$gene == gene)], 2)), sep ='\n'), bty = 'n')
}

par(mfrow = c(2, 2))
lapply(examples, gen_scatter)