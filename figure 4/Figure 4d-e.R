## Figure 4: Human lung bronchoalveolar lavage fluid proteome changes correlate with clinical parameters.

## panel 4d) The heatmap illustrates the computationally derived meta lung function variable combining multiple lung function parameters.
## (e) The volcano plot shows the multivariate regression coefficients (x-axis) and the -log10 p-value (y-axis) for BALF protein abundance with the meta lung function.

# This code will perform calculation and plotting of the meta lung function parameter
# files needed:
# BALF_proteome_clinical_numerical+categorical_CM.txt
# BALF_proteome_CM.txt



# Load R libs ####
library(Matrix)
library(randomForest)
library(missForest)
library(missMDA)
library(pheatmap)
library(ggrepel)
library(xlsx)

# Load data ####
setwd('E:/Dropbox/PhD/HUMAN paper/proteomics')
dat <- read.table('BALF_proteome_clinical_numerical+categorical_CM.txt', header = T, sep = '\t', row.names = 1)
cat <- dat[1:4]
num <- dat[5:37]
prot <- read.table('BALF_proteome_CM.txt', sep = '\t', header = T)

# Restrict to samples with protein and outcome data ####
samples <- intersect(colnames(prot), rownames(cat))
cat <- cat[samples,]
num <- num[samples,]
prot_expr <- data.matrix(prot[, samples])
prot_info <- prot[,1:4]
rownames(prot_info) <- rownames(prot_expr) <- as.character(1:nrow(prot_expr))

# Impute NAs ####
imputed <- missForest(t(prot_expr))
imputed <- imputed$ximp
imputed <- t(imputed)

save(imputed, prot_info, num, cat, file = 'BALF_imputed_CM.RData')

tmp <- cbind(t(prot_expr), lung_function)
tmp <- missForest(tmp)$ximp

imputed <- t(tmp[,1:nrow(prot_expr)])
lung_function <- tmp[,(nrow(prot_expr)+1) : ncol(tmp)]

# Calculate meta lung function variable ####
lung_function <- num[, colnames(num)[c( 5:8, 13)]]
mean(is.na(lung_function))
lung_function_imputed <- imputePCA(lung_function)$completeObs
pca <- prcomp(lung_function_imputed, scale. = T)

anno <- data.frame(meta_lung_function = pca$x[,1], gender = cat$sex, age = num$age..years.)
anno$gender[which(anno$gender == '')] <- NA
rownames(anno) <- rownames(lung_function_imputed)
pheatmap(t(lung_function_imputed[order(pca$x[,1]),]),
         cluster_cols = F, annotation_col = anno, scale = 'row',
         show_colnames = F, breaks = seq(-2, 2, length = 101))

# Associate protein expression with lung function accounting for gender and age ####
lung_function <- pca$x[,1]
res <- t(apply(imputed, 1, function(x){
  afit <- lm(x ~ anno$gender + anno$age + anno$meta_lung_function) 
  tmp <- coefficients(summary(afit))
  tmp[nrow(tmp),]
}))
colnames(res) <- c('estimate', 'std_error', 't_value', 'p_value')
res <- data.frame(prot_info[,c('Protein.IDs', 'Gene.names')], res)
res$adjusted_pval <- p.adjust(res$p_value, method = 'BH')

sig <- which(res$adjusted_pval < 0.25)

aframe <- data.frame(coef = res$estimate, pval = -log10(res$p_value), gene = prot_info$Gene.names)
ggplot(aframe, aes(x = coef, y = pval, color = pval)) + geom_point(alpha = 0.5) +
  xlab('Coefficient') + ylab('Significance [-log10 p-value]') +
  scale_colour_gradient2(low = "black", mid = 'grey', high = "red") + 
  geom_text_repel(data = aframe[sig,], aes(label = gene))

resids <- t(apply(imputed[sig,], 1, function(x) residuals(lm(x ~ anno$gender + anno$age))))

load('BlueYellowColormaps_V1.RData')

rownames(anno) <- colnames(imputed)
pheatmap(imputed[sig,order(anno$meta_lung_function)],
         annotation_col = anno, scale = 'row', labels_row = res$Gene.names[sig],
         cluster_cols = F, breaks = seq(-3, 3, length = 256), show_colnames = F,
         color = yellow2blue)

# Plot CRTAC1 association with meta lung function ####
plot(imputed[328,], anno$meta_lung_function,
     ylab = 'Meta lung function', xlab = 'CRTAC1 expression', pch = 19)
abline(coefficients(aov(anno$meta_lung_function ~ imputed[328,])))
legend('topleft', paste(paste('Coefficient:', signif(res$estimate[328], 2)),
                        paste('P-value:', signif(res$p_value[328], 2)), sep ='\n'), bty = 'n')

# Plot CRTAC1 association with meta lung function ####
plot(imputed[29,], anno$meta_lung_function,
     ylab = 'Meta lung function', xlab = 'TNC expression', pch = 19)
abline(coefficients(aov(anno$meta_lung_function ~ imputed[29,])))
legend('topleft', paste(paste('Coefficient:', signif(res$estimate[29], 2)),
                        paste('P-value:', signif(res$p_value[29], 2)), sep ='\n'), bty = 'n')


# Plot CRTAC1 association with meta lung function ####
plot(imputed[234,], anno$meta_lung_function,
     ylab = 'Meta lung function', xlab = 'CFHR1 expression', pch = 19)
abline(coefficients(aov(anno$meta_lung_function ~ imputed[234,])))
legend('topleft', paste(paste('Coefficient:', signif(res$estimate[234], 2)),
                        paste('P-value:', signif(res$p_value[234], 2)), sep ='\n'), bty = 'n')

write.xlsx(res, file = 'Protein_lung_function_correlation_regression_CM_CFHR1.xlsx')
