## Figure 2: Multi-cohort single cell data reveals transcriptional changes in >40 cell types and altered cell type frequencies in diseases progression
## for panels a, b, and c check the scanpy code for the single cell objects

## panel 2d ##
## files needed:
## SchiReyBan_final_mapped_ids_metadata.txt
## 191210_SchiReyBan_counts.mtx
## 191204_SchiReyBan_merged_raw_genes.csv


## Get likelihood of detection from merged single cell object first##

# Load R libs ####
library(Matrix)
library(slam)
setwd("/home/agando/Documents/Christoph/R on redPC/LUKAS files/")

# Load objects #### For metacelltypes combine cell types to meta categories###
meta <- read.delim("SchiReyBan_final_mapped_ids_metadata.txt")
conversion <- rbind(
  c("Aberrant basaloid cells",	"Aberrant basaloid cells"),
  c("AdvF (PI16+)",	"AdvF (PI16+)"),
  c("AdvF (SFRP2+)",	"AdvF (SFRP2+)"),
  c("AM",	"AM"),
  c("AM activated",	"AM activated"),
  c("AM groundstate",	"AM"),
  c("Art_EC",	"Art_EC"),
  c("AT-1 cells",	"AT-1 cells"),
  c("AT-2 cells",	"AT-2 cells"),
  c("B-cells",	"B cells"),
  c("B cells",	"B cells"),
  c("Basal cells",	"Basal cells"),
  c("Bro_EC",	"Bro_EC"),
  c("Cap-A_EC",	"Cap-A_EC"),
  c("Cap-B_EC",	"Cap-B_EC"),
  c("CD4 M/E",	"CD4 M/E"),
  c("CD4 Na",	"CD4 Na"),
  c("CD8 M/E",	"CD8 M/E"),
  c("CD8 Na",	"CD8 Na"),
  c("Ciliated cells",	"Ciliated cells"),
  c("Club cells",	"Club cells"),
  c("DC EREG",	"DC EREG"),
  c("DC IGSF21",	"DC IGSF21"),
  c("Differentiating Ciliated cells",	"Differentiating Ciliated cells"),
  c("empty", ""),
  c("Goblet cells",	"Goblet cells"),
  c("Inflammatory fibroblast 1",	"Inflammatory fibroblast 1"),
  c("Inflammatory fibroblast 3",	"Inflammatory fibroblast 3"),
  c("Lipofibroblast",	"Lipofibroblast"),
  c("Lym_EC",	"Lym_EC"),
  c("Mast cells",	"Mast cells"),
  c("mDC1",	"mDC1"),
  c("mDC2",	"mDC2"),
  c("MegaK", ""),
  c("Mesothelial cells",	"Mesothelial cells"),
  c("MKI67+ cells", ""),
  c("Mono class",	"Mono class"),
  c("Mono NC",	"Mono NC"),
  c("Myofibroblasts",	"Myofibroblasts"),
  c("Neu",	"Neu"),
  c("NK",	"NK"),
  c("NKT",	"NK"),
  c("pDC",	"pDC"),
  c("Pericytes",	"Pericytes"),
  c("Pericytes activated",	"Pericytes activated"),
  c("Plasma cells",	"Plasma cells"),
  c("preAT-1 cells",	"preAT-1 cells"),
  c("SMCs",	"SMCs"),
  c("Transitional Lipo/Myo",	""),
  c("Vein_EC",	"Vein_EC"))

meta$metacelltype <- NA
lapply(1:nrow(conversion), function(x){
  meta$metacelltype[which(meta$cell_type == conversion[x, 1])] <<- conversion[x, 2]
})

cts <- readMM("/mnt/smb/Christoph/191210_SchiReyBan_counts.mtx")
cts2 <- as.simple_triplet_matrix(cts)
genes <- read.csv("/mnt/smb/Christoph/191204_SchiReyBan_merged_raw_genes.csv")
genes <- as.character(genes$colummn)

# Split by patient and cell type ####
groups <- paste(meta$metacelltype, meta$patient_id, sep = "|")
asplit <- split(1:nrow(meta), groups)

# Calculate likelihood of detection ####
expr <- do.call(cbind, lapply(asplit, function(x) colapply_simple_triplet_matrix(cts2[x, ], function(x) mean(x > 0))))
rownames(expr) <- genes

# Calculate number of cells ####
num_cells <- unlist(lapply(asplit, length))

# Calculate the number of total reads ####
total_reads <- unlist(lapply(asplit, function(x) sum(cts2[x, ])))

# Get phenotype info on a sample-level ####
samples <- unique(unlist(lapply(colnames(expr), function(x) strsplit(x, "|", fixed = T)[[1]][2])))
meta_tmp <- meta[, c("Age", "patient_id", "Sex", "Smoking.Status", "Smoking.History", "health_status", "data_set", "disease")]
meta_tmp <- meta_tmp[match(samples, meta_tmp$patient_id),]
sample_info <- meta_tmp

save(expr, num_cells, total_reads, sample_info,
     file = "LikeliDetec_metacelltypes_new_CM.RData")
#save(expr, num_cells, total_reads, sample_info,
#     file = "LikeliDetec_allcelltypes_new_CM.RData")



## Calculate DGE based on modalities in sample_info like health_status (healthy vs ILD), smoking_history (yes oder no) or Age#### 
## change in line 153, 156 and 162 accordingly e.g. smoking_history, smoking_history, smoking_historyYes

## files needed that were generated above:
## LikeliDetec_allcelltypes_new_CM.RData for all cell types
## LikeliDetec_metacelltypes.RData for metacell types

# Load data ####
library(readxl)
library(pheatmap)
library(BiocManager)
library(preprocessCore)
library(dplyr)
library(plotrix)

setwd("/home/agando/Documents/Christoph/R on redPC/LUKAS files/")
#load("LikeliDetec_allcelltypes_new_CM.RData")
load("LikeliDetec_metacelltypes.RData")

# Define objects ####
cellgroups <- unlist(lapply(colnames(expr), function(x) strsplit(x, "|", fixed = T)[[1]][1]))
metacelltypes <- setdiff(unique(cellgroups), "")
samples <- unlist(lapply(colnames(expr), function(x) strsplit(x, "|", fixed = T)[[1]][2]))

# Run regression ####
runRegression <- function(celltype){
  print(celltype)
  
  expr_subset <- expr[,which(cellgroups == celltype)]
  expr_subset <- expr_subset[which(apply(expr_subset, 1, function(x) sum(x > 0.1)) > 10), ]
  expr_subset <- sqrt(expr_subset)
  
  if(ncol(expr_subset) < 10) return(NULL)
  
  total_reads_subset <- total_reads[which(cellgroups == celltype)]
  num_cells_subset <- num_cells[which(cellgroups == celltype)]
  
  samples_subset <- samples[which(cellgroups == celltype)]
  meta_tmp_subset <- sample_info[match(samples_subset, sample_info$patient_id),]
  
  res <- t(apply(expr_subset, 1, function(x){
    aframe <- data.frame(expr = x, meta_tmp_subset, total_reads_subset, num_cells_subset)
    
    if(length(unique(aframe$data_set)) > 1){
      afit <- try(lm(expr ~ total_reads_subset + num_cells_subset + data_set + health_status, data = aframe))
    }
    if(length(unique(aframe$data_set)) == 1){
      afit <- try(lm(expr ~ total_reads_subset + num_cells_subset + health_status, data = aframe))
    }
    
    if(class(afit) == "try-error") return(NULL)
    
    coefs <- coefficients(summary(afit))
    coefs["health_statusILD", c(3, 4)]
  }))
  res[sort.list(res[,2]),]
}
res <- lapply(metacelltypes, runRegression)

names(res) <- metacelltypes

# Add average detection likelihoods ####
res2 <- lapply(metacelltypes, function(x){
  avg_detection <- rowMeans(expr[,which(cellgroups == x)])
  tmp <- data.frame(res[[x]])
  tmp$avg_detection <- avg_detection[rownames(tmp)]
  colnames(tmp) <- c("t.value", "p.value", "avg.detection")
  tmp <- cbind("genes" = rownames(tmp), tmp)
  rownames(tmp) <- NULL
  tmp
})
names(res2) <- names(res)


# Save files ####
path <- "/home/agando/Documents/Christoph/R on redPC/LUKAS files/dge/metacelltypes_health-status/"
lapply(names(res2), function(x) write.csv(res2[[x]], file = paste0(path, x, ".csv")))

## Combine all files into one file
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files/dge/metacelltypes_health-status/")
files <- list.files()
celltypes <- gsub(".csv", "", files)

tab <- data.frame(gene = NA)
for(celltype in celltypes){
  print(paste("Adding", celltype))
  tmp <- read.delim(paste0(celltype, ".csv"), sep = ",")
  colnames(tmp) <- c("gene", paste(celltype, sep = "_", colnames(tmp)[-1]))
  tab <- merge(tab, tmp, by = "gene", all = TRUE)
}

## Sort Data Frame
order <- c(1, grep("t.value", colnames(tab)), grep("p.value", colnames(tab)), grep("detect", colnames(tab)))
tab <- tab[-which(is.na(tab$gene)), order]

## Write output
write.table(tab, file = "metacelltypes_health-status_merged.txt", quote = F, sep = "\t", row.names = F)


## Final figure panels 2d and 2e were generated by clustering the t.values of the corresponding metacelltypes in the Perseus software

## Figure 2f, 2g, 2h
## Generate DGE vox plots showing %detection
## Files needed:
## LikeliDetec_metacelltypes.RData

# Load DGE results ####
setwd("N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files")
path <- "N:/AGSchiller/Christoph/ANALYSIS/191118_ASK-MLT_multiomics/deconvolution/LUKAS files/dge/metacelltypes_health-status/"

files <- list.files(path, full.names = T)
res2 <- lapply(files, function(x) read.csv(x, row.names = 1))
nom <- gsub(path, "", files)
nom <- gsub(".csv", "", fixed = T, nom)
names(res2) <- nom
load("LikeliDetec_metacelltypes.RData")
#load("LikeliDetec_allcelltypes_new_CM.RData")

# Define objects ####
cellgroups <- unlist(lapply(colnames(expr), function(x) strsplit(x, "|", fixed = T)[[1]][1]))
metacelltypes <- setdiff(unique(cellgroups), "")
samples <- unlist(lapply(colnames(expr), function(x) strsplit(x, "|", fixed = T)[[1]][2]))

# Define plotting function ####
genPlot <- function(gene = "CRTAC1", celltype = "EPCAM+ Epithelium", treat = "health_status"){
  gene_expr <- expr[gene, ]
  ok <- which(cellgroups == celltype)
  gene_expr <- gene_expr[ok]
  aframe <- data.frame(gene_expr, sample_info[match(samples[ok], sample_info$patient_id),])
  par(mar = c(5, 5,5,5))
  boxplot(split(aframe$gene_expr, paste(aframe$data_set, aframe[,treat])),
          ylab = paste(gene, "detection [%]"), main = celltype, las = 2, col = c("#440154FF", "#FDE725FF"))
}

# Visualize some examples ####
par(mfrow = c(2,2))

genPlot("KRT17", "alveolar epithelium", "health_status")
genPlot("CDKN2A", "alveolar epithelium", "health_status")

genPlot("DIO2", "fibroblast", "health_status")
genPlot("CXCL14", "fibroblast", "health_status")

genPlot("SPP1", "AM", "health_status")
genPlot("CCL7", "AM", "health_status")

genPlot("CRTAC1", "alveolar epithelium", "health_status")
