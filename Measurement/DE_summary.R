if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")

library(DESeq2)
library(tidyverse)
library(Matrix)
# library(affy)
library(limma)
# library(DEqMS)

## SRTP2021 Differential Expression Code Summary

# A
setwd("D:\\SRTP2021\\1")
# read & preprocess
matrix <- read.csv("GSE145642_Amor_RNASeq.cnt.csv")
id_name <- matrix[,c(1,10,17)]
rownames(matrix) <- matrix[,1]
matrix <- matrix[,2:9]

# Hep DESeq2
matrix_Hep <- matrix[,1:4]
condition_Hep <- factor(c("a", "y", "a", "y"),
                   levels = c("y", "a"))
colData_Hep <- data.frame(row.names = colnames(matrix_Hep), condition_Hep)

dds <- DESeqDataSetFromMatrix(matrix_Hep, colData_Hep, design = ~condition_Hep)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\1\\DEGseq_result_Hep.csv')

# KP DESeq2
matrix_KP <- matrix[,5:8]
condition_KP <- factor(c("y", "y", "a", "a"),
                        levels = c("y", "a"))
colData_KP <- data.frame(row.names = colnames(matrix_KP), condition_KP)

dds <- DESeqDataSetFromMatrix(matrix_KP, colData_KP, design = ~condition_KP)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\1\\DEGseq_result_KP.csv')

# add id & go term
diff_gene_Group2 <- read.csv("D:\\SRTP2021\\1\\DEGseq_result_Hep.csv")
diff_gene_Group2 <- left_join(diff_gene_Group2, id_name, by = c("X" = "ID"))
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\1\\DEGseq_result_Hep_with_name.csv')

diff_gene_Group2 <- read.csv("D:\\SRTP2021\\1\\DEGseq_result_KP.csv")
diff_gene_Group2 <- left_join(diff_gene_Group2, id_name, by = c("X" = "ID"))
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\1\\DEGseq_result_KP_with_name.csv')

# A GSE152329
matrix <- read.table("D:\\SRTP2021\\10\\GSE152329_ccl4_liver_fibrosis.kallisto.counts.txt", header = T)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
matrix <- matrix[,order(colnames(matrix))]
matrix <- round(matrix)
condition <- factor(c(rep("a", 99), rep("y", 99)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.01 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\10\\DEGseq_result.csv')



# B

# df_selected_final liver
df_selected_final <- read.csv("D:/SRTP2021/12/liver/df_liver_selected_final_11.30.csv")
matrix <- df_selected_final %>%
  subset(select=c('Ly6e', 'Cd52', 'Ifitm3', 'Itm2b', 'Ccl5', 'Rpl3','Ier3','Ifitm2',
                  'Fxyd5','Cyba', 'target', 'B2m'))

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, -11]
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
# write to csv
write.csv(tT_up, file='D:\\SRTP2021\\12\\limma_result_liver_more.csv')


# heart
df_selected_final <- read.csv("D:/SRTP2021/12/heart/df_selected_final_heart.csv")
markers <- str_to_title(c('Prg4', 'Mdk', 'Mgll', 'Ifitm1', 'Plvap'))
matrix <- df_selected_final %>% select(c(markers,'target'))
  # Columns `Vsir`, `Ap2m1`, `Anxa3`, `Selp`, `Atp2b1`, etc. don't exist

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, colnames(matrix)!='target']
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
# write to csv
write.csv(tT_up, file='D:\\SRTP2021\\12\\heart\\limma_result_heart.csv')


# heart_full
df_selected_final <- read.csv("D:/SRTP2021/12/heart/df_heart.csv")
matrix <- df_selected_final # %>% select(c(markers,'target'))

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, colnames(matrix)!='target']
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
tT_up_sig <- tT_up[tT_up$adj.P.Val<0.05,]
# write to csv
write.csv(tT_up_sig, file='D:\\SRTP2021\\12\\heart\\limma_result_heart_sig.csv')


# marrow
df_selected_final <- read.csv("D:/SRTP2021/12/model_marrow/df_selected_final_11.30.csv")
matrix <- df_selected_final %>%
  subset(select=c('Ly6e', 'Cd52', 'Ifitm3', 'Itm2b', 'Ccl5', 'Rpl3','Ier3','Ifitm2',
                  'Fxyd5','Cyba', 'target', 'B2m'))

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, -11]
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
# write to csv
write.csv(tT, file='D:\\SRTP2021\\12\\limma_result_more.csv')

# lung 
df_selected_final <- read.csv("D:/SRTP2021/12/lung/df_lung_selected_final_11.30.csv")
matrix <- df_selected_final %>%
  subset(select=c('Ly6e', 'Cd52', 'Ifitm3', 'Itm2b', 'Ccl5', 'Rpl3','Ier3','Ifitm2',
                  'Fxyd5','Cyba', 'target', 'B2m'))

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, -11]
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
# write to csv
write.csv(tT_up, file='D:\\SRTP2021\\12\\limma_result_lung_more.csv')

# spleen
matrix <- read.csv('D:/SRTP2021/12/spleen/df_spleen_short.csv')
matrix <- matrix[,-1]

# create 'design'
ages <- matrix$target
aged <- ages > 3
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[, -1]
expr <- t(expr)
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]
# write to csv
write.csv(tT, file='D:\\SRTP2021\\12\\limma_result_spleen_more.csv')



# C
matrix <- read.table("D:\\SRTP2021\\7\\GSE100906_scRNA_rawCounts.txt", header = T)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
condition <- factor(c(rep("y", 80), rep("a", 69)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.01 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\7\\DEGseq_result.csv')



# E
matrix302 <- read.table("D:\\SRTP2021\\14_human_proteome\\HTMGNBGX3_Donor_302_GENE_LEVEL_COUNTS_withGeneNames.txt", header = T)
matrix302_f <- matrix302[rowSums(matrix302>10)>5,colSums(matrix302>1000)>10]

matrix304 <- read.table("D:\\SRTP2021\\14_human_proteome\\HTMGNBGX3_Donor_304_GENE_LEVEL_COUNTS_withGeneNames.txt", header = T)
matrix304_f <- matrix304[rowSums(matrix304>10)>5,colSums(matrix304>1000)>10]

matrix305 <- read.table("D:\\SRTP2021\\14_human_proteome\\HTMGNBGX3_Donor_305_GENE_LEVEL_COUNTS_withGeneNames.txt", header = T)
matrix305_f <- matrix305[rowSums(matrix305>10)>5,colSums(matrix305>1000)>10]

matrix309 <- read.table("D:\\SRTP2021\\14_human_proteome\\HTMGNBGX3_Donor_309_GENE_LEVEL_COUNTS_withGeneNames.txt", header = T)
matrix309_f <- matrix309[rowSums(matrix309>10)>5,colSums(matrix309>1000)>10]


  # merge
aged <- merge(matrix302_f, matrix305_f, by=0)
colnames(aged) <- c('name', seq(1,240))

young <- merge(matrix304_f, matrix309_f, by=0)
colnames(young) <- c('name', seq(241, 542))

matrix <- merge(aged, young, by='name')
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]

round_matrix <- round(matrix)

condition <- factor(c(rep("a", 240), rep("y", 302)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(round_matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\14_human_proteome\\DEGseq_result.csv')
write.csv(matrix, file='D:\\SRTP2021\\14_human_proteome\\merged_data.csv')

  # limma for normalized counts
matrix <- read.csv("D:/SRTP2021/14_human_proteome/merged_data.csv")
# create 'design'
aged <- c(rep(T, 240), rep(F, 302))
design <- model.matrix(~0+aged)
# create 'expr'
expr <- matrix[,-1]
rownames(expr) <- expr[,1]
expr <- expr[,-1]
# fit model
fit <- lmFit(expr, design)
cts <- paste("agedTRUE","agedFALSE", sep = "-")
cont.matrix <- makeContrasts(contrasts = cts, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.05)
volcanoplot(fit2)
tT <- topTable(fit2, number = Inf, sort.by = 'logFC')
tT$name <- rownames(tT)

# only keep upregulated genes
tT_up <- tT[tT$logFC>0,]

# write to csv
write.csv(tT_up, file='D:\\SRTP2021\\14_human_proteome\\limma_result.csv')




# G
construct_matrix_rat <- function(n,name){
  mm <- readMM(paste('D:/SRTP2021/16_rat/GSM433182',n,'_BM-', name,'_matrix.mtx',sep = ''))
  matrix <- as.matrix(mm)
  barcode <- read.table(paste('D:/SRTP2021/16_rat/GSM433182',n,'_BM-', name,'_barcodes.tsv',sep = ''))
  genes <- read.table(paste('D:/SRTP2021/16_rat/GSM433182',n,'_BM-', name,'_genes.tsv',sep = ''))
  df <- as.data.frame(matrix)
  row.names(df) <- genes$V1
  colnames(df) <- barcode$V1
  return(df)
}

M_O <- construct_matrix_rat('3','M-O')
M_O_f <- M_O[rowSums(M_O>2)>5,colSums(M_O>2)>1000]

M_Y <- construct_matrix_rat('2','M-Y')
M_Y_f <- M_Y[rowSums(M_Y>2)>5,colSums(M_Y>2)>1000]

F_O <- construct_matrix_rat('6','F-O')
F_O_f <- F_O[rowSums(F_O>2)>5,colSums(F_O>2)>1000]

F_Y <- construct_matrix_rat('5','F-Y')
F_Y_f <- F_Y[rowSums(F_Y>2)>5,colSums(F_Y>2)>1000]

aged <- merge(M_O_f, F_O_f, by=0)
colnames(aged) <- c('name', seq(1,727))

young <- merge(M_Y_f, F_Y_f, by=0)
colnames(young) <- c('name', seq(728,1416))

matrix <- merge(aged, young, by='name')
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]

  # export for machine learning
matrix_con <- read.csv('D:\\SRTP2021\\DESeq2_results\\G_converted.csv')
matrix_con <- matrix_con[,-c(1,2,1420)]
matrix_con_noNA <- na.omit(matrix_con)
matrix_con_noNA_noDup <- matrix_con_noNA[!duplicated(matrix_con_noNA$name),]
rownames(matrix_con_noNA_noDup) <- matrix_con_noNA_noDup$name
matrix_con_noNA_noDup <- matrix_con_noNA_noDup[,-1417]
matrix_t <- as.data.frame(t(matrix_con_noNA_noDup))
matrix_t$target <- c(rep(27, 727), rep(5, 689))
write.csv(matrix_t, 'D:\\SRTP2021\\16_rat\\df_G.csv')

condition <- factor(c(rep("a", 727), rep("y", 689)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0.5)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\16_rat\\DEGseq_result.csv')
write.csv(matrix, file='D:\\SRTP2021\\16_rat\\matrixG.csv')


# F
construct_matrix_human <- function(gse, name){
  mm <- readMM(paste('D:/SRTP2021/15_human_rnaseq/GSM33961',gse,'_matrix_', name,'.mtx', sep = ''))
  matrix <- as.matrix(mm)
  # barcode <- read.table(paste('D:/SRTP2021/15_human_rnaseq/GSM33961',gse,'_barcodes_', name,'.tsv', sep = ''))
  genes <- read.table(paste('D:/SRTP2021/15_human_rnaseq/GSM33961',gse,'_genes_', name,'.tsv', sep = ''))
  df <- as.data.frame(matrix)
  row.names(df) <- genes$V1
  # colnames(df) <- barcode$V1
  return(df)
}

A_59 <- construct_matrix_human('61', 'A')
B_47 <- construct_matrix_human('62', 'B')
C1_60 <- construct_matrix_human('63', 'C1')
Ck_59 <- construct_matrix_human('64', 'Ck')
C2_60 <- construct_matrix_human('65', 'C2')
E_30 <- construct_matrix_human('66', 'E')
F_41 <- construct_matrix_human('67', 'F')
G_58 <- construct_matrix_human('68', 'G')
H_50 <- construct_matrix_human('69', 'H')
J_43 <- construct_matrix_human('70', 'J')
K_84 <- construct_matrix_human('71', 'K')
L_57 <- construct_matrix_human('72', 'L')
M_60 <- construct_matrix_human('73', 'M')
N_67 <- construct_matrix_human('74', 'N')
O_50 <- construct_matrix_human('75', 'O')
P_58 <- construct_matrix_human('76', 'P')
Q_66 <- construct_matrix_human('77', 'Q')
R_31 <- construct_matrix_human('78', 'R')
S1_56 <- construct_matrix_human('79', 'S1')
Sk1_55 <- construct_matrix_human('80', 'Sk1')
S2_56 <- construct_matrix_human('81', 'S2')
Sk2_55 <- construct_matrix_human('82', 'Sk2')
T_24 <- construct_matrix_human('83', 'T')
U_46 <- construct_matrix_human('84', 'U')
W_28 <- construct_matrix_human('85', 'W')

  # young: TW, aged: KN
colnames(T_24) <- seq(1,234)
colnames(W_28) <- seq(235,298)
young <- merge(T_24, W_28, by=0)
rownames(young) <- young[,1]
young <- young[,-1]

colnames(K_84) <- seq(299,365)
colnames(N_67) <- seq(366,459)
aged <- merge(K_84, N_67, by=0)
rownames(aged) <- aged[,1]
aged <- aged[,-1]

matrix <- merge(young, aged, by=0)
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]

condition <- factor(c(rep("y", 234), rep("a", 67)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

write.csv(matrix, file='D:\\SRTP2021\\15_human_rnaseq\\four_sampleXX.csv')

  # try shorter
T_24_filtered <- T_24[rowSums(T_24>0)>5,colSums(T_24>1)>1000]
K_84_filtered <- K_84[rowSums(K_84>0)>5,colSums(K_84>1)>1000]
colnames(T_24_filtered) <- seq(1,234)
colnames(K_84_filtered) <- seq(235,301)

matrix <- merge(T_24_filtered, K_84_filtered, by=0)
row.names(matrix) <- matrix[,1]
matrix <- matrix[,-1]

condition <- factor(c(rep("y", 234), rep("a", 67)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\15_human_rnaseq\\DEGseq_result.csv')
write.csv(matrix, file='D:\\SRTP2021\\15_human_rnaseq\\K_84_T_24.csv')

# rat revisited
construct_matrix_rat <- function(name){
  mm <- readMM(paste('D:/SRTP2021/22_rat_revisited/',name,'_matrix.mtx', sep = ''))
  matrix <- as.matrix(mm)
  barcode <- read.table(paste('D:/SRTP2021/22_rat_revisited/',name,'_barcodes.tsv', sep = ''))
  genes <- read.table(paste('D:/SRTP2021/22_rat_revisited/',name,'_genes.tsv', sep = ''))
  df <- as.data.frame(matrix)
  row.names(df) <- genes$V1
  # colnames(df) <- barcode$V1
  return(df)
}

  # liver
MY <- construct_matrix_rat("Liver-M-Y")
MO <- construct_matrix_rat("Liver-M-O")
FY <- construct_matrix_rat("Liver-F-Y")
FO <- construct_matrix_rat("Liver-F-O")

filter_matrix <- function(df) {
  return(df[rowSums(df>0)>100,colSums(df>1)>500])
}

MY <- filter_matrix(MY)
MO <- filter_matrix(MO)
FY <- filter_matrix(FY)
FO <- filter_matrix(FO)

aged <- merge(MO, FO, by=0)
aged_size <- length(colnames(aged))-1
colnames(aged) <- c('name', seq(1,aged_size))

young <- merge(MY, FY, by=0)
young_size <- length(colnames(young))-1
colnames(young) <- c('name', seq(aged_size+1,aged_size + young_size))

matrix <- merge(aged, young, by='name')
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]
  
condition <- factor(c(rep("a", aged_size), rep("y", young_size)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\22_rat_revisited\\DEGseq_result_liverX.csv')
write.csv(matrix, file='D:\\SRTP2021\\22_rat_revisited\\liver_matrix.csv')

cyba <- matrix[row.names(matrix)=="ENSRNOG00000013014",]
cyba_t <- as.data.frame(t(cyba))
cyba_t$aged <- as.vector(c(rep("a",4778), rep("y",4704)))

# matrix to plot df
matrix_to_plotdf(matrix, "G_liver_8")
  
  # heart
MY <- construct_matrix_rat("Aorta-M-Y")
MO <- construct_matrix_rat("Aorta-M-O")
FY <- construct_matrix_rat("Aorta-F-Y")
FO <- construct_matrix_rat("Aorta-F-O")

MY <- filter_matrix(MY)
MO <- filter_matrix(MO)
FY <- filter_matrix(FY)
FO <- filter_matrix(FO)

aged <- merge(MO, FO, by=0)
aged_size <- length(colnames(aged))-1
colnames(aged) <- c('name', seq(1,aged_size))

young <- merge(MY, FY, by=0)
young_size <- length(colnames(young))-1
colnames(young) <- c('name', seq(aged_size+1,aged_size + young_size))

matrix <- merge(aged, young, by='name')
rownames(matrix) <- matrix[,1]
matrix <- matrix[,-1]

condition <- factor(c(rep("a", aged_size), rep("y", young_size)), levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)

# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\22_rat_revisited\\DEGseq_result_heart.csv')
write.csv(matrix, file='D:\\SRTP2021\\22_rat_revisited\\heart_matrix.csv')

cyba <- matrix[row.names(matrix)=="ENSRNOG00000013014",]
cyba_t <- as.data.frame(t(cyba))
cyba_t$aged <- as.vector(c(rep("a",4778), rep("y",4704)))

# subset markers heart
matrix_to_plotdf <- function(matrix, filename) {
  G_names <- read.csv('D:/SRTP2021/22_rat_revisited/11_genes_name_conversion.csv')
  matrix$ID <- row.names(matrix)
  G_heart_8 <- matrix %>% filter(ID %in% G_names$converted_alias)
  row.names(G_heart_8) <- c("Cdkn1a", "Ly6e", "Tp53", "Cyba", "Ifitm2", "Gapdh", 
                            "Icam1","Anxa3","Tspan8","Vista","Fxyd5")
  G_heart_8 <- G_heart_8[,which(colnames(G_heart_8)!="ID")]
  G_heart_8 <- as.data.frame(t(G_heart_8))
  G_heart_8$target <- c(rep("aged", aged_size), rep("young", young_size))
  write.csv(G_heart_8, file=paste('D:\\SRTP2021\\22_rat_revisited\\',filename,'.csv', sep=''))
}
matrix_to_plotdf(matrix, 'heart_11_markers')



# Monkey artery
matrix <- read_table('D:/SRTP2021/23_monkey/GSE117715_Cynomolgus_monkey_aging_artery_count.txt')
matrix <- matrix[-6447,]

filter_monkey_matrix <- function(df) {
  return(df[rowSums(df>0)>100,colSums(df>1)>500])
}

matrix.fil <- filter_monkey_matrix(matrix)
  
  # separate aorta and coronary arteries
AA <- matrix.fil[,startsWith(colnames(matrix.fil), 'AA')]
AA$Gene <- matrix.fil$Gene

CA <- matrix.fil[,startsWith(colnames(matrix.fil), 'CA')]
CA$Gene <- matrix.fil$Gene

  # specify the Young & Aged samples
colnames(AA) <- gsub('.+Y.+','Y',colnames(AA))
colnames(AA) <- gsub('.+O.+','O',colnames(AA))
AA <- as.data.frame(AA)
row.names(AA) <- AA$Gene
AA <- AA[,!colnames(AA) %in% c('Gene')]


colnames(CA) <- gsub('.+Y.+','Y',colnames(CA))
colnames(CA) <- gsub('.+O.+','O',colnames(CA))
CA <- as.data.frame(CA)
row.names(CA) <- CA$Gene
CA <- CA[,!colnames(CA) %in% c('Gene')]

# DESeq2 for AA
matrix <- AA
tmp_aged <- str_sub(colnames(matrix),1,1)
condition <- factor(tmp_aged, 
                    levels = c("Y", "O"))
colData <- data.frame(row.names = colnames(matrix), condition)
  # DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\23_monkey\\DEGseq_result_AA.csv')

# DESeq2 for CA
matrix <- CA
tmp_aged <- str_sub(colnames(matrix),1,1)
condition <- factor(tmp_aged, 
                    levels = c("Y", "O"))
colData <- data.frame(row.names = colnames(matrix), condition)
# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\23_monkey\\DEGseq_result_CA.csv')

  # test expr
anxa3 <- CA[row.names(CA)=='ANXA3',]
y.anxa3 <- anxa3[,str_sub(colnames(anxa3),1,1)=='Y']
o.anxa3 <- anxa3[,str_sub(colnames(anxa3),1,1)=='O']
hist(as.matrix(o.anxa3[1,]))
hist(as.matrix(y.anxa3[1,]))
df.anxa <- data.frame(expr=c(as.matrix(o.anxa3[1,]), as.matrix(y.anxa3[1,])),
                      age=c(rep('o',length(as.matrix(o.anxa3[1,]))),rep('y',length(as.matrix(y.anxa3[1,])))))
ggplot(df.anxa) +
  geom_boxplot(aes(x=age, color=age, y=expr))

  # expr
markers <- toupper(c("Cdkn1a", "Ly6e", "Tp53", "Cyba", "Ifitm2", "Gapdh", 
                     "Icam1","Anxa3","Tspan8","Vsir","Fxyd5"))
matrix_to_plotdf_monkey <- function(matrix, filename, markers) {
  matrix$ID <- row.names(matrix)
  G_heart_8 <- matrix %>% filter(ID %in% markers)
  G_heart_8 <- G_heart_8[,which(colnames(G_heart_8)!="ID")]
  aged <- colnames(G_heart_8)
  G_heart_8.t <- as.data.frame(t(G_heart_8))
  G_heart_8.t$target <- aged
  write.csv(G_heart_8.t, file=paste('D:\\SRTP2021\\23_monkey\\',filename,'.csv', sep=''))
}
matrix_to_plotdf_monkey(CA, "CA_11_markers", markers)
matrix_to_plotdf_monkey(AA, "AA_11_markers", markers)

# 24 intestine
filter_matrix_intestine <- function(df) {
  return(df[rowSums(df>0)>100,colSums(df>1)>300])
}
ctrl <- read.table('D:/SRTP2021/24_intestine/matrix/c_con_matrix.tsv', header = T)
ctrl.fil <- filter_matrix_intestine(ctrl)
exp <- read.table('D:/SRTP2021/24_intestine/matrix/c_exp_matrix.tsv', header = T)
exp.fil <- filter_matrix_intestine(exp)

intestine <- merge(ctrl.fil, exp.fil, by=0)
colnames(intestine) <- c('name', rep('y',length(colnames(ctrl.fil))),
                         rep('a',length(colnames(exp.fil))))
# DESeq2
matrix <- intestine
rownames(matrix) <- matrix$name
matrix <- matrix[,-1]
tmp_aged <- str_sub(colnames(matrix),1,1)
condition <- factor(tmp_aged, 
                    levels = c("y", "a"))
colData <- data.frame(row.names = colnames(matrix), condition)
# DESeq2
dds <- DESeqDataSetFromMatrix(matrix, colData, design = ~condition)
dds <- DESeq(dds)
dds
res <- results(dds)
res
summary(res)
diff_gene_Group2 <- subset(res, padj < 0.05 & log2FoldChange > 0)
dim(diff_gene_Group2)
head(diff_gene_Group2)
write.csv(diff_gene_Group2, file='D:\\SRTP2021\\24_intestine\\DEGseq_result_intestine.csv')









### LINEAR CORRELATION WITH KNOWN MARKERS? TRY B HEART
matrix <- read.csv("D:/SRTP2021/12/heart/df_heart.csv")
matrix.less <- matrix %>%
  select(c("target", "Anxa3", "Icam1", "Ly6e", "Fxyd5", "Tspan8",
           "Cdkn1a", "Trp53", "Plaur", "Cdkn2a", "Gapdh"))
# model with age R
  # ctrl
age_gapdh <- lm(data = matrix.less, target ~ Gapdh)
summary(age_gapdh) # 0.003
  # positive
age_cdkn2a <- lm(data = matrix.less, target ~ Cdkn2a)
summary(age_cdkn2a) # 0.005
age_cdkn1a <- lm(data = matrix.less, target ~ Cdkn1a)
summary(age_cdkn1a) # 0.046
age_tp53 <- lm(data = matrix.less, target ~ Trp53)
summary(age_tp53) # -0.000
age_plaur <- lm(data = matrix.less, target ~ Plaur)
summary(age_plaur) # 0.014
  # our markers
age_icam1 <- lm(data = matrix.less, target ~ Icam1)
summary(age_icam1) # 0.016
age_anxa3 <- lm(data = matrix.less, target ~ Anxa3)
summary(age_anxa3) # 0.028
age_ly6e <- lm(data = matrix.less, target ~ Ly6e)
summary(age_ly6e) # 0.021
age_fxyd5 <- lm(data = matrix.less, target ~ Fxyd5)
summary(age_fxyd5) # 0.009

# positive ~ marker
cdkn1a_icam1<- lm(data = matrix.less, Cdkn1a ~ Icam1)
summary(cdkn1a_icam1) # 0.004
cdkn1a_anxa3<- lm(data = matrix.less, Cdkn1a ~ Anxa3)
summary(cdkn1a_anxa3) # 0.074
cdkn1a_ly6e<- lm(data = matrix.less, Cdkn1a ~ Ly6e)
summary(cdkn1a_ly6e) # 0.111
cdkn1a_icam1<- lm(data = matrix.less, Cdkn1a ~ Icam1)
summary(cdkn1a_icam1) # 0.004

corrplot::corrplot(cor(matrix.less, method = "spearman"), method = "ellipse")
