# Cell Subgroups Analysis

# Import Packages
library(tidyverse)
library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(SingleR)
library(celldex)
# library(harmony)

## Import Datasets & Combine
# B Heart
  # Tabula Muris Senis
B.young <- read.csv('D:/SRTP2021/12/heart/young_t.csv')
row.names(B.young) <- B.young$index
B.young <- B.young[,-1]
B.young.seurat <- raw_to_seurat(B.young, 'B.heart.y', 'y')
B.old <- read.csv('D:/SRTP2021/12/heart/old_t.csv')
row.names(B.old) <- B.old$index
B.old <- B.old[,-1]
B.old.seurat <- raw_to_seurat(B.old, 'B.heart.o', 'o')
  # combine by anchors
B.heart.anchors <- FindIntegrationAnchors(object.list = list(B.young.seurat, B.old.seurat), dims = 1:20)
B.heart <- IntegrateData(anchorset = B.heart.anchors, dims = 1:20)
# saveRDS(B.heart, 'D:/SRTP2021/12/heart/B_heart.rds')
B.heart <- readRDS('D:/SRTP2021/12/heart/B_heart.rds')


# B Liver
B.young <- read.csv('D:/SRTP2021/12/liver/out/young_t.csv')
row.names(B.young) <- B.young$index
B.young <- B.young[,-1]
B.young.seurat <- raw_to_seurat(B.young, 'B.liver.y', 'y')
B.old <- read.csv('D:/SRTP2021/12/liver/out/old_t.csv')
row.names(B.old) <- B.old$index
B.old <- B.old[,-1]
B.old.seurat <- raw_to_seurat(B.old, 'B.liver.o', 'o')
# combine by anchors
B.liver.anchors <- FindIntegrationAnchors(object.list = list(B.young.seurat, B.old.seurat), 
                                          dims = 1:20)
B.liver <- IntegrateData(anchorset = B.liver.anchors, dims = 1:20)
# saveRDS(B.liver, 'D:/SRTP2021/12/liver/B_liver.rds')
B.liver <- readRDS('D:/SRTP2021/12/liver/B_liver.rds')

# B Lung
B.young <- read.csv('D:/SRTP2021/12/lung/out/young_t.csv')
row.names(B.young) <- B.young$index
B.young <- B.young[,-1]
B.young.seurat <- raw_to_seurat(B.young, 'B.lung.y', 'y')
B.old <- read.csv('D:/SRTP2021/12/lung/out/old_t.csv')
row.names(B.old) <- B.old$index
B.old <- B.old[,-1]
B.old.seurat <- raw_to_seurat(B.old, 'B.lung.o', 'o')
# combine by anchors
B.lung.anchors <- FindIntegrationAnchors(object.list = list(B.young.seurat, B.old.seurat), 
                                          dims = 1:20)
B.lung <- IntegrateData(anchorset = B.lung.anchors, dims = 1:20)
saveRDS(B.lung, 'D:/SRTP2021/12/lung/B_lung.rds')
B.lung <- readRDS('D:/SRTP2021/12/lung/B_lung.rds')

# G Aorta
dir_to_seurat <- function(dir, dataname, age) {
  data <- Read10X(dir)
  data <- CreateSeuratObject(counts = data, project = dataname, min.cells = 5)
  data$name <- dataname
  data$age <- age
  data <- subset(data, subset = nFeature_RNA > 500)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  return(data)
}
G.amo <- dir_to_seurat('D:/SRTP2021/22_rat_revisited/22_sub/AMO/', 'G.amo', 'o')
G.afo <- dir_to_seurat('D:/SRTP2021/22_rat_revisited/22_sub/AFO/', 'G.afo', 'o')
G.afy <- dir_to_seurat('D:/SRTP2021/22_rat_revisited/22_sub/AFY/', 'G.afy', 'y')
G.amy <- dir_to_seurat('D:/SRTP2021/22_rat_revisited/22_sub/AMY/', 'G.amy', 'y')
  # combine by anchors
G.anchors <- FindIntegrationAnchors(object.list = list(G.amo, G.afo, G.afy, G.amy), 
                                    dims = 1:20)
G.combined <- IntegrateData(anchorset = G.anchors, dims = 1:20)
# saveRDS(G.combined, 'D:/SRTP2021/22_rat_revisited/G_sub/G_combined.rds')
G.combined <- readRDS('D:/SRTP2021/22_rat_revisited/G_sub/G_combined.rds')


# Monkey
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
AA.o <- AA[,startsWith(colnames(AA),'O')]
AA.y <- AA[,startsWith(colnames(AA),'Y')]

colnames(CA) <- gsub('.+Y.+','Y',colnames(CA))
colnames(CA) <- gsub('.+O.+','O',colnames(CA))
CA <- as.data.frame(CA)
row.names(CA) <- CA$Gene
CA <- CA[,!colnames(CA) %in% c('Gene')]
CA.o <- CA[,startsWith(colnames(CA),'O')]
CA.y <- CA[,startsWith(colnames(CA),'Y')]
  # create seurat objects
raw_to_seurat <- function(data, dataname, age) {
  data <- CreateSeuratObject(counts = data, project = dataname, min.cells = 5)
  data$name <- dataname
  data$age <- age
  data <- subset(data, subset = nFeature_RNA > 500)
  data <- NormalizeData(data, verbose = FALSE)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  return(data)
}
ca.o <- raw_to_seurat(CA.o, 'CA.o', 'o')
ca.y <- raw_to_seurat(CA.y, 'CA.y', 'y')
aa.o <- raw_to_seurat(AA.o, 'AA.o', 'o')
aa.y <- raw_to_seurat(AA.y, 'AA.y', 'y')

monkey.anchors <- FindIntegrationAnchors(object.list = list(ca.o, ca.y, aa.o, aa.y), 
                                         dims = 1:20)
monkey.combined <- IntegrateData(anchorset = monkey.anchors, dims = 1:20)
# saveRDS(monkey.combined, 'D:/SRTP2021/23_monkey/monkey_combined.rds')
monkey.combined <- readRDS('D:/SRTP2021/23_monkey/monkey_combined.rds')

# intestine
con <- read.table('D:/SRTP2021/24_intestine/matrix/c_con_matrix.tsv',header = T)
exp <- read.table('D:/SRTP2021/24_intestine/matrix/c_exp_matrix.tsv',header = T)
filter_matrix_intestine <- function(df) {
  return(df[rowSums(df>0)>100,colSums(df>1)>300])
}
con.fil <- filter_matrix_intestine(con)
exp.fil <- filter_matrix_intestine(exp)
con.fil.obj <- raw_to_seurat(con.fil,'intestine.y','y')
exp.fil.obj <- raw_to_seurat(exp.fil,'intestine.o','o')

intestine.anchors <- FindIntegrationAnchors(object.list = list(con.fil.obj, exp.fil.obj), 
                                            dims = 1:20)
intestine.combined <- IntegrateData(anchorset = intestine.anchors, dims = 1:20)
# saveRDS(intestine.combined, 'D:/SRTP2021/24_intestine/intestine_combined.rds')
intestine.combined <- readRDS('D:/SRTP2021/24_intestine/intestine_combined.rds')


# Additional Young Data from Tabula Muris
  # read in Robj
old_seurat_to_integrated_new <- function(path, name, age) {
  load(path)
  tiss.new <- UpdateSeuratObject(tiss)
  tiss.new$name <- name
  tiss.new$age <- age
  return(tiss.new)
}

B.heart.addi <- old_seurat_to_integrated_new('D:/SRTP2021/12/heart/droplet_Heart_and_Aorta_seurat_tiss.Robj', 
                                             'B.heart.add', 'y')
B.liver.addi <- old_seurat_to_integrated_new('D:/SRTP2021/12/liver/droplet_Liver_seurat_tiss.Robj', 
                                             'B.liver.add', 'y')
B.lung.addi <- old_seurat_to_integrated_new('D:/SRTP2021/12/lung/droplet_Lung_seurat_tiss.Robj', 
                                             'B.lung.add', 'y')
  # integration
integrate_datasets <- function(data.list) {
  data.anchors <- FindIntegrationAnchors(object.list = data.list, dims = 1:20)
  data.combined <- IntegrateData(anchorset = data.anchors, dims = 1:20)
  return(data.combined)
}
B.heart <- integrate_datasets(list(B.heart, B.heart.addi))
B.liver <- integrate_datasets(list(B.liver, B.liver.addi))
B.lung <- integrate_datasets(list(B.lung, B.lung.addi))
B.liver.heart <- integrate_datasets(list(B.heart, B.liver))
B.liver.heart.analysis <- integrated_analysis(B.liver.heart)


# Integrated Analysis
integrated_analysis <- function(data) {
  #DefaultAssay(data) <- "integrated"
  # Run the standard workflow for visualization & clustering
  data <- ScaleData(data, verbose = F)
  data <- RunPCA(data, npcs = 30, verbose = F)
  # t-SNE and Clustering
  data <- RunUMAP(data, reduction = "pca", dims=1:20)
  data <- FindNeighbors(data, reduction = "pca", dims = 1:20)
  data <- FindClusters(data, resolution = 0.5)
  return(data)
}
  # G plots
G.combined <- integrated_analysis(G.combined)
DimPlot(G.combined, reduction = "umap", group.by = "name", pt.size=0.5) +
  ggtitle("G Rat Aorta UMAP (Sample)")
DimPlot(G.combined, reduction = "umap", group.by = "age", pt.size=0.5) +
  ggtitle("G Rat Aorta UMAP (Aged)")
DimPlot(G.combined, reduction = "umap", label = TRUE, pt.size=0.5) +
  ggtitle("G Rat Aorta UMAP (Cluster)")
DimPlot(G.combined, reduction = "umap", split.by = "name", pt.size=0.5, ncol=2) +
  ggtitle("G Rat Aorta UMAP (Sample VS Cluster)")

  # B Heart plots
B.heart.analysis <- integrated_analysis(B.heart)
DimPlot(B.heart.analysis, reduction = "umap", group.by = "name", pt.size=0.5) +
  ggtitle("B Heart UMAP (Sample)")
DimPlot(B.heart.analysis, reduction = "umap", group.by = "age", pt.size=0.5) +
  ggtitle("B Heart UMAP (Aged)")
DimPlot(B.heart.analysis, reduction = "umap", label = TRUE, pt.size=0.5) +
  ggtitle("B Heart UMAP (Cluster)")
DimPlot(B.heart.analysis, reduction = "umap", split.by = "name", pt.size=0.5, ncol=2) +
  ggtitle("B Heart UMAP (Sample VS Cluster)")

  # B Liver plots
B.liver.analysis <- integrated_analysis(B.liver)
DimPlot(B.liver.analysis, reduction = "umap", group.by = "name", pt.size=0.5) +
  ggtitle("B Liver UMAP (Sample)")
DimPlot(B.liver.analysis, reduction = "umap", group.by = "age", pt.size=0.5) +
  ggtitle("B Liver UMAP (Aged)")
DimPlot(B.liver.analysis, reduction = "umap", label = TRUE, pt.size=0.5) +
  ggtitle("B Liver UMAP (Cluster)")
DimPlot(B.liver.analysis, reduction = "umap", split.by = "name", pt.size=0.5, ncol=1) +
  ggtitle("B Liver UMAP (Sample VS Cluster)")

  # monkey plots
monkey.combined <- integrated_analysis(monkey.combined)
DimPlot(monkey.combined, reduction = "umap", group.by = "name", pt.size=0.5) +
  ggtitle("Monkey Arteries UMAP (Sample)")
DimPlot(monkey.combined, reduction = "umap", group.by = "age", pt.size=0.5) +
  ggtitle("Monkey Arteries UMAP (Aged)")
DimPlot(monkey.combined, reduction = "umap", label = TRUE, pt.size=0.5) +
  ggtitle("Monkey Arteries UMAP (Cluster)")
DimPlot(monkey.combined, reduction = "umap", split.by = "name", pt.size=0.5, ncol=2) +
  ggtitle("Monkey Arteries UMAP (Sample VS Cluster)")

# intestine
intestine.analysis <- integrated_analysis(intestine.combined)
DimPlot(intestine.analysis, reduction = "umap", group.by = "name", pt.size=0.5) +
  ggtitle("Intestine UMAP (Sample)")
DimPlot(intestine.analysis, reduction = "umap", group.by = "age", pt.size=0.5) +
  ggtitle("Intestine UMAP (Aged)")
DimPlot(intestine.analysis, reduction = "umap", label = TRUE, pt.size=0.5) +
  ggtitle("Intestine UMAP (Cluster)")
DimPlot(intestine.analysis, reduction = "umap", split.by = "name", pt.size=0.5, ncol=2) +
  ggtitle("Intestine UMAP (Sample VS Cluster)")


# Identify conserved cell type markers
best_markers_cluster <- function(data, start, end) {
  DefaultAssay(data) <- "RNA"
  best_markers <- c()
  for (i in start:end) {
    current_markers <- FindConservedMarkers(data, ident.1 = i, 
                                            grouping.var = "name", verbose = FALSE)
    best_marker <- row.names(current_markers)[1]
    best_markers <- c(best_markers, best_marker)
    print(paste(best_marker,": ",i,"/",(end)," markers found...",sep = ''))
  }
  return(best_markers)
}

  # G markers plot
# G.markers <- best_markers_cluster(G.combined, 0, 13)
G.markers <- c('Alpl','Nov','Cadm3','Myl12a','Cybb','Myh11','Emcn','Tfap2b',
                   'Cd3g','Itgb7','Kif20b','Cd163','Gja5','Cst6')
DefaultAssay(G.combined) <- "RNA"
FeaturePlot(G.combined, features = G.markers, min.cutoff = "q9", ncol = 3)

    # Violin Plot
G.markers.violin <- c('Alpl','Nov','Cadm3','Myl12a','Cybb','Myh11','Emcn',
                      'Cd3g','Itgb7','Kif20b','Cd163','Cst6','Anxa3','Icam1')
VlnPlot(G.combined, features = G.markers.violin, stack = T, flip = T) +
  ggtitle("G Rat Marker Expression")

  # B markers plot
# B.markers <- best_markers_cluster(B.heart.analysis, start = 17, end = 17)
B.markers <- c('Gpihbp1','Cd68','Mfap4','Itga6','Smoc2','Gpm6a','Tbx20',
               'Ptprcap','Ltbp2','Alox12','Ryr2','Myh11','Car4','Vtn',
               'Lmod2','Sox10')
DefaultAssay(B.heart.analysis) <- "RNA"
FeaturePlot(B.heart.analysis, features = B.markers, min.cutoff = "q9", ncol = 4)
    # Violin Plot
B.markers.violin <- c('Gpihbp1','Cd68','Mfap4','Itga6','Smoc2','Tbx20',
                      'Ptprcap','Ryr2','Car4','Lmod2','Anxa3','Icam1')
VlnPlot(B.heart.analysis, features = B.markers.violin, stack = T, flip = T) +
  ggtitle("B Heart Marker Expression")

  # monkey markers plot
DefaultAssay(monkey.combined) <- "RNA"
# monkey.markers <- best_markers_cluster(monkey.combined, start = 11, end = 13)
monkey.markers <- c('BMP4','FCN3','MYOCD','SCARA5','PLVAP','LMOD1','ERCC-00074',
                    'GPIHBP1','C1orf162','RELN','CD6','SLC45A3') # 10 missing
FeaturePlot(monkey.combined, features = monkey.markers, min.cutoff = "q9", ncol = 3)
    # Violin Plot
monkey.markers.violin <- c('BMP4','FCN3','MYOCD','SCARA5','LMOD1','ERCC-00074',
                           'C1orf162','CD6','ANXA3','ICAM1')
VlnPlot(monkey.combined, features = monkey.markers.violin, stack = T, flip = T) +
  ggtitle("Monkey Arteries Marker Expression")



# SingleR for Cell Annotation
# references
hpca.se <- HumanPrimaryCellAtlasData()
mouse.ref <- MouseRNAseqData()

# manual function to check cell markers:
    # clusters to check:12,13
    # find all markers
data <- B.heart.analysis
current_marker <- FindConservedMarkers(data, ident.1 = 12, 
                                       verbose = FALSE, grouping.var = 'name')
FeaturePlot(data, rownames(current_marker)[9:12], min.cutoff = "q9")

# B.heart cell annotation
    # store cluster indices
smooth.index <- which(Idents(B.heart.analysis) %in% c(11,13,15))
endo.gpm6a.index <- which(Idents(B.heart.analysis) %in% c(5))
    # prediction
pred.B.heart <- SingleR(test = B.heart.analysis@assays$RNA@data, ref = mouse.ref, assay.type.test = 1,
                          labels = mouse.ref$label.main, clusters = B.heart.analysis@active.ident)
new.cluster.ids.B.heart <- pred.B.heart$pruned.labels
names(new.cluster.ids.B.heart) <- levels(B.heart.analysis)
B.heart.analysis <- RenameIdents(B.heart.analysis, new.cluster.ids.B.heart)
    # rename wrong labels
Idents(B.heart.analysis, cells=smooth.index) <- "Smooth Muscle Cell"
Idents(B.heart.analysis, cells=endo.gpm6a.index) <- "Gpm6a+ Endothelial Cell"
    # plot
DimPlot(B.heart.analysis, reduction = "umap", pt.size=0.5) +
  ggtitle("B Heart UMAP (Cell Type)")

  # B.liver cell annotation
    # prediction
pred.B.liver <- SingleR(test = B.liver.analysis@assays$RNA@data, ref = mouse.ref, assay.type.test = 1,
                        labels = mouse.ref$label.fine, clusters = B.liver.analysis@active.ident)
new.cluster.ids.B.liver <- pred.B.liver$pruned.labels
names(new.cluster.ids.B.liver) <- levels(B.liver.analysis)
B.liver.analysis <- RenameIdents(B.liver.analysis, new.cluster.ids.B.liver)
    # plot
DimPlot(B.liver.analysis, reduction = "umap", pt.size=0.5) +
  ggtitle("B Liver UMAP (Cell Type)")

  # G cell annotation
    # store cluster indices
smooth.index <- which(Idents(G.combined) %in% c(1,5))
ery.index <- which(Idents(G.combined) %in% c(10))
endo.gja5.index <- which(Idents(G.combined) %in% c(12))
endo.mmrn1.index <- which(Idents(G.combined) %in% c(13))
    # prediction
pred.G.combined <- SingleR(test = G.combined@assays$RNA@data, ref = mouse.ref, assay.type.test = 1,
                           labels = mouse.ref$label.main, clusters = G.combined@active.ident)
new.cluster.ids.G.combined <- pred.G.combined$pruned.labels
names(new.cluster.ids.G.combined) <- levels(G.combined)
G.combined <- RenameIdents(G.combined, new.cluster.ids.G.combined)
    # rename wrong labels
Idents(G.combined, cells=smooth.index) <- "Smooth Muscle Cell"
Idents(G.combined, cells=ery.index) <- "Erythroid Cell"
Idents(G.combined, cells=endo.gja5.index) <- "Gja5+ Endothelial Cell"
Idents(G.combined, cells=endo.mmrn1.index) <- "Mmrn1+ Endothelial Cell"
    # plot
DimPlot(G.combined, reduction = "umap", pt.size=0.5, label=T) +
  ggtitle("G Rat Aorta UMAP (Cell Type)")

  # monkey cell annotation
    # store cluster indices
adipo.index <- which(Idents(monkey.combined) %in% c(7))
neu.index <- which(Idents(monkey.combined) %in% c(9))
endo.ackr1.index <- which(Idents(monkey.combined) %in% c(4))
    # prediction
pred.monkey.combined <- SingleR(test = monkey.combined@assays$RNA@data, ref = hpca.se, assay.type.test = 1,
                                labels = hpca.se$label.main, clusters = monkey.combined@active.ident)
new.clusters.ids.monkey.combined <- pred.monkey.combined$pruned.labels
names(new.clusters.ids.monkey.combined) <- levels(monkey.combined)
monkey.combined <- RenameIdents(monkey.combined, new.clusters.ids.monkey.combined)
    # rename wrong labels
Idents(monkey.combined, cells=adipo.index) <- "Adipocyte"
Idents(monkey.combined, cells=neu.index) <- "Neuron"
Idents(monkey.combined, cells=endo.ackr1.index) <- "ACKR1+ Endothelial Cell"
    # plot
DimPlot(monkey.combined, reduction = "umap", pt.size=0.5, label = T) +
  ggtitle("Monkey Arteries UMAP (Cell Type)")



# Membrane Marker Visualization
cart.markers <- c('Icam1','Anxa3','Tspan8','Vsir','Ly6e','Fxyd5', 'Plaur')
  # B heart
    # localization
DefaultAssay(B.heart.analysis) <- "RNA"
FeaturePlot(B.heart.analysis, features = cart.markers, min.cutoff = "q9", ncol = 3)
    # upregulation
VlnPlot(B.heart.analysis, features = cart.markers, stack = T, flip = T, split.by = "age") +
  ggtitle("B Heart Marker Expression")

  # B liver
cart.markers.more <- c(cart.markers, "Ifitm2", "Cyba", "Ier3", "Capn2", "Gbp2", "Arl6ip5", 
                       "Arl4c", "Itm2b", "Ly6e", "H2-eb1", "Fxyd5", "Ifitm3", "Cd27", "Icam2", "Rpsa", "Traf3ip3", "S100a6", "Cd69", "Hcst", "Slc14a1", "Cd52")
    # localization
DefaultAssay(B.liver.analysis) <- "RNA"
FeaturePlot(B.liver.analysis, features = cart.markers.more, min.cutoff = "q9", ncol = 3)
    # upregulation
VlnPlot(B.liver.analysis, features = cart.markers.more, stack = T, flip = T, split.by = "age") +
  ggtitle("B Liver Marker Expression")

  # G Aorta
    # localization
DefaultAssay(G.combined) <- "RNA"
FeaturePlot(G.combined, features = cart.markers, min.cutoff = "q9", ncol = 3)
    # upregulation
VlnPlot(G.combined, features = cart.markers, stack = T, flip = T, split.by = "age") +
  ggtitle("G Aorta Marker Expression")

  # Monkey Arteries
    # localization
DefaultAssay(monkey.combined) <- "RNA"
FeaturePlot(monkey.combined, features = toupper(cart.markers), min.cutoff = "q9", ncol = 3)
    # upregulation
VlnPlot(monkey.combined, features = toupper(cart.markers), stack = T, flip = T, split.by = "age") +
  ggtitle("Monkey Arteries Marker Expression")




