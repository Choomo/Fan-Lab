library(dplyr)
library(Seurat)
library(patchwork)
library("DoubletFinder")
setwd("/share/home/zhuxm/liaochh/data/Matrix/Read10X/A/")
Afolders <- list.files("/share/home/zhuxm/liaochh/data/Matrix/Read10X/A/")

#Load AB dataset
OrganoidA.list <- lapply(Afolders, function(Afolder){
  CreateSeuratObject(counts = Read10X(Afolder), 
                     project = Afolder)
})
Atime <- c("D104","D110","D205","D45","D60","D60", "D90")
Aorig <- c("A-D104","A-D110","A-D205","A-D45","A-D60","A-D60", "A-D90")
for (i in c(1:7)){
  OrganoidA.list[[i]] <- AddMetaData(object = OrganoidA.list[[i]], metadata = Atime[i], col.name = "Time")
  OrganoidA.list[[i]] <- AddMetaData(object = OrganoidA.list[[i]], metadata = Aorig[i], col.name = "orig")
}
for (i in c(1:7)){
  levels(OrganoidA.list[[i]]@meta.data[["orig.ident"]]) <- "A"
}
for (i in c(1:7)){
  OrganoidA.list[[i]][["percent.mt"]] <- PercentageFeatureSet(OrganoidA.list[[i]], pattern = "^MT-")
}
OrganoidA.list[[1]] <- subset(OrganoidA.list[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)
OrganoidA.list[[2]] <- subset(OrganoidA.list[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
OrganoidA.list[[3]] <- subset(OrganoidA.list[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)
OrganoidA.list[[4]] <- subset(OrganoidA.list[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 7.5)
OrganoidA.list[[5]] <- subset(OrganoidA.list[[5]], subset = nFeature_RNA > 200 & nFeature_RNA < 1000 & percent.mt < 10)
OrganoidA.list[[6]] <- subset(OrganoidA.list[[6]], subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 5)
OrganoidA.list[[7]] <- subset(OrganoidA.list[[7]], subset = nFeature_RNA > 200 & nFeature_RNA < 500 & percent.mt < 8)

for (i in c(1:7)){
  OrganoidA.list[[i]] <- NormalizeData(OrganoidA.list[[i]])
  OrganoidA.list[[i]] <- FindVariableFeatures(OrganoidA.list[[i]], selection.method = "vst", nfeatures = 2000)
  OrganoidA.list[[i]] <- ScaleData(OrganoidA.list[[i]])
  OrganoidA.list[[i]] <- RunPCA(OrganoidA.list[[i]])
  OrganoidA.list[[i]] <- RunUMAP(OrganoidA.list[[i]], dims = 1:15)
  OrganoidA.list[[i]] <- FindNeighbors(OrganoidA.list[[i]], dims = 1:15)
  OrganoidA.list[[i]] <- FindClusters(OrganoidA.list[[i]], resolution = 0.5)
}

OrganoidA.list <- lapply( X = OrganoidA.list, FUN = function(x){
  sweep.res.list_organoidA <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
  sweep.stats_organoidA <- summarizeSweep(sweep.res.list_organoidA, GT = FALSE)
  bcmvn_organoid <- find.pK(sweep.stats_organoidA)
  homotypic.prop <- modelHomotypic(x$seurat_clusters)
  DoubletRate = ncol(x)*8*1e-6
  nExp_poi <- round(DoubletRate*nrow(x@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
})

OrganoidA.list[[1]] <- subset(OrganoidA.list[[1]], DF.classifications_0.25_0.09_5 == "Singlet")
OrganoidA.list[[2]] <- subset(OrganoidA.list[[2]], DF.classifications_0.25_0.09_30 == "Singlet")
OrganoidA.list[[3]] <- subset(OrganoidA.list[[3]], DF.classifications_0.25_0.09_283 == "Singlet")
OrganoidA.list[[4]] <- subset(OrganoidA.list[[4]], DF.classifications_0.25_0.09_36 == "Singlet")
OrganoidA.list[[5]] <- subset(OrganoidA.list[[5]], DF.classifications_0.25_0.09_27 == "Singlet")
OrganoidA.list[[6]] <- subset(OrganoidA.list[[6]], DF.classifications_0.25_0.09_66 == "Singlet")
OrganoidA.list[[7]] <- subset(OrganoidA.list[[7]], DF.classifications_0.25_0.09_3 == "Singlet")

load("/share/home/zhuxm/liaochh/data/Matrix/Read10X/B/OrganoidB.RData")
Borig <- c("B-D24", "B-D30","B-D42", "B-D59")
for (i in c(1:4)){
  OrganoidB[[i]] <- AddMetaData(object = OrganoidB[[i]], metadata = Borig[i], col.name = "orig")
}
for (i in c(1:4)){
  OrganoidB[[i]][["percent.mt"]] <- PercentageFeatureSet(OrganoidB[[i]], pattern = "^MT-")
}
OrganoidB[[1]] <- subset(OrganoidB[[1]], subset = nFeature_RNA > 200 & nFeature_RNA < 3300 & percent.mt < 5)
OrganoidB[[2]] <- subset(OrganoidB[[2]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 3)
OrganoidB[[3]] <- subset(OrganoidB[[3]], subset = nFeature_RNA > 200 & nFeature_RNA < 3500 & percent.mt < 5)
OrganoidB[[4]] <- subset(OrganoidB[[4]], subset = nFeature_RNA > 200 & nFeature_RNA < 3000 & percent.mt < 4)
for (i in c(1:4)){
  OrganoidB[[i]] <- NormalizeData(OrganoidB[[i]])
  OrganoidB[[i]] <- FindVariableFeatures(OrganoidB[[i]], selection.method = "vst", nfeatures = 2000)
  OrganoidB[[i]] <- ScaleData(OrganoidB[[i]])
  OrganoidB[[i]] <- RunPCA(OrganoidB[[i]])
  OrganoidB[[i]] <- RunUMAP(OrganoidB[[i]], dims = 1:15)
  OrganoidB[[i]] <- FindNeighbors(OrganoidB[[i]], dims = 1:15)
  OrganoidB[[i]] <- FindClusters(OrganoidB[[i]], resolution = 0.5)
}

OrganoidB <- lapply( X = OrganoidB, FUN = function(x){
  sweep.res.list_organoidB <- paramSweep_v3(x, PCs = 1:15, sct = FALSE)
  sweep.stats_organoidB <- summarizeSweep(sweep.res.list_organoidB, GT = FALSE)
  bcmvn_organoid <- find.pK(sweep.stats_organoidB)
  homotypic.prop <- modelHomotypic(x$seurat_clusters)
  DoubletRate = ncol(x)*8*1e-6
  nExp_poi <- round(DoubletRate*nrow(x@meta.data))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  x <- doubletFinder_v3(x, PCs = 1:15, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
})

OrganoidB[[1]] <- subset(OrganoidB[[1]], DF.classifications_0.25_0.09_762 == "Singlet")
OrganoidB[[2]] <- subset(OrganoidB[[2]], DF.classifications_0.25_0.09_39 == "Singlet")
OrganoidB[[3]] <- subset(OrganoidB[[3]], DF.classifications_0.25_0.09_442 == "Singlet")
OrganoidB[[4]] <- subset(OrganoidB[[4]], DF.classifications_0.25_0.09_245 == "Singlet")

#¶ÁÈ¡ÆäËûDataset
OrganoidE <- CreateSeuratObject(counts = Read10X("/share/home/zhuxm/liaochh/data/Matrix/Read10X/E-M8/"), project = "E")
levels(OrganoidE@meta.data[["orig.ident"]]) <- "E"
OrganoidE <- AddMetaData(object = OrganoidE, metadata = "D240", col.name = "Time")
OrganoidE <- AddMetaData(object = OrganoidE, metadata = "E-D240", col.name = "orig")
Organoid3 <- read.table("/share/home/zhuxm/liaochh/data/Matrix/C-GSE122783/GSE122783_organoid.txt", header = TRUE,row.names = 1)
Organoid3 <- Organoid3[,-c(1,2,3,4,5)]
Organoid3 <- as.sparse(Organoid3)
Organoid3 <- CollapseSpeciesExpressionMatrix(Organoid3)
OrganoidC <- CreateSeuratObject(counts = Organoid3, project = "C")
rm(Organoid3)
OrganoidC <- AddMetaData(object = OrganoidC, metadata = "C-D25-D35", col.name = "orig")
OrganoidC <- AddMetaData(object = OrganoidC, metadata = "D25-35", col.name = "Time")
Organoid4 <- read.csv("/share/home/zhuxm/liaochh/data/Matrix/D-GSE119893/GSE119893-rawCounts.csv", header = TRUE,row.names = 1)
Organoid4 <- as.sparse(Organoid4)
OrganoidD <- CreateSeuratObject(counts = Organoid4, project = "D")
rm(Organoid4)
Dtime1 <- rep("D60", times = 578)
Dtime2 <- rep("D90", times = 661)
Dtime3 <- rep("D200",times = 737)
Dtime <- c(Dtime1, Dtime2, Dtime3)
OrganoidD <- AddMetaData(object = OrganoidD, metadata = Dtime, col.name = "Time")
OrganoidD <- AddMetaData(object = OrganoidD, metadata = "D", col.name = "orig")

OrganoidC[["percent.mt"]] <- PercentageFeatureSet(OrganoidC, pattern = "^MT-")
OrganoidD[["percent.mt"]] <- PercentageFeatureSet(OrganoidD, pattern = "^MT-")
OrganoidE[["percent.mt"]] <- PercentageFeatureSet(OrganoidE, pattern = "^MT-")
OrganoidC <- subset(OrganoidC, subset = nFeature_RNA > 200 & nFeature_RNA < 10000)
OrganoidD <- subset(OrganoidD, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)
OrganoidE <- subset(OrganoidE, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10)

Organoid.merged <- merge(OrganoidA.list[[1]],
                  y = c(OrganoidA.list[[2]],OrganoidA.list[[3]],OrganoidA.list[[4]],OrganoidA.list[[5]],OrganoidA.list[[6]],OrganoidA.list[[7]],OrganoidB[[1]],OrganoidB[[2]],OrganoidB[[3]],OrganoidB[[4]],OrganoidC, OrganoidD,OrganoidE),
                      add.cell.ids = c("A-D104","A-D110","A-D205","A-D45","A-D60a","A-D60b","A-D90","B-D24","B-D30", "B-D42","B-D59","C","D","E"),
                  project = "Organoid.merged")
rm(OrganoidA.list)
rm(OrganoidB)
rm(OrganoidC)
rm(OrganoidD)
rm(OrganoidE)
Organoid.list <- SplitObject(Organoid.merged, split.by = "orig")
rm(Organoid.merged)
Organoid.list <- lapply(X = Organoid.list, FUN = function(x){
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
features <- SelectIntegrationFeatures(object.list = Organoid.list)

#integration
Organoid.anchors <- FindIntegrationAnchors(object.list = Organoid.list, anchor.features = features)
Organoid.combined <- IntegrateData(anchorset = Organoid.anchors)
DefaultAssay(Organoid.combined) <- "integrated"
Organoid.combined <- ScaleData(Organoid.combined, verbose = TRUE)
Organoid.combined <- RunPCA(Organoid.combined, verbose = TRUE)
Organoid.combined <- FindNeighbors(Organoid.combined, dims = 1:18)
Organoid.combined <- FindClusters(Organoid.combined, resolution = 0.5)
Organoid.combined <- RunUMAP(Organoid.combined, dims = 1:18)
save(Organoid.combined, file = "/share/home/zhuxm/liaochh/data/Matrix/Organoid_combined_new3.RData")
