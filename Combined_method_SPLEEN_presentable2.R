library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
future::plan("multisession")

spln <- readRDS("C:/Users/kiera/OneDrive/Desktop/Ames Research/LABELLED_ERYS_AND_PHAGOS_reduced.rds")


names(splcelltypes) <- levels(spln)
spln <- RenameIdents(spln, splcelltypes)

Idents(spln) <- "Celltype"

#saveRDS(spln, "Spleen_Erys_and_Phagos.rds") 1/27


Spleen <- readRds("SPLEEN_WHOLE_LABELLED.rds")

spln1.data <- ReadMtx(
    mtx = "all_data/Spleen/GO14/filtered/matrix.mtx",
    features = "all_data/Spleen/GO14/filtered/features.tsv",
    cells = "all_data/Spleen/GO14/filtered/barcodes.tsv")

spln1 <- CreateSeuratObject(
    counts = spln1.data)
spln1
spln1[["percent.mt"]] <- PercentageFeatureSet(spln1, pattern = "^mt-")
spln1 <- subset(spln1, subset = nFeature_RNA > 10 & nFeature_RNA < 3000 & 
		percent.mt < 20 & nCount_RNA < 10000 &  nCount_RNA > 30)
spln1 <- subset(spln1, subset = nFeature_RNA < 200 & 
		(`Hba-a1` + `Hba-a2` +`Hbb-bs` +`Hbb-bt` + Hemgn < 10), invert = TRUE)
spln1


Spleen[["percent.mt"]] <- PercentageFeatureSet(Spleen, pattern = "^mt-")
Spleen <- subset(Spleen, subset = nFeature_RNA > 10 & nFeature_RNA < 3000 & 
		percent.mt < 20 & nCount_RNA < 10000 &  nCount_RNA > 30)
Spleen <- subset(spln, subset = nFeature_RNA < 200 & 
		(`Hba-a1` + `Hba-a2` +`Hbb-bs` +`Hbb-bt` + Hemgn < 10), invert = TRUE)





spln2.data <- ReadMtx(
    mtx = "all_data/Spleen/GO16/filtered/matrix.mtx",
    features = "all_data/Spleen/GO16/filtered/features.tsv",
    cells = "all_data/Spleen/GO16/filtered/barcodes.tsv")

spln2 <- CreateSeuratObject(
    counts = spln2.data)
spln2
spln2[["percent.mt"]] <- PercentageFeatureSet(spln2, pattern = "^mt-")
spln2 <- subset(spln2, subset = nFeature_RNA > 10 & nFeature_RNA < 3000 & 
		percent.mt < 20 & nCount_RNA < 10000 &  nCount_RNA > 30)
spln2 <- subset(spln2, subset = nFeature_RNA < 200 & 
		(`Hba-a1` + `Hba-a2` +`Hbb-bs` +`Hbb-bt` + Hemgn < 10), invert = TRUE)
spln2


spln3.data <- ReadMtx(
    mtx = "all_data/Spleen/GO19/filtered/matrix.mtx",
    features = "all_data/Spleen/GO19/filtered/features.tsv",
    cells = "all_data/Spleen/GO19/filtered/barcodes.tsv")

spln3 <- CreateSeuratObject(
    counts = spln3.data)
spln3
spln3[["percent.mt"]] <- PercentageFeatureSet(spln3, pattern = "^mt-")
spln3 <- subset(spln3, subset = nFeature_RNA > 10 & nFeature_RNA < 3000 & 
		percent.mt < 20 & nCount_RNA < 10000 &  nCount_RNA > 30)
spln3 <- subset(spln3, subset = nFeature_RNA < 200 & 
		(`Hba-a1` + `Hba-a2` +`Hbb-bs` +`Hbb-bt` + Hemgn < 10), invert = TRUE)
spln3


spln4.data <- ReadMtx(
    mtx = "all_data/Spleen/GO20/filtered/matrix.mtx",
    features = "all_data/Spleen/GO20/filtered/features.tsv",
    cells = "all_data/Spleen/GO20/filtered/barcodes.tsv")

spln4 <- CreateSeuratObject(
    counts = spln4.data)
spln4
spln4[["percent.mt"]] <- PercentageFeatureSet(spln4, pattern = "^mt-")
spln4 <- subset(spln4, subset = nFeature_RNA > 10 & nFeature_RNA < 3000 & 
		percent.mt < 20 & nCount_RNA < 10000 &  nCount_RNA > 30)
spln4 <- subset(spln4, subset = nFeature_RNA < 200 & 
		(`Hba-a1` + `Hba-a2` +`Hbb-bs` +`Hbb-bt` + Hemgn < 10), invert = TRUE)
spln4




spln <- readRDS(file = "SPLEEN_ALL_ERYS_UMAPED2.rds")


(VlnPlot(spln1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3) |
VlnPlot(spln2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)) /
(VlnPlot(spln3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3) | 
VlnPlot(spln4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0, ncol = 3)) 

spln1$subject <- "GO14"
spln2$subject <- "GO16"
spln3$subject <- "GO19"
spln4$subject <- "GO20"

spln <- merge(spln1, y = c(spln2, spln3, spln4))
table(spln$subject)

FeatureScatter(spln, group.by="subject",feature1 = "Hba-a1", feature2 = "Hbb-bt", raster=T)

VlnPlot(spln, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by="subject",pt.size = 0, ncol = 3)

#saveRDS(spln, file = "filtered_spln_cells.rds") 

#==========================================================================
spln[["RNA.diversity"]] <- spln$nFeature_RNA / spln$nCount_RNA
		#Values: 0.8-0.95
spln[["RNA.homogeneity"]] <- spln$nCount_RNA / spln$nFeature_RNA

spln[["percent.mt"]] <- PercentageFeatureSet(spln, pattern = "^mt-")
		#Mito genes
spln[["percent.cd"]] <- PercentageFeatureSet(spln, pattern = "^Cd[ck]")
		#Cell cycle genes, CDCs and CDKs
spln[["percent.rb"]] <- PercentageFeatureSet(spln, pattern =  "^Rb[sl]")
		#Ribo small/large segment genes

# Counts for alpha/beta Hb poduction
spln$HBax <- PercentageFeatureSet(spln, pattern = "^Hba-") * spln$nCount_RNA / 100
spln$HBbx <- PercentageFeatureSet(spln, pattern = "^Hbb-") * spln$nCount_RNA / 100

# Not v helpful
spln$BAratio <- ifelse(spln$HBbx > 0 & spln$HBax > 0, (spln$HBbx / spln$HBax), 0)
spln$BAspread <- ifelse(spln$HBbx > 0 & spln$HBax > 0, (spln$HBbx - spln$HBax), 0)



spln[["Erythroidness"]] <- PercentageFeatureSet(spln, features = c(
		"Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt",
		"Hemgn", "Ermap", "Slc25a21", "Slc4a1"))
		#Values: 0.2-40

spln[["B.cellness"]] <- PercentageFeatureSet(spln, features = c(
		"Pax5", "Ebf4", "Bank1", "Bach2",
		"Vpreb1", "Lef1", "Cecr2", "Vpreb3"))
		#Values: 0.5-1.5

spln[["T.cellness"]] <- PercentageFeatureSet(spln, features = c(
		"Ccl5", "Nkg7", "Bcl11b", "Skap1",
		"Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
		#Values: 1-2

spln[["Macrophageness"]] <- PercentageFeatureSet(spln, features = c(
  "Cd68", "Adgre1", "Itgam", "Cd14",
  "Lyz2", "Mertk", "Csf1r", "Marco"))
#Values:

spln[["Neutrophilness"]] <- PercentageFeatureSet(spln, features = c(
		"Ltf", "Ngp", "Mmp8", "Retnlg",
		"S100a8", "Ltf", "S100a9", "Elane"))
		#Values: 0-4

spln[["Monocyteness"]] <- PercentageFeatureSet(spln, features = c(
		"F13a1", "Gria3", "Prtn3", "Mpo", 
		"Crip1", "Lgals1", "Tmsb10", "Ly6a"))
		#Values: 0.5-2
#===========================================================================
spln[["RNA"]] <- JoinLayers(spln[["RNA"]])

Spleen <- spln
spln <- Spleen
################ Standard processing ###############
spln <- NormalizeData(spln)

# User input
num_variable_features <- 4000
neighbors_dims <- 20
cluster_resolution <- .4
umap_dims <- 20


spln <- FindVariableFeatures(spln, selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(spln)
spln <- ScaleData(spln, features = all.genes)#VariableFeatures(object = spln))
spln <- RunPCA(spln, features = all.genes)#VariableFeatures(object = spln))
ElbowPlot(spln, 50)

spln <- FindNeighbors(spln, dims = 1:neighbors_dims)
spln <- FindClusters(spln, resolution = cluster_resolution)
spln <- RunUMAP(spln, dims = 1:umap_dims, return.model = F)
dim <- DimPlot(spln, group.by="ident",reduction = "umap",label = T, pt.size = .75, shuffle = T, raster = T) + NoLegend()
dim2 <- DimPlot(spln, group.by="subject",reduction = "umap",label = T, pt.size = .75, shuffle = T, raster = T) + NoLegend()
dim|dim2
feats <- FeaturePlot(spln, features = c(
		"Hba-a1", "Ighm", 
		"nCount_RNA", "Macrophageness"), 
		min.cutoff = c(1,0,20,0), max.cutoff = c(100,2,1500,2), 
		ncol = 2, pt.size = 2, raster=T)
dim + feats


all_erysfeats <- FeaturePlot(spln, features = c(
  "Cd34", "EPs", "Hmgb2", "Bpgm",
  "Hemgn", "Hba-a1", "nCount_RNA", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
  min.cutoff = c(0,0,3,0,
                 0,0,0,.05),
  max.cutoff = c(1,2,5,1,
                 2,20,10000,.9))
all_erysfeats


bfeats <- FeaturePlot(spln, features = c(
  "Pax5", "Ly6c1", "Bank1", "Bach2",
  "Vpreb1", "Lef1", "Cecr2", "Vpreb3"),	 pt.size = 1.5,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(5,1,5,5,
                 1,1,1,1))
bfeats


macfeats <- FeaturePlot(spln, features = c(
  "Slc40a1", "Fcgr1", "Gfra2", "C1qa",
  "Spic", "Marco", "Hmox1", "Cd68"),	 pt.size = 1.5,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(1,1,1,1,
                 1,1,1,1))
macfeats

monocytefeats <- FeaturePlot(spln, features = c(
  "Ly6c2", "Ccr2", "Fcgr3", "Cx3cr1",
  "Csf1r", "Cd14", "Fcnb", "Lyz2"),	 pt.size = 1.5,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(1,1,1,1,
                 1,1,1,1))
monocytefeats


erycelltypes <- FeaturePlot(spln, features = c(
  "Ighm", "Bpgm", "Erythroidness", "B.cellness",
  "Hba-a1", "nCount_RNA", "Monocyteness", "Macrophageness"),	 pt.size = 1,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(7,1,20,1.5,
                 20,2000,1,1))
erycelltypes

VlnPlot(spln, features = "Macrophageness", pt.size=0)


clustmarks <- FindMarkers(spln, ident.1 = 6, ident.2=1, min.pct = 0.1, only.pos = F, logfc.threshold = .5)
topmarks <- clustmarks %>% slice_max(n = 40, order_by = avg_log2FC)
head(clustmarks, n = 40)


erys <- subset(spln, idents = c(1,4,6))
erys <- subset(erys, subset = Ighm == 0)



num_variable_features <- 4000
neighbors_dims <- 20
cluster_resolution <- .4
umap_dims <- 20
erys <- FindVariableFeatures(erys, selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(erys)
erys <- ScaleData(erys, features = all.genes)#VariableFeatures(object = erys))
erys <- RunPCA(erys, features = all.genes)#VariableFeatures(object = erys))
ElbowPlot(erys, 50)

erys <- FindNeighbors(erys, dims = 1:neighbors_dims)
erys <- FindClusters(erys, resolution = cluster_resolution)
erys <- RunUMAP(erys, dims = 1:umap_dims, return.model = F)



dim <- DimPlot(erys, group.by="ident",reduction = "umap",label = T, pt.size = .75, shuffle = T, raster = T) + NoLegend()

feats <- FeaturePlot(erys, features = c(
  "Hba-a1", "B.cellness", 
  "nCount_RNA", "Macrophageness"), 
  min.cutoff = c(1,0,20,0), max.cutoff = c(100,2,1500,2), 
  ncol = 2, pt.size = 2, raster=T)
dim + feats

#saveRDS(spln, file = "SPLEEN_WHOLE_LABELLED.rds")
clustmarks <- FindMarkers(spln, ident.1 = 9, min.pct = 0.1, only.pos = F, logfc.threshold = .5)
topmarks <- clustmarks %>% slice_max(n = 40, order_by = avg_log2FC)
head(topmarks, n = 40)



spln$Celltype <- ifelse(Idents(spln) == 13 | Idents(spln) == 22 | Idents(spln) == 24 | Idents(spln) == 37 | Idents(spln) == 30, "Erythroid Cells", "Other")


spln$Celltype <- ifelse(spln$Celltype == "Other" & spln@assays$RNA$data["Pf4",] > 1.9, "Platelets", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 19, "Plasma cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 8 | Idents(spln) == 49 | Idents(spln) == 53 | Idents(spln) == 40 , "Plasmacytoid Dendritic Cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 41 | Idents(spln) == 4 | Idents(spln) == 33 | Idents(spln) == 32 | Idents(spln) == 5 | Idents(spln) == 50 | Idents(spln) == 46 | Idents(spln) == 52 | Idents(spln) == 28, "Neutrophils",  spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 0 | Idents(spln) == 12 | Idents(spln) == 18 | Idents(spln) == 15 | Idents(spln) == 38 | Idents(spln) == 42 | Idents(spln) == 6 | Idents(spln) == 21 | Idents(spln) == 48 | Idents(spln) == 2 | Idents(spln) == 51 | Idents(spln) == 47 | Idents(spln) == 40  | Idents(spln) == 20 | Idents(spln) == 44 | Idents(spln) == 16 | Idents(spln) == 27 | Idents(spln) == 29, "Mononuclear Phagocytes", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 1 | Idents(spln) == 0 | Idents(spln) == 3 | Idents(spln) == 25 | Idents(spln) == 26 | Idents(spln) == 36 | Idents(spln) == 23 | Idents(spln) == 14, "B Cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 31 | Idents(spln) == 7 | Idents(spln) == 9 | Idents(spln) == 39 | Idents(spln) == 43 | Idents(spln) == 10 | Idents(spln) == 17 | Idents(spln) == 11 | Idents(spln) == 35, "T cells", spln$Celltype)



all_erys$Celltype <- ifelse(Idents(all_erys) == 21 | Idents(all_erys) == 17 | Idents(all_erys) == 15 | Idents(all_erys) == 6 | Idents(all_erys) == 2 | Idents(all_erys) == 0, "Phagocytic B Cells", "Other")
all_erys$Celltype <- ifelse(Idents(all_erys) == 5, "Macrophages", all_erys$Celltype)
all_erys$Celltype <- ifelse(Idents(all_erys) == 9, "High-RNA Phagocyte", all_erys$Celltype)
all_erys$Celltype <- ifelse(Idents(all_erys) == 4 | Idents(all_erys) == 10 | Idents(all_erys) == 19 | Idents(all_erys) == 22, "Bpgm Type Phagocytes", all_erys$Celltype)
all_erys$Celltype <- ifelse(all_erys$Celltype == "Other" & all_erys@assays$RNA$data["Ighm",] == 0, "Erythroid Cells", all_erys$Celltype)



spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 19, "Plasma cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 8 | Idents(spln) == 49 | Idents(spln) == 53 | Idents(spln) == 40 , "Plasmacytoid Dendritic Cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 41 | Idents(spln) == 4 | Idents(spln) == 33 | Idents(spln) == 32 | Idents(spln) == 5 | Idents(spln) == 50 | Idents(spln) == 46 | Idents(spln) == 52 | Idents(spln) == 28, "Neutrophils",  spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 0 | Idents(spln) == 12 | Idents(spln) == 18 | Idents(spln) == 15 | Idents(spln) == 38 | Idents(spln) == 42 | Idents(spln) == 6 | Idents(spln) == 21 | Idents(spln) == 48 | Idents(spln) == 2 | Idents(spln) == 51 | Idents(spln) == 47 | Idents(spln) == 40  | Idents(spln) == 20 | Idents(spln) == 44 | Idents(spln) == 16 | Idents(spln) == 27 | Idents(spln) == 29, "Mononuclear Phagocytes", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 1 | Idents(spln) == 0 | Idents(spln) == 3 | Idents(spln) == 25 | Idents(spln) == 26 | Idents(spln) == 36 | Idents(spln) == 23 | Idents(spln) == 14, "B Cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & Idents(spln) == 31 | Idents(spln) == 7 | Idents(spln) == 9 | Idents(spln) == 39 | Idents(spln) == 43 | Idents(spln) == 10 | Idents(spln) == 17 | Idents(spln) == 11 | Idents(spln) == 35, "T cells", spln$Celltype)
spln$Celltype <- ifelse(spln$Celltype == "Other" & spln@assays$RNA$data["Pf4",] > 1.9, "Platelets", spln$Celltype)


dim | DimPlot(spln, group.by="Celltype",label = T, pt.size = 1, shuffle = T, raster = T)














# Identitify clusters of interest and subset
########## Subset erythroid cells ###############


clustmarks <- FindMarkers(spln, ident.1 = 8, min.pct = 0.5, only.pos = F, logfc.threshold = .5)
head(clustmarks, n = 40)


# Identitify clusters of interest and subset
all_erys <- subset(spln, idents = c(48,52,33,57,10,22,38,47,39,31,51,35,26,56,27,23,16))
all_erys

########## Reprocess Isolated erythroid cells ############

# User input
num_variable_features_erys <- 750
neighbors_dims_erys <- 8
cluster_resolution_erys <- .6
umap_dims_erys <- 8


all_erys <- FindVariableFeatures(all_erys, selection.method = "vst", nfeatures = num_variable_features_erys)
all_erys.genes <- rownames(all_erys)
all_erys <- ScaleData(all_erys, features = VariableFeatures(all_erys))
all_erys <- RunPCA(all_erys, features = VariableFeatures(all_erys))
ElbowPlot(all_erys, 50)

all_erys <- FindNeighbors(all_erys, dims = 1:neighbors_dims_erys)
all_erys <- FindClusters(all_erys, resolution = cluster_resolution_erys)
all_erys <- RunUMAP(all_erys, dims = 1:umap_dims_erys, reduction = "pca", return.model = T)

all_erysdim <- DimPlot(all_erys,raster=T,group.by="ident",reduction = "umap", repel = F, label = T, shuffle = T, pt.size = .5) + NoLegend()
all_erysdim2 <- DimPlot(all_erys,raster=T,group.by="cohort",reduction = "umap", repel = T, label = T, shuffle = T, pt.size = .5) + NoLegend()
all_erysdim|all_erysdim2
all_erysfeats <- FeaturePlot(all_erys, features = c(
		, "Hba-a1", "nCount_RNA", "Tfrc",
		"Bpgm", "Ighm", "Gypa", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
		min.cutoff = c(0,0,0,0,
				0,0,0,.05),
		max.cutoff = c(1,15,1000,1,
				1,1,1,.9))
all_erysfeats


########## Compare handpicked to computed celltypes and clean up dataset ##############

saveRDS(all_erys, file = "SPLEEN_FO_ery_cells.rds")


table(all_erys$handpicked_celltype,Idents(all_erys))


clustmarks <- FindMarkers(all_erys, ident.1 = 4,ident.2=3, min.pct = 0.5, only.pos = F, 
	logfc.threshold = .2)
head(clustmarks, n = 40)


all_markers <- FindAllMarkers(Spleen, min.pct = 0.2, logfc.threshold = .2)
write.csv(all_markers, "WHOLE_SPLEEN_MARKERS.csv")



topmarks <- all_erys.markers %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC)

heat = DoHeatmap(all_erys, features = topmarks$gene, size = 4, angle = 10) + NoLegend()
all_erysdim | heat

all_erys <- subset(all_erys, idents = c(1,4,7), invert =T)

saveRDS(all_erys, "SPLEEN_ALL_ERYS_UMAPED2.rds")









###
