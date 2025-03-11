library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

blood1.data <- ReadMtx(
    mtx = "all_data/Blood/GY10/raw/matrix.mtx",
    features = "all_data/Blood/GY10/raw/features.tsv",
    cells = "all_data/Blood/GY10/raw/barcodes.tsv")

blood1 <- CreateSeuratObject(
    counts = blood1.data)


blood1
blood1[["percent.mt"]] <- PercentageFeatureSet(blood1, pattern = "^mt-")

blood1 <- subset(blood1, subset = nFeature_RNA < 2000 & 
		nCount_RNA < 4000 &  nCount_RNA > 30 & percent.mt < 20)
blood1 <- subset(blood1, subset = nCount_RNA < 200 & 
		(Bpgm + `Hba-a1` + `Hba-a2` + `Hbb-bs` +`Hbb-bt`< 100), invert = T)
blood1
VlnPlot(blood1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)


blood2.data <- ReadMtx(
    mtx = "all_data/Blood/GY4/raw/matrix.mtx",
    features = "all_data/Blood/GY4/raw/features.tsv",
    cells = "all_data/Blood/GY4/raw/barcodes.tsv")

blood2 <- CreateSeuratObject(
    counts = blood2.data)

blood2
blood2[["percent.mt"]] <- PercentageFeatureSet(blood2, pattern = "^mt-")

blood2 <- subset(blood2, subset = nFeature_RNA < 2000 & 
		nCount_RNA < 4000 &  nCount_RNA > 30 & percent.mt < 20)
blood2 <- subset(blood2, subset = nCount_RNA < 200 & 
		(Bpgm + `Hba-a1` + `Hba-a2` + `Hbb-bs` +`Hbb-bt`< 100), invert = T)
blood2


blood3.data <- ReadMtx(
    mtx = "all_data/Blood/GY6/raw/matrix.mtx",
    features = "all_data/Blood/GY6/raw/features.tsv",
    cells = "all_data/Blood/GY6/raw/barcodes.tsv")

blood3 <- CreateSeuratObject(
    counts = blood3.data)
blood3
blood3[["percent.mt"]] <- PercentageFeatureSet(blood3, pattern = "^mt-")

blood3 <- subset(blood3, subset = nFeature_RNA < 2000 & 
		nCount_RNA < 4000 &  nCount_RNA > 30 & percent.mt < 20)
blood3 <- subset(blood3, subset = nCount_RNA < 200 & 
		(Bpgm + `Hba-a1` + `Hba-a2` + `Hbb-bs` +`Hbb-bt`< 100), invert = T)
blood3



blood4.data <- ReadMtx(
    mtx = "all_data/Blood/GY9/raw/matrix.mtx",
    features = "all_data/Blood/GY9/raw/features.tsv",
    cells = "all_data/Blood/GY9/raw/barcodes.tsv")

blood4 <- CreateSeuratObject(
    counts = blood4.data)
blood4
blood4[["percent.mt"]] <- PercentageFeatureSet(blood4, pattern = "^mt-")

blood4 <- subset(blood4, subset = nFeature_RNA < 2000 & 
		nCount_RNA < 4000 &  nCount_RNA > 30 & percent.mt < 20)
blood4 <- subset(blood4, subset = nCount_RNA < 200 & 
		(Bpgm + `Hba-a1` + `Hba-a2` + `Hbb-bs` +`Hbb-bt`< 100), invert = T)
blood4



blood1$subject <- "GY10"
blood2$subject <- "GY4"
blood3$subject <- "GY6"
blood4$subject <- "GY9"

bld <- merge(blood1, y = c(blood2, blood3, blood4))
bld <- blood1


VlnPlot(bld, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "subject", pt.size = 0)

VlnPlot(rawbld, features = c("Cd34"), pt.size = 0)

FeatureScatter(bld, feature1 = "Hba-a1", feature2 = "nFeature_RNA", raster=T) + 
FeatureScatter(bld, feature1 = "Hba-a1", feature2 = "nCount_RNA", raster=T) +
FeatureScatter(bld, feature1 = "Erythroidness", feature2 = "RNA.diversity") + 
FeatureScatter(bld, feature1 = "Erythroidness", feature2 = "percent.mt")

bldsave <- #bld




#==========================================================================
bld[["RNA.diversity"]] <- bld$nFeature_RNA / bld$nCount_RNA
		#Values: 0.8-0.95
bld[["RNA.homogeneity"]] <- bld$nCount_RNA / bld$nFeature_RNA

bld[["percent.mt"]] <- PercentageFeatureSet(bld, pattern = "^mt-")
		#Mito genes
bld[["percent.cd"]] <- PercentageFeatureSet(bld, pattern = "^Cd[ck]")
		#Cell cycle genes
bld[["percent.rb"]] <- PercentageFeatureSet(bld, pattern =  "^Rb[sl]")
		#Ribo genes


# Counts for alpha/beta Hb poduction
bld$HBax <- PercentageFeatureSet(bld, pattern = "^Hba-") * bld$nCount_RNA / 100
bld$HBbx <- PercentageFeatureSet(bld, pattern = "^Hbb-") * bld$nCount_RNA / 100

# Not v helpful
#bld$BAratio <- ifelse(bld$HBbx > 0 & bld$HBax > 0, (bld$HBbx / bld$HBax), 0)
#bld$BAspread <- ifelse(bld$HBbx > 0 & bld$HBax > 0, (bld$HBbx - bld$HBax), 0)



bld[["Erythroidness"]] <- PercentageFeatureSet(bld, features = c(
		"Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt",
		"Hemgn", "Ermap", "Slc25a21", "Slc4a1"))
		#Values: 0.2-40

bld[["B.cellness"]] <- PercentageFeatureSet(bld, features = c(
		"Pax5", "Ebf4", "Bank1", "Bach2",
		"Vpreb1", "Lef1", "Cecr2", "Vpreb3"))
		#Values: 0.5-1.5

bld[["T.cellness"]] <- PercentageFeatureSet(bld, features = c(
		"Ccl5", "Nkg7", "Bcl11b", "Skap1",
		"Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
		#Values: 1-2

bld[["Macrophageness"]] <- PercentageFeatureSet(bld, features = c(
		"Ccl5", "Nkg7", "Bcl11b", "Skap1",
		"Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
		#Values:

bld[["Neutrophilness"]] <- PercentageFeatureSet(bld, features = c(
		"Ltf", "Ngp", "Mmp8", "Retnlg",
		"S100a8", "Ltf", "S100a9", "Elane"))
		#Values: 0-4

bld[["Monocyteness"]] <- PercentageFeatureSet(bld, features = c(
		"F13a1", "Gria3", "Prtn3", "Mpo", 
		"Crip1", "Lgals1", "Tmsb10", "Ly6a"))
		#Values: 0.5-2
#===========================================================================




################ Standard processing ###############


bld <- FindVariableFeatures(bld, selection.method = "vst", nfeatures = 6000)
all.genes <- rownames(bld)
bld <- ScaleData(bld, features = VariableFeatures(object = bld))
bld <- RunPCA(bld, features = VariableFeatures(object = bld))
ElbowPlot(bld, 50)

bld <- FindNeighbors(bld, dims = 1:35, k.param=120)
bld <- FindClusters(bld, resolution = 0.1)
bld <- RunUMAP(bld, dims = 1:35, return.model = F, n.neighbors=120)
dim <- DimPlot(bld, reduction = "umap",label = TRUE, pt.size = .75, shuffle = T, raster = T) + NoLegend()


feats <- FeaturePlot(bld, features = c(
		"Erythroidness", "Hemgn", 
		"Bpgm", "Ighm"), 
		min.cutoff = c(2,0,0,0), max.cutoff = c(80,1,100,100), 
		ncol = 2, pt.size = 2, raster=T)
dim + feats









################# Naive marker-based cell selection #################


# Common Myeloid Progenitor:

CMPdata <- FetchData(bld, vars = c('Kit', "Gata1", "Slamf1", "Hmgb3", "Cd34", "Flt3", "Ly6a"))

CMP_scatter <- ggplot(data = CMPdata) +
		geom_point(mapping = aes(x = Cd34, y = Kit, color = Ly6a, size = Slamf1, alpha = .5), 
		position = 'jitter') + scale_color_gradientn(colors = c("blue", "red"))

CMPs_highlighted <- FeaturePlot(bld, features = c(
		"Cd34", "Kit", 
		"CMPs", "Ly6a"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(1,.5,10,.5), 
		ncol = 2, pt.size = 2, raster=T)


# Erythroid Progenitors:

EPs_data <- FetchData(bld, vars = c('Kit',"Gata1", "Hemgn", "Epor", "Sox6", "Tal1", "Gypc","EPs","Ermap"))

EPs_scatter <- ggplot(data = blastdata) +
		geom_point(mapping = aes(x = Tal1 + Gypc + Sox6, y = Gata1 + Kit, color = EPs, size = Epor + Hemgn + Ermap, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))

EPs_highlighted <- FeaturePlot(bld, features = c(
		"EPs", "Hemgn", 
		"Tall", "Sox6"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,.5,.5,.5), 
		ncol = 2, pt.size = 2, raster=T)


# Proerythroblasts:

PEBs_data <- FetchData(bld, vars = c("EBs","Slc4a1","Epb42","Tnfrsf1a",'Tfrc', 'Bpgm', "Hemgn", "Cd44", "Epor", "Ermap","Kit", "PEBs", "Klf1","Tal1", "Gata1"))

PEBs_scatter <- ggplot(data = PEBs_data) +
		geom_point(mapping = aes(x = Tfrc, y = Cd44, color = Tnfrsf1a, size = PEBs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
PEBs_scatter

PEBs_highlighted <- FeaturePlot(bld, features = c(
		"EBs", "Epb42", 
		"Tnfrsf1a", "Cd44"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,1.6,1,100), 
		ncol = 2, pt.size = 2, raster=T)
dim + PEBs_highlighted


# Erythroblasts (Baso/Poly/Ortho):

EBs_data <- FetchData(bld, vars = c("EBs","Slc4a1","Epb42","Tnfrsf1a",'Tfrc', 'Bpgm', "Hemgn", "Cd44", "Epor", "Ermap","Kit", "PEBs", "Klf1","Tal1", "Gata1"))

EBs_scatter <- ggplot(data = EBs_data) +
		geom_point(mapping = aes(x = Tfrc, y = Cd44, color = Tnfrsf1a, size = EBs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
EBs_scatter


EBs_highlighted <- FeaturePlot(bld, features = c(
		"EBs", "Epb42", 
		"Tnfrsf1a", "Cd44"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,1.6,1,100), 
		ncol = 2, pt.size = 2, raster=T)
dim + EBs_highlighted


# Reticulocytes:

RCs_data <- FetchData(bld, vars = c('Tfrc', 'Kit', 'Bpgm', "Tnfrsf1a", "Gypa", "Bcl11a", "Slc11a2", "Ldb1","Sox6", "Myb", "Cd24a",  "HBbx", "HBax","RCs","Bnip3l"))


RCs_scatter <- ggplot(data = RCs_data) +
		geom_point(mapping = aes(x = Tfrc, y = log(HBbx+HBax), color = Bnip3l, size = RCs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
RCs_scatter


RCs_highlighted <- FeaturePlot(bld, features = c(
		"RCs", "Elk4", 
		"Bnip3l", "Ejection"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,.5,.5,.5), 
		ncol = 2, pt.size = 2, raster=T)
dim+RCs_highlighted


# Erythrocytes:

ECs_data <- FetchData(bld, vars = c('Tfrc', 'Kit', 'Bpgm', "Gata1", "Hemgn", "Tnfrsf1a", "Slamf1",
		"Eng", "Gypa", "Hmgb3", "Bcl11a", "Cd36", "Ldb1","Sox6", "Myb", "Cd24a", "Klf1","Tal1","Hebp1","Hmbs", "Stat1", "HBbx", "EPs","HBax"))


ECs_scatter <- ggplot(data = ECs_data) +
		geom_point(mapping = aes(x = HBax, y = Tfrc, color = Cd24a, size = Bpgm, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
ECs_scatter


ECs_highlighted <- FeaturePlot(bld, features = c(
		"Bpgm", "HBax", 
		"HBbx", "Hemgn"), 
		min.cutoff = c(0,20,20,0), max.cutoff = c(.5,1000,1000,.5), 
		ncol = 2, pt.size = 2, raster=T)
ECs_highlighted






########## Subset erythroid cells (both methods) ###############

clustmarks <- FindMarkers(bld, ident.1 = 20, min.pct = 0.25, only.pos = T, logfc.threshold = 3)
head(clustmarks, n = 40)


# Identitify clusters of interest and subset
ery_clusters <- subset(bld, idents = c(30,33,18))

# Group handpicked cells and subset
bld$handpicked_celltype <- ifelse(bld$RCs > 0.5, "Reticulocyte",
		ifelse(bld$EBs > 0.5, "Erythroblast",
		ifelse(bld$PEBs > 0.5, "Proerythroblasts",
		ifelse(bld$EPs > 0.5, "Erythroid Progenitor",
		ifelse(bld$CMPs > 0.5, "Common Myeloid Progenitor",
		"0")))))



handpicked_cells <- subset(bld, subset = handpicked_celltype == "0", invert = T)
table(handpicked_cells$handpicked_celltype)
handpicked_cells

handpicked_cells$selection_method <- "handpicked"
ery_clusters$selection_method <- "from cluster"

all_erys <- merge(handpicked_cells, ery_clusters)
all_erys




########## Reprocess Isolated erythroid cells ############

all_erys <- NormalizeData(all_erys)
all_erys <- FindVariableFeatures(all_erys, selection.method = "vst", nfeatures = 10000)
all_erys.genes <- rownames(all_erys)
all_erys <- ScaleData(all_erys, features = all_erys.genes)
all_erys <- RunPCA(all_erys, features = all_erys.genes)
ElbowPlot(all_erys, 50)


all_erys <- FindNeighbors(all_erys, dims = 1:5)
all_erys <- FindClusters(all_erys, resolution = 0.5)
all_erys <- RunUMAP(all_erys, dims = 1:5, reduction = "pca", return.model = TRUE)

all_erysdim <- DimPlot(all_erys,raster=F,group.by="ident",reduction = "umap", repel = T, label = T, shuffle = T, pt.size = 2) + NoLegend()
all_erysdim

all_erysfeats <- FeaturePlot(all_erys, features = c(
		"Cd34", "Cd44", "Hemgn", "Epor",
		"HBax", "Ejection", "Bpgm", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
		min.cutoff = c(0,0,0,0,
				0,0,0,.05),
		max.cutoff = c(1,1,1,1,
				1000,1,3,.9))
all_erysfeats




########## Compare handpicked to computed celltypes and clean up dataset ##############

saveRDS(all_erys, file = "Combined_Method_FEMUR_cells.rds")


table(all_erys$handpicked_celltype,Idents(all_erys))


clustmarks <- FindMarkers(all_erys, ident.1 = 1, min.pct = 0.25, only.pos = T, logfc.threshold = 2)
head(clustmarks, n = 40)


all_erys.markers <- FindAllMarkers(all_erys, min.pct = 0.25, logfc.threshold = 1)
topmarks <- all_erys.markers %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC)

heat = DoHeatmap(all_erys, features = topmarks$gene, size = 4, angle = 10) + NoLegend()
all_erysdim | heat

all_erys <- subset(all_erys, idents = c(8,9,10), invert =T)











###
