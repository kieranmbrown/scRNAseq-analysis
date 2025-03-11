library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

humerus1.data <- ReadMtx(
    mtx = "all_data/Humerus/HYG1/matrix.mtx",
    features = "all_data/Humerus/HYG1/features.tsv",
    cells = "all_data/Humerus/HYG1/barcodes.tsv")

hmrs1 <- CreateSeuratObject(
    counts = humerus1.data)

humerus2.data <- ReadMtx(
    mtx = "all_data/Humerus/HYG4/matrix.mtx",
    features = "all_data/Humerus/HYG4/features.tsv",
    cells = "all_data/Humerus/HYG4/barcodes.tsv")

hmrs2 <- CreateSeuratObject(
    counts = humerus2.data)

humerus3.data <- ReadMtx(
    mtx = "all_data/Humerus/HYG6/matrix.mtx",
    features = "all_data/Humerus/HYG6/features.tsv",
    cells = "all_data/Humerus/HYG6/barcodes.tsv")

hmrs3 <- CreateSeuratObject(
    counts = humerus3.data)

humerus4.data <- ReadMtx(
    mtx = "all_data/Humerus/HYG9/matrix.mtx",
    features = "all_data/Humerus/HYG9/features.tsv",
    cells = "all_data/Humerus/HYG9/barcodes.tsv")

hmrs4 <- CreateSeuratObject(
    counts = humerus4.data)

humerus1$subject <- "HYG1"
humerus2$subject <- "HYG4"
humerus3$subject <- "HYG6"
humerus4$subject <- "HYG9"

hmrs <- merge(humerus1, y = c(humerus2, humerus3, humerus4))
hmrs




#==========================================================================
hmrs[["RNA.diversity"]] <- hmrs$nFeature_RNA / hmrs$nCount_RNA
		#Values: 0.8-0.95
hmrs[["RNA.homogeneity"]] <- hmrs$nCount_RNA / hmrs$nFeature_RNA

hmrs[["percent.mt"]] <- PercentageFeatureSet(hmrs, pattern = "^mt-")
		#Mito genes
hmrs[["percent.cd"]] <- PercentageFeatureSet(hmrs, pattern = "^Cd[ck]")
		#Cell cycle genes
hmrs[["percent.rb"]] <- PercentageFeatureSet(hmrs, pattern =  "^Rb[sl]")
		#Ribo genes


VlnPlot(hmrs, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0)

FeatureScatter(hmrs, feature1 = "percent.mt", feature2 = "nFeature_RNA",raster=T) + 
FeatureScatter(hmrs, feature1 = "percent.mt", feature2 = "nCount_RNA",raster=T)

hmrs <- subset(hmrs, subset = nFeature_RNA > 30 & nFeature_RNA < 2000 & 
		percent.mt < 10 & nCount_RNA < 10000)

hmrs <- subset(hmrs, subset = nFeature_RNA > 1200 & subject=="HYG4", invert=T)




# Counts for alpha/beta Hb poduction
hmrs$HBax <- PercentageFeatureSet(hmrs, pattern = "^Hba-") * hmrs$nCount_RNA / 100
hmrs$HBbx <- PercentageFeatureSet(hmrs, pattern = "^Hbb-") * hmrs$nCount_RNA / 100

# Not v helpful
hmrs$BAratio <- ifelse(hmrs$HBbx > 0 & hmrs$HBax > 0, (hmrs$HBbx / hmrs$HBax), 0)
hmrs$BAspread <- ifelse(hmrs$HBbx > 0 & hmrs$HBax > 0, (hmrs$HBbx - hmrs$HBax), 0)



hmrs[["Erythroidness"]] <- PercentageFeatureSet(hmrs, features = c(
		"Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt",
		"Hemgn", "Ermap", "Slc25a21", "Slc4a1"))
		#Values: 0.2-40

hmrs[["B.cellness"]] <- PercentageFeatureSet(hmrs, features = c(
		"Pax5", "Ebf4", "Bank1", "Bach2",
		"Vpreb1", "Lef1", "Cecr2", "Vpreb3"))
		#Values: 0.5-1.5

hmrs[["T.cellness"]] <- PercentageFeatureSet(hmrs, features = c(
		"Ccl5", "Nkg7", "Bcl11b", "Skap1",
		"Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
		#Values: 1-2

hmrs[["Macrophageness"]] <- PercentageFeatureSet(hmrs, features = c(
		"Ccl5", "Nkg7", "Bcl11b", "Skap1",
		"Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
		#Values:

hmrs[["Neutrophilness"]] <- PercentageFeatureSet(hmrs, features = c(
		"Ltf", "Ngp", "Mmp8", "Retnlg",
		"S100a8", "Ltf", "S100a9", "Elane"))
		#Values: 0-4

hmrs[["Monocyteness"]] <- PercentageFeatureSet(hmrs, features = c(
		"F13a1", "Gria3", "Prtn3", "Mpo", 
		"Crip1", "Lgals1", "Tmsb10", "Ly6a"))
		#Values: 0.5-2
#===========================================================================



############## Specific celltype filters for erythroid lineage ###################

hmrs[["CMPs"]] <- ifelse((hmrs$RNA@data["Ly6a",] < 0.5 & hmrs$RNA@data["Kit",] > 0.5), 
		hmrs$RNA@data["Cd34",], 0)

hmrs[["EPs"]] <- ifelse((
		(hmrs$RNA@data["Tal1",] + hmrs$RNA@data["Gypc",] + hmrs$RNA@data["Sox6",] > 0.5) &
		(hmrs$RNA@data["Epor",] + hmrs$RNA@data["Hemgn",]  + hmrs$RNA@data["Ermap",] + hmrs$RNA@data["Eng",] > 0.5) &
		(hmrs$RNA@data["Gata1",]  + hmrs$RNA@data["Kit",] > 0.5)),
		(hmrs$RNA@data["Epor",] + hmrs$RNA@data["Hemgn",]  + hmrs$RNA@data["Ermap",] +
		hmrs$RNA@data["Tal1",] + hmrs$RNA@data["Gypc",] + hmrs$RNA@data["Sox6",]), 0)

hmrs[["PEBs"]] <- ifelse((
		((hmrs$RNA@data["Epor",] + hmrs$RNA@data["Ermap",] + hmrs$RNA@data["Hemgn",] + hmrs$RNA@data["Gypa",] > 0.5) | 
		(hmrs$RNA@data["Tal1",] > .5) | 
		(hmrs$RNA@data["Tfrc",] > 3.5)) &
		(hmrs$RNA@data["Cd44",] > 0.5) & 
		(hmrs$RNA@data["Tfrc",] > 0.5) & 
		(hmrs$RNA@data["Epb42",] + hmrs$RNA@data["Bpgm",] < 1)),

		(hmrs$RNA@data["Epor",] + hmrs$RNA@data["Hemgn",]  + 
		hmrs$RNA@data["Tfrc",] + hmrs$RNA@data["Tal1",]), 0)

hmrs[["EBs"]] <- ifelse((
		((hmrs$RNA@data["Epor",] + hmrs$RNA@data["Ermap",] + hmrs$RNA@data["Hemgn",] > 0.5) | 
		(hmrs$RNA@data["Tnfrsf1a",] < 1)) &
		(hmrs$RNA@data["Bpgm",] + hmrs$RNA@data["Slc4a1",] + hmrs$RNA@data["Epb42",] + hmrs$RNA@data["Ppox",] > 0.5) &
		(hmrs$RNA@data["Cd44",] > 0.5) & (hmrs$RNA@data["Tfrc",] > 1)),

		(hmrs$RNA@data["Epb42",] + hmrs$RNA@data["Hemgn",] + 
		hmrs$RNA@data["Bpgm",] + hmrs$RNA@data["Tfrc",]), 0)

hmrs[["RCs"]] <- ifelse(
		(log(hmrs$HBbx + hmrs$HBax)/7 + hmrs$RNA@data["Tfrc",]/4 > 1) & 
		(hmrs$RNA@data["Tfrc",] > 0.5),

		(hmrs$RNA@data["Epb42",] +
		hmrs$RNA@data["Bpgm",] + hmrs$RNA@data["Tfrc",]  + 
		hmrs$RNA@data["Hba-a1",] + hmrs$RNA@data["Ppox",]), 0)

hmrs[["Ejection"]] <- ifelse(hmrs$HBbx + hmrs$HBax > 25,
		(hmrs$RNA@data["Bnip3l",]),0)





################ Standard processing ###############
for (set in list(hmrs1, hmrs2, hmrs3, hmrs4)){
	print(set)

set <- NormalizeData(set)
set <- FindVariableFeatures(set, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(set)
set <- ScaleData(set, features = all.genes)
set <- RunPCA(set, features = VariableFeatures(object = set))

hmrs <- FindNeighbors(hmrs, dims = 1:30)
hmrs <- FindClusters(hmrs, resolution = 1.5)
hmrs <- RunUMAP(hmrs, dims = 1:30, return.model = F)
dim <- DimPlot(hmrs, reduction = "umap",label = TRUE, pt.size = .75, shuffle = T, raster = T) + NoLegend()

feats <- FeaturePlot(hmrs, features = c(
		"Erythroidness", "Hba-a1", 
		"Hemgn", "Ermap"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(10,1,20,20), 
		ncol = 2, pt.size = 2, raster=T)
dim + feats

plot_name <- paste0("UMAP_", hmrs1,".png")
ggsave(plot_name, marker_feats)

}


datasets <- list(hmrs1 = hmrs1, hmrs2 = hmrs2, hmrs3 = hmrs3, hmrs4 = hmrs4)

for (name in names(datasets)) {
  
  set <- datasets[[name]]
  print(name)

  set <- NormalizeData(set)
  set <- FindVariableFeatures(set, selection.method = "vst", nfeatures = 2500)
  all.genes <- rownames(set)
  set <- ScaleData(set, features = all.genes)
  set <- RunPCA(set, features = VariableFeatures(object = set))

  set <- FindNeighbors(set, dims = 1:30)
  set <- FindClusters(set, resolution = 1.5)
  set <- RunUMAP(set, dims = 1:30, return.model = F)
  dim <- DimPlot(set, reduction = "umap",label = TRUE, pt.size = .75, shuffle = T, raster = T) + NoLegend()

  feats <- FeaturePlot(set, features = c(
    	"Gypa", "Hba-a1", 
    	"Hemgn", "Ermap"), 
    	min.cutoff = c(0,0,0,0), max.cutoff = c(1,1,20,20), 
    	ncol = 2, pt.size = 2, raster=T)
	output <- dim + feats

  plot_name <- paste0("UMAP_", name, ".png")
  ggsave(plot_name, output)
}



hmrs <- NormalizeData(hmrs)
hmrs <- FindVariableFeatures(hmrs, selection.method = "vst", nfeatures = 10000)
all.genes <- rownames(hmrs)
hmrs <- ScaleData(hmrs, features = all.genes)
hmrs <- RunPCA(hmrs, features = VariableFeatures(object = hmrs))


hmrs <- FindNeighbors(hmrs, dims = 1:45)
hmrs <- FindClusters(hmrs, resolution = 1.5)
hmrs <- RunUMAP(hmrs, dims = 1:15, return.model = F)
dim <- DimPlot(hmrs, reduction = "umap",label = TRUE, pt.size = .75, shuffle = T, raster = T) + NoLegend()


feats <- FeaturePlot(hmrs, features = c(
		"Erythroidness", "Hba-a1", 
		"Hemgn", "Ermap"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(10,1,20,20), 
		ncol = 2, pt.size = 2, raster=T)
dim + feats









################# Naive marker-based cell selection #################


# Common Myeloid Progenitor:

CMPdata <- FetchData(hmrs, vars = c('Kit', "Gata1", "Slamf1", "Hmgb3", "Cd34", "Flt3", "Ly6a","CMPs"))

CMP_scatter <- ggplot(data = CMPdata) +
		geom_point(mapping = aes(x = Cd34, y = Kit, color = CMPs, size = Ly6a, alpha = .5), 
		position = 'jitter') + scale_color_gradientn(colors = c("blue", "red"))
CMP_scatter

CMPs_highlighted <- FeaturePlot(hmrs, features = c(
		"Cd34", "Kit", 
		"CMPs", "Ly6a"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(1,.5,1,.5), 
		ncol = 2, pt.size = 2, raster=T)
CMPs_highlighted

# Erythroid Progenitors:

EPs_data <- FetchData(hmrs, vars = c('Kit',"Gata1", "Hemgn", "Epor", "Sox6", "Tal1", "Gypc","EPs","Ermap"))

EPs_scatter <- ggplot(data = EPs_data) +
		geom_point(mapping = aes(x = Tal1 + Gypc + Sox6, y = Gata1 + Kit, color = EPs, size = Epor + Hemgn + Ermap, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
EPs_scatter

EPs_highlighted <- FeaturePlot(hmrs, features = c(
		"EPs", "Hemgn", 
		"Tal1", "Sox6"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,.5,.5,.5), 
		ncol = 2, pt.size = 2, raster=T)
EPs_highlighted

# Proerythroblasts:

PEBs_data <- FetchData(hmrs, vars = c("EBs","Slc4a1","Epb42","Tnfrsf1a",'Tfrc', 'Bpgm', "Hemgn", "Cd44", "Epor", "Ermap","Kit", "PEBs", "Klf1","Tal1", "Gata1"))

PEBs_scatter <- ggplot(data = PEBs_data) +
		geom_point(mapping = aes(x = Tfrc, y = Cd44, color = Tnfrsf1a, size = PEBs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
PEBs_scatter

PEBs_highlighted <- FeaturePlot(hmrs, features = c(
		"PEBs", "Tal1", 
		"Gypa", "Cd44"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,1.6,1,100), 
		ncol = 2, pt.size = 2, raster=T)
dim + PEBs_highlighted


# Erythroblasts (Baso/Poly/Ortho):

EBs_data <- FetchData(hmrs, vars = c("EBs","Slc4a1","Epb42","Tnfrsf1a",'Tfrc', 'Bpgm', "Hemgn", "Cd44", "Epor", "Ermap","Kit", "PEBs", "Klf1","Tal1", "Gata1"))

EBs_scatter <- ggplot(data = EBs_data) +
		geom_point(mapping = aes(x = Tfrc, y = Cd44, color = Tnfrsf1a, size = EBs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
EBs_scatter


EBs_highlighted <- FeaturePlot(hmrs, features = c(
		"EBs", "Epb42", 
		"Tnfrsf1a", "Cd44"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,1.6,1,100), 
		ncol = 2, pt.size = 2, raster=T)
dim + EBs_highlighted


# Reticulocytes:

RCs_data <- FetchData(hmrs, vars = c('Tfrc', 'Kit', 'Bpgm', "Tnfrsf1a", "Gypa", "Bcl11a", "Slc11a2", "Ldb1","Sox6", "Myb", "Cd24a",  "HBbx", "HBax","RCs","Bnip3l"))


RCs_scatter <- ggplot(data = RCs_data) +
		geom_point(mapping = aes(x = Tfrc, y = log(HBbx+HBax), color = Bnip3l, size = RCs, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
RCs_scatter


RCs_highlighted <- FeaturePlot(hmrs, features = c(
		"RCs", "Elk4", 
		"percent.mt", "percent.rb"), 
		min.cutoff = c(0,0,3,0), max.cutoff = c(.5,.5,6,.35), 
		ncol = 2, pt.size = 2, raster=T)
dim+RCs_highlighted


# Erythrocytes:

ECs_data <- FetchData(hmrs, vars = c('Tfrc', 'Kit', 'Bpgm', "Gata1", "Hemgn", "Tnfrsf1a", "Slamf1",
		"Eng", "Gypa", "Hmgb3", "Bcl11a", "Cd36", "Ldb1","Sox6", "Myb", "Cd24a", "Klf1","Tal1","Hebp1","Hmbs", "Stat1", "HBbx", "EPs","HBax"))


ECs_scatter <- ggplot(data = ECs_data) +
		geom_point(mapping = aes(x = HBax, y = Tfrc, color = Cd24a, size = Bpgm, alpha = .5), position = 'jitter') +
		scale_color_gradientn(colors = c("blue", "red"))
ECs_scatter


ECs_highlighted <- FeaturePlot(hmrs, features = c(
		"Alas2", "Klf1", 
		"Gypa", "Hemgn"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(.5,1,1,.5), 
		ncol = 2, pt.size = 2, raster=T)
ECs_highlighted






########## Subset erythroid cells (both methods) ###############

clustmarks <- FindMarkers(hmrs, ident.1 = 4, min.pct = 0.25, only.pos = T, logfc.threshold = 1)
head(clustmarks, n = 40)


# Identitify clusters of interest and subset
ery_clusters <- subset(hmrs, idents = c(25,26,33))

# Group handpicked cells and subset
hmrs$handpicked_celltype <- ifelse(hmrs$RCs > 0.5, "Reticulocyte",
		ifelse(hmrs$EBs > 0.5, "Erythroblast",
		ifelse(hmrs$PEBs > 0.5, "Proerythroblasts",
		ifelse(hmrs$EPs > 0.5, "Erythroid Progenitor",
		ifelse(hmrs$CMPs > 0.5, "Common Myeloid Progenitor",
		"0")))))



handpicked_cells <- subset(hmrs, subset = handpicked_celltype == "0", invert = T)
table(handpicked_cells$handpicked_celltype)
handpicked_cells

cells1 <- colnames(handpicked_cells)
cells2 <- colnames(ery_clusters)

shared_cells <- intersect(x = cells1, y = cells2)

ery_clusters <- subset(ery_clusters, cells = shared_cells, invert = T)


handpicked_cells$selection_method <- "handpicked"
ery_clusters$selection_method <- "from cluster"



all_erys <- merge(handpicked_cells, ery_clusters)
table(all_erys$selection_method)
all_erys


########## Reprocess Isolated erythroid cells ############

all_erys <- NormalizeData(all_erys)
all_erys <- FindVariableFeatures(all_erys, selection.method = "vst", nfeatures = 10000)
all_erys.genes <- rownames(all_erys)
all_erys <- ScaleData(all_erys, features = all_erys.genes)
all_erys <- RunPCA(all_erys, features = all_erys.genes)
ElbowPlot(all_erys, 50)


all_erys <- FindNeighbors(all_erys, dims = 1:6)
all_erys <- FindClusters(all_erys, resolution = 1.5)
all_erys <- RunUMAP(all_erys, dims = 1:6, reduction = "pca", return.model = TRUE)

all_erysdim <- DimPlot(all_erys,raster=F,group.by="subject",reduction="umap",repel=T,label=T, shuffle = T, pt.size = 2) + NoLegend()
quick_feats <- FeaturePlot(all_erys, features = c(
		"Cd34", "Cd44", "Hemgn", "HBax"),	 pt.size = 4,raster = T,ncol=2, 
		min.cutoff = c(0,0,0,0),max.cutoff = c(1,1,1,500))
all_erysdim|quick_feats

all_erysfeats <- FeaturePlot(all_erys, features = c(
		"Cd34", "Cd44", "Hemgn", "Tfrc",
		"HBax", "Bpgm", "Ejection", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
		min.cutoff = c(0,0,0,0,
				0,0,0,.05),
		max.cutoff = c(1,1,1,1,
				1000,3,3,.9))
all_erysfeats




########## Compare handpicked to computed celltypes and clean up dataset ##############


saveRDS(all_erys, file = "Combined_Method_FEMUR_cells.rds")



table(all_erys$handpicked_celltype,Idents(all_erys))


clustmarks <- FindMarkers(all_erys, ident.1 = 5, min.pct = 0.25, only.pos = T, logfc.threshold = 1)
head(clustmarks, n = 40)


all_erys.markers <- FindAllMarkers(all_erys, min.pct = 0.25, logfc.threshold = 1)
topmarks <- all_erys.markers %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC)

heat = DoHeatmap(all_erys, features = topmarks$gene, size = 4, angle = 10) + NoLegend()
all_erysdim | heat




VlnPlot(all_erys, features = c("RNA.diversity", "nFeature_RNA","nCount_RNA","HBax"),group.by="subject", ncol=4))

VlnPlot(all_erys, features = c("Hba-a1", "Gypa","Alas2"), ncol=3)/VlnPlot(testing, features = c("Hba-a1", "Gypa","Alas2"), ncol=3)

VlnPlot(all_erys, features = c("Hba-a1", "Gypa","Alas2"), ncol=3)/VlnPlot(testing2, features = c("Hba-a1", "Gypa","Alas2"), ncol=3)


all_erys <- subset(all_erys, idents = c(0,1,4,5,7), invert =T)

all_erys <- subset(all_erys, idents = c(5,6), subset = Gypa+Alas2 < 0.5, invert =T)

testing2 <- subset(all_erys, idents = c(1,3,6,7), subset = Cd44 > 0.5, invert =F)




all_erysfeats <- FeaturePlot(testing, features = c(
		"Cd34", "Cd44", "Hemgn", "Tfrc",
		"HBax", "Slamf1", "Eng", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
		min.cutoff = c(0,0,0,0,
				0,0,0,.05),
		max.cutoff = c(1,1,1,1,
				1000,3,3,.9))
all_erysfeats




###
