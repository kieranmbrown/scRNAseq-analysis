
####################################
# Processing script for Femur data #
####################################


# Read in files and inspect cell counts

fmr1.data <- ReadMtx(
  mtx = "all_data/RR10/FLT_data/2cWTFLT/matrix.mtx",
  features = "all_data/RR10/FLT_data/2cWTFLT/features.tsv",
  cells = "all_data/RR10/FLT_data/2cWTFLT/barcodes.tsv")

fmr1 <- CreateSeuratObject(
  counts = fmr1.data)
fmr1

fmr2.data <- ReadMtx(
  mtx = "all_data/RR10/FLT_data/2eWTFLT/matrix.mtx",
  features = "all_data/RR10/FLT_data/2eWTFLT/features.tsv",
  cells = "all_data/RR10/FLT_data/2eWTFLT/barcodes.tsv")

fmr2 <- CreateSeuratObject(
  counts = fmr2.data)
fmr2

fmr3.data <- ReadMtx(
  mtx = "all_data/RR10/FLT_data/6dWTFLT/matrix.mtx",
  features = "all_data/RR10/FLT_data/6dWTFLT/features.tsv",
  cells = "all_data/RR10/FLT_data/6dWTFLT/barcodes.tsv")

fmr3 <- CreateSeuratObject(
  counts = fmr3.data)
fmr3

# Assign metadata fields, merge to cohort-level object

fmr1$subject <- "SO-01"
fmr2$subject <- "SO-02"
fmr3$subject <- "SO-03"

fmr <- merge(fmr1, c(fmr2,fmr3))

fmr$Cohort <- "SO"
fmr$age <- "Old"
fmr$experiment <- "Space"





fmr1.data <- ReadMtx(
  mtx = "all_data/RR10/GC_data/4cWTGC/matrix.mtx",
  features = "all_data/RR10/GC_data/4cWTGC/features.tsv",
  cells = "all_data/RR10/GC_data/4cWTGC/barcodes.tsv")

fmr1 <- CreateSeuratObject(
  counts = fmr1.data)
fmr1

fmr2.data <- ReadMtx(
  mtx = "all_data/RR10/GC_data/4dWTGC/matrix.mtx",
  features = "all_data/RR10/GC_data/4dWTGC/features.tsv",
  cells = "all_data/RR10/GC_data/4dWTGC/barcodes.tsv")

fmr2 <- CreateSeuratObject(
  counts = fmr2.data)
fmr2

fmr3.data <- ReadMtx(
  mtx = "all_data/RR10/GC_data/8dWTGC/matrix.mtx",
  features = "all_data/RR10/GC_data/8dWTGC/features.tsv",
  cells = "all_data/RR10/GC_data/8dWTGC/barcodes.tsv")

fmr3 <- CreateSeuratObject(
  counts = fmr3.data)
fmr3

# Assign metadata fields, merge to cohort-level object

fmr1$subject <- "RR10GC-01"
fmr2$subject <- "RR10GC-02"
fmr3$subject <- "RR10GC-03"

fmrgc <- merge(fmr1, c(fmr2,fmr3))

fmrgc$Cohort <- "RR10GC"
fmrgc$age <- "Old"
fmrgc$experiment <- "Ground"


fmrall <- merge(fmr,fmrgc)

fmr <- fmrall

# Filter cells for standardization and inspect for problems

fmr[["RNA.diversity"]] <- fmr$nFeature_RNA / fmr$nCount_RNA
fmr[["percent.mt"]] <- PercentageFeatureSet(fmr, pattern = "^mt-")

fmr <- subset(fmr, subset = percent.mt<5 & nFeature_RNA<6000 & nCount_RNA<30000 & nCount_RNA>30)

VlnPlot(fmr, features= c("nFeature_RNA", "nCount_RNA", "RNA.diversity", "percent.mt"),
        ncol=4, pt.size=0,group.by="subject")
table(fmr$subject)


# Specifically check 'Sample Order Bias' which reduces # of low-RNA cells 
# captured for the samples listed first on the sample sheet for the sequencing run 

FeatureScatter(fmr, feature1="nFeature_RNA", feature2="nCount_RNA", split.by="subject")

VlnPlot(fmr, features="nFeature_RNA", pt.size=0,group.by="subject")

fmr$sample_number <- as.numeric(str_extract(fmr$subject, "\\d+$"))

fmr$sample_rank <- fmr$sample_number
#library(DescTools)
for(cohort in unique(fmr$Cohort)){
  fmr$sample_rank <- ifelse(fmr$Cohort == cohort, Rank(subset(fmr,subset=Cohort==cohort)$sample_number, ties.method="dense"), fmr$sample_rank)}
  

for(cohort in unique(fmr$Cohort)){
  subset_cohort <- subset(fmr, subset=Cohort==cohort)
  
  for(sample in unique(subset_cohort$subject))
    subject_rank <- 
  fmr$sample_rank <- ifelse(fmr$Cohort == cohort, rank(fmr$sample_number, ties.method="dense"), fmr$sample_rank)}

VlnPlot(fmr, features="nFeature_RNA", pt.size=0,group.by="subject", sort="sample_rank")



# Calculate some useful fields for annotation
#==========================================================================
fmr[["RNA.homogeneity"]] <- fmr$nCount_RNA / fmr$nFeature_RNA
fmr[["percent.cd"]] <- PercentageFeatureSet(fmr, pattern = "^Cd[ck]")
#Cell cycle genes, CDCs and CDKs
fmr[["percent.rb"]] <- PercentageFeatureSet(fmr, pattern =  "^Rb[sl]")*fmr$nCount_RNA
#Ribo small/large segment genes

# Counts for alpha/beta Hb poduction
fmr$HBax <- PercentageFeatureSet(fmr, pattern = "^Hba-") * fmr$nCount_RNA / 100
fmr$HBbx <- PercentageFeatureSet(fmr, pattern = "^Hbb-") * fmr$nCount_RNA / 100

# Not v helpful
fmr$BAratio <- ifelse(fmr$HBbx > 0 & fmr$HBax > 0, (fmr$HBbx / fmr$HBax), 0)
fmr$BAspread <- ifelse(fmr$HBbx > 0 & fmr$HBax > 0, (fmr$HBbx - fmr$HBax), 0)

fmr[["Erythroidness"]] <- PercentageFeatureSet(fmr, features = c(
  "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt",
  "Hemgn", "Ermap", "Slc25a21", "Slc4a1"))
#Values: 0.2-40

fmr[["B.cellness"]] <- PercentageFeatureSet(fmr, features = c(
  "Pax5", "Ebf4", "Bank1", "Bach2",
  "Vpreb1", "Lef1", "Cecr2", "Vpreb3"))
#Values: 0.5-1.5

fmr[["T.cellness"]] <- PercentageFeatureSet(fmr, features = c(
  "Ccl5", "Nkg7", "Bcl11b", "Skap1",
  "Ms4a4b", "Tox", "Cd3e", "Ikzf2"))
#Values: 1-2

fmr[["Macrophageness"]] <- PercentageFeatureSet(fmr, features = c(
  "Cd68", "Adgre1", "Itgam", "Cd14",
  "Lyz2", "Mertk", "Csf1r", "Marco"))
#Values:

fmr[["Neutrophilness"]] <- PercentageFeatureSet(fmr, features = c(
  "Ltf", "Ngp", "Mmp8", "Retnlg",
  "S100a8", "Ltf", "S100a9", "Elane"))
#Values: 0-4

fmr[["Monocyteness"]] <- PercentageFeatureSet(fmr, features = c(
  "F13a1", "Gria3", "Prtn3", "Mpo", 
  "Crip1", "Lgals1", "Tmsb10", "Ly6a"))
#Values: 0.5-2
#===========================================================================

# Specific celltype filters for erythroid lineage
fmr <- NormalizeData(fmr)

# if Seurat is using v5/Assay5: use JoinLayers() and fmr$RNA$data
# otherwise use fmr$RNA@data
fmr[["RNA"]] <- JoinLayers(fmr[["RNA"]])


fmr[["CMPs"]] <- ifelse((fmr$RNA$data["Ly6a",] < 0.5 & fmr$RNA$data["Kit",] > 0.5), 
                        fmr$RNA$data["Cd34",], 0)
fmr[["EPs"]] <- ifelse((
  ((fmr$RNA$data["Ly6a",] + fmr$RNA$data["Itga2b",] + fmr$RNA$data["Cd34",]) < 2) &
    (fmr$RNA$data["Kit",] > 0) & (fmr$RNA$data["Eng",] > 0) &
    ((fmr$RNA$data["Gata1",]  + fmr$RNA$data["Ldb1",] + fmr$RNA$data["Sox6",] + fmr$RNA$data["Eng",]) > 0)),
  (fmr$RNA$data["Eng",] + fmr$RNA$data["Gata1",]  + fmr$RNA$data["Ldb1",] + fmr$RNA$data["Sox6",]), 0)


################ Standard processing ###############

# User input 
num_variable_features <- 15000
cluster_resolution <- 2


fmr <- FindVariableFeatures(fmr, selection.method = "vst", nfeatures = num_variable_features)
fmr.genes <- rownames(fmr)
fmr <- ScaleData(fmr, features = VariableFeatures(fmr))
fmr <- RunPCA(fmr, features = VariableFeatures(fmr))
# Range nVFs from 1500-7000 or use fmr.genes depending on object
ElbowPlot(fmr, 50)

# Use elbowplot crook to decide dims, then adjust to separate clusters meaningfully
neighbors_dims <- 8
umap_dims <- 8
fmr <- FindNeighbors(fmr, dims = 1:neighbors_dims)
fmr <- FindClusters(fmr, resolution = cluster_resolution)
fmr <- RunUMAP(fmr, dims = 1:umap_dims, return.model = F)
dim <- DimPlot(fmr, group.by="ident",reduction = "umap",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
dim2 <- DimPlot(fmr, group.by="subject",reduction = "umap",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
#dim+dim2
feats <- FeaturePlot(fmr, features = c(
  "Hemgn", "Erythroidness", 
  "Hba-a1", "Pf4"), 
  min.cutoff = c(0,1,2,0), max.cutoff = c(1,10,200,1), 
  ncol = 2, pt.size = 2, raster=T)
dim + feats



# Check handpicked erythroid progenitor cells 
handpicked_cell_feats <- FeaturePlot(fmr, features = c(
  "CMPs", "EPs"), 
  min.cutoff = c(0,0), max.cutoff = c(1,1), 
  ncol = 2, pt.size = .5, raster=F)
handpicked_cell_feats


###########################################
# saveRDS(fmr, file = ".rds") #
###########################################


########## Subset erythroid cells (both methods) ###############

clustmarks <- FindMarkers(fmr, ident.1 = 24, min.pct = 0.2, only.pos = T, logfc.threshold = .3)
head(clustmarks, n = 40)



# Identitify clusters of interest and subset 
ery_clusters <- subset(fmr, idents = c())

# Group handpicked cells and subset
fmr$handpicked_celltype <- ifelse(fmr$EPs > 0, "Erythroid Progenitor",
                          ifelse(fmr$CMPs > 0, "Common Myeloid Progenitor","Other"))
handpicked_cells <- subset(fmr, subset = handpicked_celltype == "Other", invert = T)
table(handpicked_cells$handpicked_celltype)
handpicked_cells
cells1 <- colnames(handpicked_cells)
cells2 <- colnames(ery_clusters)
shared_cells <- intersect(x = cells1, y = cells2)
ery_clusters <- subset(ery_clusters, cells = shared_cells, invert = T)
handpicked_cells$selection_method <- "GeneExpression"
ery_clusters$selection_method <- "Cluster"
all_erys <- merge(handpicked_cells, ery_clusters)
table(all_erys$selection_method)
all_erys

# again, JoinLayers() if Seurat object is Assay5 (if all_erys has many layers)
all_erys[["RNA"]] <- JoinLayers(all_erys[["RNA"]])


########## Reprocess Isolated erythroid cells ############

# User input
num_variable_features_erys <- 5000

cluster_resolution_erys <- 0.5


all_erys <- FindVariableFeatures(all_erys, selection.method = "vst", nfeatures = num_variable_features_erys)
all_erys.genes <- rownames(all_erys)
all_erys <- ScaleData(all_erys, features = VariableFeatures(all_erys))# all_erys.genes)
all_erys <- RunPCA(all_erys, features = VariableFeatures(all_erys))# all_erys.genes)
ElbowPlot(all_erys, 50)

neighbors_dims_erys <- 7
umap_dims_erys <- 7
all_erys <- FindNeighbors(all_erys, dims = 1:neighbors_dims_erys)
all_erys <- FindClusters(all_erys, resolution = cluster_resolution_erys)
all_erys <- RunUMAP(all_erys, dims = 1:umap_dims_erys, reduction = "pca", return.model = T)
all_erysdim <- DimPlot(all_erys,raster=F,group.by="ident",reduction = "umap", repel = T, label = T, shuffle = T, pt.size = 2) + NoLegend()
#all_erysdim
# View spectrum of erythoid development
all_erysfeats <- FeaturePlot(Femur, features = c(
  "Cd34", "EPs", "percent.mt", "Bnip3l",
  "Hemgn", "Hba-a1", "nCount_RNA", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
  min.cutoff = c(0,0,2.5,1,
                 0,0,0,.05),
  max.cutoff = c(1,2,3,1.2,
                 2,20,10000,.9))
all_erysfeats


clustmarks <- FindMarkers(all_erys, ident.1 = 1, min.pct = 0.4, only.pos = F, logfc.threshold = .5)
head(clustmarks, n = 40)


########## Compare handpicked to computed celltypes and clean up dataset ##############

# save erythoid subset
# saveRDS(all_erys, file = "")



clustmarks <- FindMarkers(Femur2, ident.1 = 5, min.pct = 0.3, only.pos = F, 
                          logfc.threshold = .5)
head(clustmarks, n = 40)

# Generate heatmap and investigate subpopulations
all_erys.markers <- FindAllMarkers(all_erys, min.pct = 0.25, logfc.threshold = 1)
topmarks <- all_erys.markers %>% group_by(cluster) %>% slice_max(n = 6, order_by = avg_log2FC)

heat = DoHeatmap(all_erys, features = topmarks$gene, size = 4, angle = 10) + NoLegend()
all_erysdim | heat

