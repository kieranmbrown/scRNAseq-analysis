#################################
## SCT NORMALIZATION PROCEDURE ##
#################################

#spln <- readRDS("spln_WHOLE_LABELLED.rds")
#spln <- readRDS("spln_ALL_ERYS_LABELLED.rds")
spln <- readRDS("SPLEEN_ALL_ERYS_UMAPED2.rds")




cohort_split <- SplitObject(spln2, split.by = "Cohort")

output_cohorts <- list()

for (i in seq_along(cohort_split)) {
  
  cohort <- cohort_split[[i]]
  
  cohort[["RNA"]] <- split(cohort[["RNA"]], f = cohort$subject)
  
  cohort <- SCTransform(cohort)
  cohort <- RunPCA(cohort)
  DefaultAssay(cohort) <- "SCT"
  
  cohort <- IntegrateLayers(object = cohort, method = CCAIntegration, orig.reduction = "pca",
                            new.reduction = "integrated.cca", normalization.method = "SCT", verbose=F, k.weight=43)
  
  # re-join layers after integration
  cohort[["RNA"]] <- JoinLayers(cohort[["RNA"]])
  
  output_cohorts[[i]] <- cohort
}



cohort <- spln

cohort[["RNA"]] <- split(cohort[["RNA"]], f = cohort$Cohort)

cohort <- SCTransform(cohort)
cohort <- RunPCA(cohort)
DefaultAssay(cohort) <- "SCT"

options(future.globals.maxSize=734003200)
cohort <- IntegrateLayers(object = cohort, method = CCAIntegration, orig.reduction = "pca",
            new.reduction = "integrated.cca", normalization.method = "SCT", verbose=F, k.weight=43)

# re-join layers after integration
cohort[["RNA"]] <- JoinLayers(cohort[["RNA"]])


spln2 <- cohort

cluster_resolution <- .4
neighbors_dims <- 5
umap_dims <- 5
spln2 <- FindNeighbors(spln2, dims = 1:neighbors_dims, reduction = "integrated.cca")
spln2 <- FindClusters(spln2, resolution = cluster_resolution)
spln2 <- RunUMAP(spln2, dims = 1:umap_dims, reduction = "integrated.cca")
#ElbowPlot(spln2, ndims = 50)

dim <- DimPlot(spln, group.by="ident",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
dim2 <- DimPlot(spln2, group.by="Cohort",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
feats <- FeaturePlot(spln2, features = c(
  "Cd34", "Hmgb2", 
  "Hba-a1", "Bpgm"), 
  min.cutoff = c(0,1,1,1), max.cutoff = c(1,3,10,4), 
  ncol = 2, pt.size = 2, raster=T, reduction = "umap")
dim + feats



all_erysfeats <- FeaturePlot(spln2, features = c(
  "Cd34", "Hmgb2", "Epor", "Bpgm",
  "Hemgn", "Hba-a1", "Gypa", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
  min.cutoff = c(0,1,0,0,
                 0,0,0,.05),
  max.cutoff = c(1,3,2,4,
                 4,20,2,.9))
all_erysfeats

clustmarks <- FindMarkers(spln3,assay = "RNA", ident.2 = 2, ident.1=c(1,0), 
                          min.pct = 0.2, only.pos = T, logfc.threshold = .3)
head(clustmarks, n = 100)


spln3 <- subset(spln3, idents = c(6), invert = T)

num_variable_features <- 1500

spln3 <- FindVariableFeatures(spln3, assay = "SCT",selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(spln3)
spln3 <- ScaleData(spln3, features = VariableFeatures(spln3)``)
spln3 <- RunPCA(spln3, features = VariableFeatures(spln3))

neighbors_dims = 5
umap_dims = 5
cluster_resolution = .2

spln3 <- FindNeighbors(spln3, dims = 1:neighbors_dims, reduction = "integrated.cca")
spln3 <- FindClusters(spln3, resolution = cluster_resolution)
spln3 <- RunUMAP(spln3, dims = 1:umap_dims, reduction = "integrated.cca")
#ElbowPlot(spln3, ndims = 50)

dim <- DimPlot(spln3, group.by="ident",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
dim2 <- DimPlot(spln3, group.by="Cohort",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
feats <- FeaturePlot(spln3, features = c(
  "Cd34", "Hmgb2", 
  "Hba-a1", "Bnip3l"), 
  min.cutoff = c(0,1,1,1), max.cutoff = c(1,3,10,2), 
  ncol = 2, pt.size = 2, raster=T, reduction = "umap")
dim + feats

all_erysfeats <- FeaturePlot(spln3, features = c(
  "Cd34", "Hmgb2", "Epor", "Bpgm",
  "Hemgn", "Hba-a1", "Gypa", "RNA.diversity"),	 pt.size = 2,raster = T,ncol=4, 
  min.cutoff = c(0,1,.4,0,
                 0,0,0,.05),
  max.cutoff = c(1,2,.7,4,
                 1.5,20,2,.8))
all_erysfeats

all_erysfeats <- FeaturePlot(spln3, features = c(
  "Cd34", "Klf1", "Epor", "Gypa",
  "Alas2", "Slc4a1", "percent.rp", "Bnip3l"),	 pt.size = 2,raster = T,ncol=4, 
  min.cutoff = c(0,0,.4,0,
                 0,1.5,0.6,1),
  max.cutoff = c(1,1,.7,2,
                 3,3,1,2))
all_erysfeats

spln3 <- PrepSCTFindMarkers(spln3)

VlnPlot(spln3, pt.size = 0,features = c("percent.rp","percent.rb"))

spln3[["percent.rb"]] <- PercentageFeatureSet(spln3, pattern =  "^Rb[sl]")
spln3[["percent.rp"]] <- PercentageFeatureSet(spln3, pattern =  "^Rp[sl]")


ct_table <- table(subset(spln3, subset = age == "Young")$Cohort, Idents(subset(spln3, subset = age == "Young")))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","orange","blue","lightblue"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Cell Type Proportions in Each Cohort")


space_diff <- ct_proportions["SO",]-ct_proportions["RR10GC",]
return_diff <- ct_proportions["FY",]-ct_proportions["GY",]
diff_mat <- rbind(space_diff,return_diff)
diff_bar <- barplot(diff_mat, beside=T,col=c("blue","red"),
                    legend.text = TRUE, args.legend = list(x = "topright"))

celltypes <- c(
  "4_Ortho-Erythroblast", 
  "4_Ortho-Erythroblast",
  "4_Ortho-Erythroblast",
  "2_Baso-Erythroblast",
  "1_Pro-Erythroblast", 
  "0_Erythroid Progenitor",
  "3_Poly-Erythroblast")
names(celltypes) <- levels(spln3)
spln3 <- RenameIdents(spln3, celltypes)


heatmap.markers <- FindAllMarkers(spln2markers, min.pct = 0.20, logfc.threshold = .8, only.pos = T)
topmarks <- heatmap.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)

heat = DoHeatmap(spln2markers, features = topmarks$gene, size = 4, angle = 10) + NoLegend()
dim | heat


table(spln3$experiment)


ct_table





saveRDS(spln3, "spln_RRM2_RR210_ERYS.rds")