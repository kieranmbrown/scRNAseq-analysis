library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scCustomize)
library(stringr)
library(tidyr)
source("calc_3D_UMAPS.R")
make_3d_umap(Femur, 10)

#saveRDS(Femur, "FEMUR_ALL_ERYS_LABELLED.rds")
Femur <- readRDS("FEMUR_RRRM2_RR10_ERYS.rds")

#combined analysis

fmrallGO <- readRDS(file = "FEMUR_GO_ALL.rds")
fmrallGY <- readRDS(file = "FEMUR_GY_ALL.rds")
fmrallFO <- readRDS(file = "FEMUR_FO_ALL.rds")
fmrallFY <- readRDS(file = "FEMUR_FY_ALL.rds")
Femur <- merge(fmrallGO, c(fmrallGY, fmrallFO, fmrallFY))

fmrGO <- readRDS(file = "FEMUR_GO_redo_erys.rds")
fmrGY <- readRDS(file = "FEMUR_GY_redo_erys.rds")
fmrFO <- readRDS(file = "FEMUR_FO_redo_erys.rds")
fmrFY <- readRDS(file = "FEMUR_FY_redo_erys.rds")

fmrGO$Cohort <- "GO"
fmrGY$Cohort <- "GY"
fmrFO$Cohort <- "FO"
fmrFY$Cohort <- "FY"

Femur <- subset(Femur, subset=Sample_num==4)
Femur <- merge(fmrGO, c(fmrGY, fmrFO, fmrFY))
table(Femur$subject)



fmr1 <- readRDS("FEMUR_ALL_ERYS_LABELLED.rds")
fmr2 <- readRDS("FEMUR_IFO_redo_erys.rds")
fmr2$Cohort <- "SO"
fmr2$age <- "Old"
fmr2$experiment <- "Space"
fmr2$subject <- ifelse(fmr2$subject == "IFO3", "SO3",fmr2$subject)
Femur <- merge(fmr1, fmr2)

VlnPlot(Femur, features= c("Csf1"),
		pt.size=0,ncol = 1)

FeatureScatter(Femur, feature1="Ighm", feature2="Hba-a1",group.by="Cohort")


table(AverageExpression(Femur, group.by="Cohort",features =c("Ighm")))

AverageExpression(Femur, group.by="subject",features =c("Ighm"))
fmrsave <- Femur


num_variable_features <- 5000

Femur <- FindVariableFeatures(Femur, selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(Femur)
Femur <- ScaleData(Femur, features = VariableFeatures(Femur))
Femur <- RunPCA(Femur, features = VariableFeatures(Femur))
ElbowPlot(Femur, 50)



cluster_resolution <- .2
neighbors_dims <- 30
umap_dims <- 30
Femur <- FindNeighbors(Femur, dims = 1:neighbors_dims)
Femur <- FindClusters(Femur, resolution = cluster_resolution)
Femur <- RunUMAP(Femur, dims = 1:umap_dims)

dim <- DimPlot(Femur, group.by="ident",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
feats <- FeaturePlot(Femur, features = c(
		"Hmgb2", "Epor", 
		"Hba-a1", "Bpgm"), 
		min.cutoff = c(3,0,0,0), max.cutoff = c(5,1,8,2), 
		ncol = 2, pt.size = 2, raster=T)
dim + feats


all_erysfeats <- FeaturePlot(Femur, features = c(
  "Cd34", "Hmgb2", "Epor", "Bpgm",
  "Hemgn", "Hba-a1", "Gypa", "RNA.diversity"),	 pt.size = 2,raster = T,ncol=4, 
  min.cutoff = c(0,3,0,0,
                 0,0,0,.05),
  max.cutoff = c(1,5,1,1,
                 4,20,2,.9))
all_erysfeats

all_erysfeats <- FeaturePlot(Femur, features = c(
  "Pik3ca", "Pik3c2a", "Pik3c3", "Pik3r1",
  "Mtor", "Prkdc", "Atm", "RNA.diversity"),	 pt.size = 2,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(1,1,1,1,
                 1,1,1,1))
all_erysfeats


feats <- FeaturePlot(Femur, features = c(
  "Macrophageness", "nFeature_RNA", 
  "Monocyteness", "nCount_RNA"), 
  min.cutoff = c(0.5,50,0.5,100), max.cutoff = c(2,4000,2,15000), 
  ncol = 2, pt.size = 2, raster=T)
dim + feats


dim | DimPlot(Femur, group.by="Celltype",label = T, pt.size = 1, shuffle = T, raster = T)


head(FindMarkers(fmr1, ident.1 = "1_Pro-Erythroblast", min.pct = .2, only.pos = T, 
                          logfc.threshold = 1), n = 40)

saveRDS(Femur, "FEMUR_WHOLE_LABELLED.rds")


Femur$Celltype <- ifelse(Idents(Femur) == 11 | Idents(Femur) == 22, "Erythroids", 
                         ifelse(Idents(Femur) == 4 | Idents(Femur) == 5 | Idents(Femur) == 6, "Immature B cells", 
                                ifelse(Idents(Femur) == 21 | Idents(Femur) == 15, "Pro-B cells", 
                                ifelse(Idents(Femur) == 17, "HSCs", 
                                ifelse(Idents(Femur) == 14 | Idents(Femur) == 8, "T cells", 
                                ifelse(Idents(Femur) == 20, "Adipocytes", 
                                ifelse(Idents(Femur) == 2 | Idents(Femur) == 3 | Idents(Femur) == 12 | Idents(Femur) == 7 | Idents(Femur) == 18, "Neutrophils", 
                                       "Other")))))))


Femur$Celltype <- ifelse(Femur$CMPs > 0 | Femur$EPs > 0, "Common Myeloid Progenitor", "Other")
"""
ifelse(Femur$RNA@data["Pf4",] > 4, "Platelets",
ifelse(Femur$RNA@data["Jchain",] > 4 | Idents(Femur) == 23, "Plasma cells",
ifelse(Femur$CMPs > 0 | Femur$EPs > 0, "Common Myeloid Progenitor",
ifelse(Idents(Femur) == 1 | Idents(Femur) == 0 | Idents(Femur) == 9 | Idents(Femur) == 10, "Monocytes and Macrophages",
ifelse(Idents(Femur) == 16 | Idents(Femur) == 13 | Idents(Femur) == 19, "Plasmacytoid Dendritic cells",
       """              


Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Femur$CMPs > 0, "Common Myeloid Progenitor", Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Femur$RNA@data["Pf4",] > 4, "Platelets", Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Femur$RNA@data["Jchain",] > 4 | Idents(Femur) == 23, "Plasma cells", Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Femur$RNA@data["Siglech",] > 5, "Plasmacytoid Dendritic Cells", Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Idents(Femur) == 9 | Idents(Femur) == 10,  "Monocytes",  Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Idents(Femur) == 1, "Common Monocyte Progenitor",  Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" & Idents(Femur) == 0, "Mononuclear Phagocyte Progenitors", Femur$Celltype)
Femur$Celltype <- ifelse(Femur$Celltype == "Other" | Femur$Celltype == "NA", "Mononuclear Phagocyte Progenitors", Femur$Celltype)


Femur$Celltype <- ifelse(Idents(Femur) == 17, "HSCs", "Other")


Idents(Femur) == 1 | 



FeaturePlot(Femur, features = c(
  "Ccr2", "Cd34", "Cx3cr1", "Elane",
  "Cd68", "Mpo", "Flt3", "Sirpa"),	 pt.size = 1,raster = T,ncol=4, 
  min.cutoff = c(0,0,0,0,
                 0,0,0,0),
  max.cutoff = c(2,2,2,2,
                 2,2,2,2))















table(Femur$Cohort,Idents(Femur))

gene <- "Bank1"  

# Fetch expression data for the gene
data <- FetchData(Femur, vars = c(gene, "ident"))  # 'ident' corresponds to clusters or cell types

# Create the plot with ggplot2
ggplot(data, aes(x = ident, y = gene)) +
  geom_violin(scale = "count") +  # 'count' scales violins based on the number of observations
  labs(x = "Cell Type/Cluster", y = "Expression Level") +
  theme_minimal()



Plot_Density_Custom(seurat_object = Femur, features = c("Hemgn", "Hba-a1"))


FeaturePlot_scCustom(seurat_object = Femur, features = "Hba-a1")



Femur[["CMPs"]] <- ifelse((Femur$RNA@data["Ly6a",] < 0.5 & Femur$RNA@data["Kit",] > 0.5), 
                        Femur$RNA@data["Cd34",], 0)
Femur[["EPs"]] <- ifelse((
  ((Femur$RNA@data["Ly6a",] + Femur$RNA@data["Itga2b",] + Femur$RNA@data["Cd34",]) < 2) &
    (Femur$RNA@data["Kit",] > 0) & (Femur$RNA@data["Eng",] > 0) &
    ((Femur$RNA@data["Gata1",]  + Femur$RNA@data["Ldb1",] + Femur$RNA@data["Sox6",] + Femur$RNA@data["Eng",]) > 0)),
  (Femur$RNA@data["Eng",] + Femur$RNA@data["Gata1",]  + Femur$RNA@data["Ldb1",] + Femur$RNA@data["Sox6",]), 0)


Femur$Erythroids <- ifelse(Idents(Femur) == 11 | Idents(Femur) == 22 | Femur$CMPs > 0 | Femur$EPs > 0, "Erythroid", "Other")
table(Femur$subject, Femur$Erythroids)

erythroids_umap <- DimPlot(Femur, cols=c("red","lightgrey"),group.by="Erythroids",reduction = "umap", pt.size = 1, shuffle = T, raster = T) + NoLegend()
erythroids_umap




cohortscatter <- function(obj,gene){
  #obj <- subset(Femur,idents="3_Orthoblast")
  
  expressiondf <- as.data.frame(AverageExpression(obj, group.by="subject",features =c(gene)))
  long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
  
  long_data <- long_data %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.SO") ~ "In-Space Old",
      str_detect(Subject, "\\.FO") ~ "Flight Old",
      str_detect(Subject, "\\.FY") ~ "Flight Young",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.GY") ~ "Ground Young"
    ))
  
  ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort, size=1.5)) +
    geom_point() +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene,"Expression by Subject and Cohort"),
         x = "Subject ID",
         y = paste(gene,"Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0, NA))
}


celltypescatter <- function(obj, gene) {
  # Get average expression grouped by cell type
  expressiondf <- as.data.frame(AverageExpression(obj, group.by = "ident", features = c(gene)))
  long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "CellType", values_to = "gene")
  
  long_data$CellType <- factor(long_data$CellType, levels = c("g3_Orthoblast", "g1_Pro-Erythroblast", "g0_Erythroid Progenitor", "g2_Erythroblast"))
  
  ggplot(long_data, aes(x = CellType, y = gene, color = CellType, size = 1.5)) +
    geom_point() +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression by Cell Type"),
         x = "Cell Type",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, NA))
}



cohortscatter <- function(obj, gene) {
  # Get average expression grouped by cell type
  expressiondf <- as.data.frame(AverageExpression(obj, group.by = "ident", features = c(gene)))
  long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "CellType", values_to = "gene")
  
  # Ensure CellType is a factor
  long_data$CellType <- factor(long_data$CellType, levels = c("3_Orthoblast", "1_Pro-Erythroblast", "0_Erythroid Progenitor", "2_Erythroblast"))
  
  # Generate the plot
  ggplot(long_data, aes(x = CellType, y = gene, color = CellType)) +
    geom_point(size = 1.5) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression by Cell Type"),
         x = "Cell Type",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, NA))
}

conditionscatter <- function(obj, gene) {
  expressiondf <- as.data.frame(AverageExpression(obj, layer= "counts",assay="SCT",group.by = "subject", features = c(gene)))
  long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
  
  # define cohort by subject info
  long_data <- long_data %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "\\.SY") ~ "RR10 FLIGHT",
      str_detect(Subject, "\\.GY") ~ "RRRM2 Ground",
      str_detect(Subject, "\\.FY") ~ "RRRM2 Return",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.FO") ~ "Flight Old"
    ))
  
  # plot
  ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Subject), hjust=-.3,size = 4) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression by Subject and Cohort"),
         x = "Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0, NA))
}


conditionscatter <- function(obj, gene) {
  expressiondf <- as.data.frame(AverageExpression(obj, layer= "counts",assay="SCT",group.by = "subject", features = c(gene)))
  long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
  
  # define cohort by subject info
  long_data <- long_data %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "\\.SY") ~ "RR10 FLIGHT",
      str_detect(Subject, "\\.GY") ~ "RRRM2 Ground",
      str_detect(Subject, "\\.FY") ~ "RRRM2 Return",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.FO") ~ "Flight Old"
    ))
  
  # plot
  plot1 <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Subject), hjust=-.3,size = 4) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression by Subject and Cohort"),
         x = "Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0, NA))
  
  plot2 <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort)) +
    geom_point(size = 3) +
    geom_text_repel(aes(label = Subject), hjust=-.3,size = 4) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression by Subject and Cohort"),
         x = "Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits = c(0, NA))
  plot1 | plot2
}






table(Femur$Cohort,Idents(Femur))


all_erysfeats <- FeaturePlot(Femur, features = c(
		"Cd34", "Cd44", "Gata1", "Bpgm",
		"Ighm", "Hba-a1", "Tfrc", "RNA.diversity"),	 pt.size = 4,raster = T,ncol=4, 
		min.cutoff = c(0,0,0,0,
				          0,0,0,.05),
		max.cutoff = c(1,2,1,2,
				          4,20,2,.9))
all_erysfeats


celltypes <- c(
		"3_Orthoblast",
		"3_Orthoblast",
		"1_Pro-Erythroblast", 
		"3_Orthoblast",
		"0_Erythroid Progenitor",
		"2_Erythroblast",
		"3_Orthoblast",
		"3_Orthoblast",
		"0_Erythroid Progenitor")
names(celltypes) <- levels(Femur)
Femur <- RenameIdents(Femur, celltypes)

Femur[["celltype"]] <- Idents(Femur)
cellnums <- c(0,1,2,3,4)
names(cellnums) <- levels(Femur)
Femur <- RenameIdents(Femur, cellnums)
Idents(Femur) <- Femur$celltype

Femur$experiment <- ifelse(Femur$Cohort == "FO" | Femur$Cohort == "FY", "Flight", "Ground")
Femur$age <- ifelse(Femur$Cohort == "GY" | Femur$Cohort == "GO", "Ground", Femur$Cohort)


young <- subset(Femur, subset= age == "Old", invert = T)
youngortho <- subset(young,idents = c("4_Ortho-Erythroblast","3_Poly-Erythroblast"))

ct_table <- table(young$Cohort, Idents(young))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","orange","blue","lightblue"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Cell Type Proportions in Each Cohort")

ct_table <- table(Femur$experiment, Idents(Femur))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","blue","green"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Proportions of Cell Types in Flight/Ground")

ct_table <- table(Femur$subject, Idents(Femur))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","red","red","red","orange","orange","orange","orange","blue", "blue","blue","blue","lightblue","lightblue","lightblue","lightblue"), 
                legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Proportions of Cell Types in Flight Young/Old vs Ground")

table(Femur$Sample_num)
Femur$Sample_num <- as.numeric(str_extract(Femur$subject, "\\d{1,2}$"))

VlnPlot(Femur, features=c("nCount_RNA","Hba-a1", "RNA.diversity"),group.by="experiment",ncol=3,pt.size=0)/
VlnPlot(Spleen, features=c("nCount_RNA","Hba-a1","RNA.diversity"),group.by="experiment", ncol=3,pt.size=0)


ct_table <- table(Femur$Cohort, Idents(Femur))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","orange","blue","lightblue", "green"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Cell Type Proportions in Each Cohort")




clustmarks <- FindMarkers(Femur, ident.1 = 2, ident.2=6, min.pct = 0.2, only.pos = T, logfc.threshold = .5)
head(clustmarks, n = 40)

clustmarks <- FindMarkers(Femur, ident.1 = 5, min.pct = 0.2, only.pos = F, logfc.threshold = .5)
head(clustmarks, n = 40)




clustmarks <- FindMarkers(Femur, ident.1 = 0, ident.2=2, min.pct = 0.2, only.pos = F, logfc.threshold = .5)
head(clustmarks, n = 40)

num_variable_features <- 1250
neighbors_dims <- 10
cluster_resolution <- .2
umap_dims <- 10

Spleen <- FindVariableFeatures(Spleen, selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(Spleen)
Spleen <- ScaleData(Spleen, features = VariableFeatures(Spleen))
Spleen <- RunPCA(Spleen, features = VariableFeatures(object = Spleen))
ElbowPlot(Spleen, 50)
Spleen <- FindNeighbors(Spleen, dims = 1:neighbors_dims)
Spleen <- FindClusters(Spleen, resolution = cluster_resolution)
Spleen <- RunUMAP(Spleen, dims = 1:umap_dims, return.model = F)
spldim <- DimPlot(Spleen, group.by="ident",reduction = "umap",label = T, pt.size = .75, shuffle = T, raster = T) + NoLegend()
splfeats <- FeaturePlot(Spleen, features = c(
		"Hemgn", "Ighm", 
		"nCount_RNA", "Gypa"), 
		min.cutoff = c(0,0,0,0), max.cutoff = c(1,1,1000,1), 
		ncol = 2, pt.size = 2, raster=T)
spldim + splfeats


FeatureScatter(Spleen, feature1="Ighm",feature2="Gypa",raster=T)
VlnPlot(Spleen, features=c("Ighm","Gypa","Hba-a1","Hemgn"),ncol=2,pt.size=0)



splcelltypes <- c(
		"Marked inactive", 
		"Marked active", 
		"Unmarked active",
		"Unmarked inactive",
		"Unmarked inactive",
		"Unmarked inactive",
		"Marked inactive", 
		"Marked inactive", 
		"Unmarked active",
		"Unmarked inactive",
		"Unmarked inactive",
		"Unmarked inactive",
		"Unmarked inactive")
names(splcelltypes) <- levels(Spleen)
Spleen <- RenameIdents(Spleen, splcelltypes)



Femur$Sample_num <- ifelse((Femur$subject == "FO19" | 
                             Femur$subject == "FY14" |
                             Femur$subject == "GO19" |
                             Femur$subject == "GY9"), 4, 
                           ifelse((Femur$subject == "FO16" | 
                                    Femur$subject == "FY11" |
                                    Femur$subject == "GO16" |
                                    Femur$subject == "GY6"), 3, 
                                  ifelse((Femur$subject == "FO4" | 
                                           Femur$subject == "FY9" |
                                           Femur$subject == "GO14" |
                                           Femur$subject == "GY4"), 2, 1)))


ct_table <- table(Spleen$Cohort, Idents(Spleen))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","orange","blue","lightblue"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Cell Type Proportions in Each Cohort")

ct_table <- table(Spleen$experiment, Idents(Spleen))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","blue"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Proportions of Cell Types in Flight/Ground")

ct_table <- table(Spleen$subject, Idents(Spleen))
ct_proportions <- prop.table(ct_table, margin = 1)
bar <- barplot(ct_proportions, beside = T, col= c("red","orange","blue"), legend.text = TRUE, args.legend = list(x = "topright"))
title(main = "Proportions of Cell Types in Flight Young/Old vs Ground")

clustmarks <- FindMarkers(Femur, ident.1 = "1", min.pct = 0.2, only.pos = T, logfc.threshold = 1.7,group.by="experiment")
head(clustmarks, n = 5)


DotPlot(Femur, features = c("Hemgn", "Tfrc","Ighm","Bpgm",
                            "Hba-a2","Top2a",
                            "Gata1"), 
        group.by="Cohort",split.by="ident", cols=c("red","orange","blue","lightblue")) + RotatedAxis()


VlnPlot(Spleen, features=c("Ighm","Gypa"), ncol=2,pt.size=0)

VlnPlot(Femur, features=c("Ighm", "Tfrc"),
        split.by="Cohort",group.by="", 
        ncol=2,pt.size=0)



fmrvar[["Erythroidness"]] <- PercentageFeatureSet(fmrvar, features = c(
  "Hba-a1", "Hba-a2", "Hbb-bs", "Hbb-bt",
  "Hemgn", "Ermap", "Slc25a21", "Slc4a1"))
#Values: 0.2-40




fmrvar <- fmrGO

num_variable_features <- 3500
cluster_resolution <- 1.2



fmrvar <- FindVariableFeatures(fmrvar, selection.method = "vst", nfeatures = num_variable_features)
all.genes <- rownames(fmrvar)
fmrvar <- ScaleData(fmrvar, features = VariableFeatures(fmrvar))
fmrvar <- RunPCA(fmrvar, features = VariableFeatures(object = fmrvar))
#ElbowPlot(fmrvar, 50)


neighbors_dims <- 12
umap_dims <- 12
fmrvar <- FindNeighbors(fmrvar, dims = 1:neighbors_dims)
fmrvar <- FindClusters(fmrvar, resolution = cluster_resolution)
fmrvar <- RunUMAP(fmrvar, dims = 1:umap_dims, return.model = F)
dim <- DimPlot(fmrvar, group.by="ident",reduction = "umap",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
dim2 <- DimPlot(fmrvar, group.by="subject",reduction = "umap",label = T, pt.size = 1, shuffle = T, raster = T) + NoLegend()
#dim+dim2
feats <- FeaturePlot(fmrvar, features = c(
  "Hemgn", "Erythroidness", 
  "Hba-a1", "Pf4"), 
  min.cutoff = c(0,1,2,0), max.cutoff = c(1,20,200,1), 
  ncol = 2, pt.size = 2, raster=F)
dim + feats


splitdim <- DimPlot(fmrvar, group.by="ident",split.by="subject",
                    reduction = "umap",label = T, pt.size = .75, 
                    shuffle = T, raster = F) + NoLegend()
splitdim










fmrGY[["RNA"]] <- split(fmrGY[["RNA"]], f = fmrGY$subject)

fmrGY <- SCTransform(fmrGY)
fmrGY <- RunPCA(fmrGY)
DefaultAssay(fmrGY) <- "SCT"

fmrGY <- IntegrateLayers(object = fmrGY, method = CCAIntegration, orig.reduction = "pca",
                         new.reduction = "integrated.cca", normalization.method = "SCT", verbose=F, k.weight=43)

# re-join layers after integration
fmrGY[["RNA"]] <- JoinLayers(fmrGY[["RNA"]])

fmrGY <- FindNeighbors(fmrGY, reduction = "integrated.cca", dims = 1:8)
fmrGY <- FindClusters(fmrGY, resolution = 1)

fmrGY <- RunUMAP(fmrGY, dims = 1:15, reduction = "integrated.cca")

DimPlot(fmrGY, reduction = "umap", group.by = c("subject", "ident"))







df <- df %>%
  group_by(Cohort) %>%
  mutate(TotalCells = sum(nCells),            # Total cells per cohort
         Proportion = nCells / TotalCells)    # Proportion of each cell type

# Create the bar plot
ggplot(df, aes(x = CellType, y = Proportion, fill = Cohort)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_y_continuous(labels = scales::percent_format()) +  # Convert y-axis to percentage
  labs(title = "Proportion of Cell Types by Cohort",
       x = "Cell Type",
       y = "Proportion (%)") +
  theme_minimal() +
  guides(fill = guide_legend(title = "Cohort"))

# Optionally save the plot
ggsave("normalized_celltype_distribution.png", width = 8, height = 6)







EPs_data <- FetchData(Femur, vars = c("Hemgn", "Epor","Bpgm","Hmgb2", "Hba-a1", "ident","nCount_RNA"),layer="data")

EPs_scatter <- ggplot(data = EPs_data) +
  geom_point(mapping = aes(x = Epor, y = Hmgb2, color = ident, size = nCount_RNA, alpha = .5), position = 'jitter') #+
  #scale_color_gradientn(colors = c("blue", "red"))
EPs_scatter
















s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


s.genes <- sapply(s.genes, function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))
g2m.genes <- sapply(g2m.genes, function(x) paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x)))))


marrow <- Femur

Femur$synthesis <- PercentageFeatureSet(Femur, features = head(s.genes,n = 10))

Femur <- CellCycleScoring(Femur, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


RidgePlot(marrow, assay="RNA",features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2,group.by = "Phase")

DimPlot(marrow, group.by="Phase",reduction = "umap",label = T, pt.size = .5, shuffle = T, raster = F) + NoLegend()

s.marrow <- subset(marrow, subset = Phase == "S")


DimPlot(Femur, group.by = "Phase",split.by = "Cohort", 
        raster=T,ncol=3, pt.size = 3)

FeaturePlot(Femur, features = 'synthesis', split.by = "Cohort",ncol = 3 ,raster=T, pt.size = 2.5)
