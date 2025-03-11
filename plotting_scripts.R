library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(scCustomize)
library(stringr)
library(tidyr)
library(Biobase)
library(BiocManager)
library(BiocGenerics)
library(BiocNeighbors)
library(BiocParallel)
library(BiocVersion)
library(ggrepel)
library(cowplot)
library(presto)
library(future)
library(future.apply)
#library(DescTools)
plan("multisession")
Femur <- readRDS("FEMUR_RRRM2_RR10_ERYS.rds")
Femur$sample_number <- as.numeric(str_extract(Femur$subject, "\\d+$"))




Femur$sample_rank <- Femur$sample_number

for(cohort in unique(Femur$Cohort)){
  cohort_sub <- subset(Femur, subset = Cohort == cohort)
  number <- 0
  for(sample in unique(cohort_sub$subject)){
    number <- number + 1
    Femur$sample_rank <- ifelse(Femur$subject == sample, number, Femur$sample_rank)
  }
}


VlnPlot(Femur, features="nFeature_RNA", pt.size=0,group.by="sample_rank", split.by = "Cohort")



Femur$Cohort <- ifelse(Femur$Cohort == "SO", "RR10-FLT", Femur$Cohort)
Femur$Cohort <- ifelse(Femur$Cohort == "RR10GC", "RR10-GC", Femur$Cohort)


##############################
## Completed Plotting Tools ##
##############################


##### Average expression of Gene across cohorts by subject

expression_by_celltypes <- function(obj, layer_to_use, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  
  # gather data for pooled dataset
  expressiondf_all <- as.data.frame(AverageExpression(obj, layer= layer_to_use, assay="RNA", group.by = "subject", features = c(gene)))
  long_data_all <- pivot_longer(expressiondf_all, cols = everything(), names_to = "Subject", values_to = "gene")
  
  # define cohort for the entire dataset
  long_data_all <- long_data_all %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "\\.SY") ~ "RR10 Flight",
      str_detect(Subject, "\\.GY") ~ "Ground Young",
      str_detect(Subject, "\\.FY") ~ "Return Young",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.FO") ~ "Return Old"
    ))
  long_data_all$Cohort <- factor(long_data_all$Cohort, levels = cohort_order)
  
  # generate plot for all cell types
  plot_all <- ggplot(long_data_all, aes(x = Cohort, y = gene, color = Cohort, fill = Cohort)) +
    geom_point(size=3.5, position = position_dodge(width=.75)) +
    geom_boxplot(alpha=0.5,) +
    geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression Across All Cell Types"), subtitle="by Subject and Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
    scale_y_continuous(limits = c(0, NA))
  
  # add plot for all cell types to plots list
  plots[["All Cell Types"]] <- plot_all
  
  # loop over each cell type for specific plots
  for (cell_type in cell_types) {
    
    sub_obj <- subset(obj, idents = cell_type)
    
    # compute average expression for the subset
    expressiondf <- as.data.frame(AverageExpression(sub_obj, layer= layer_to_use, assay="RNA", group.by = "subject", features = c(gene)))
    long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
    
    # assign cohorts to subjects
    long_data <- long_data %>%
      mutate(Cohort = case_when(
        str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
        str_detect(Subject, "\\.SY") ~ "RR10 Flight",
        str_detect(Subject, "\\.GY") ~ "Ground Young",
        str_detect(Subject, "\\.FY") ~ "Return Young",
        str_detect(Subject, "\\.GO") ~ "Ground Old",
        str_detect(Subject, "\\.FO") ~ "Return Old"
      ))
    long_data$Cohort <- factor(long_data$Cohort, levels = cohort_order)
    
    # generate plot for current cell type
    plot <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort, fill = Cohort)) +
      geom_point(size=3.5, position = position_dodge(width=.75)) +
      geom_boxplot(alpha=0.5) +
      geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
      theme_minimal() + NoLegend() +
      labs(title = paste(gene, "Expression in", cell_type), subtitle="by Subject and Cohort",
           y = paste(gene, "Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
      scale_y_continuous(limits = c(0, NA))
    plots[[cell_type]] <- plot
  }
  plot_names <- sort(names(plots))
  sorted_plots <- plots[plot_names]
  
  # combine all plots
  plot_grid(plotlist = sorted_plots, align = "v", ncol = 3)
}




#### Quick expression across all cell types for each cohort by subject

expression_pooled <- function(obj, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  
  expressiondf_all <- as.data.frame(AverageExpression(obj, assay="RNA", group.by = "subject", features = c(gene)))
  long_data_all <- pivot_longer(expressiondf_all, cols = everything(), names_to = "Subject", values_to = "gene")
  
  # Assign cohorts
  long_data_all <- long_data_all %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "\\.SY") ~ "RR10 Flight",
      str_detect(Subject, "\\.GY") ~ "Ground Young",
      str_detect(Subject, "\\.FY") ~ "Return Young",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.FO") ~ "Return Old"
    ))
  long_data_all$Cohort <- factor(long_data_all$Cohort, levels = cohort_order)

  # Generate plot for all cell types
  ggplot(long_data_all, aes(x = Cohort, y = gene, color = Cohort, fill = Cohort)) +
    geom_point(size=3.5, position = position_dodge(width=.75)) +
    geom_boxplot(alpha=0.5) +
    geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression in"), subtitle="by Subject and Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
    scale_y_continuous(limits = c(0, NA))
  
}




#### Fraction of cells belonging to each cell type for each cohort by subject 

celltype_proportions <- function(obj) {
  # Get unique cell types
  cell_types <- unique(Idents(obj))
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  
  
  # Create a table with fraction of cells in each type per subject, convert to df for ggplot
  proportions_df <- as.data.frame(prop.table(table(obj$subject, Idents(obj)), margin = 1))
  names(proportions_df) <- c("Subject", "CellType", "Fraction")
  
  # Assign cohorts
  proportions_df <- proportions_df %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "SY") ~ "RR10 Flight",
      str_detect(Subject, "GY") ~ "Ground Young",
      str_detect(Subject, "FY") ~ "Return Young",
      str_detect(Subject, "GO") ~ "Ground Old",
      str_detect(Subject, "FO") ~ "Return Old"
    ))
  # Specify orxer for comparison
  proportions_df$Cohort <- factor(proportions_df$Cohort, levels = cohort_order)
  print(head(proportions_df))
  # Plot with whisker box plot 
  ggplot(proportions_df, aes(x = Cohort, y = Fraction, fill = CellType, color = CellType)) +
    geom_point(size=1.5, position = position_dodge(width=.75)) +
    geom_boxplot(alpha=0.5) +
    theme_minimal() +
    labs(title = "Relative Abundance of Each Cell Type by Cohort",
         x = "Cohort",
         y = "Fraction of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, NA))
}

celltype_proportion_comparison <- function(obj) {
  # Get unique cell types
  cell_types <- unique(Idents(obj))
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  celltype_order <- c("0_Erythroid Progenitor","1_Pro-Erythroblast","2_Baso-Erythroblast","3_Poly-Erythroblast","4_Ortho-Erythroblast")
  my_colors <- c("RR10 Ground" = "#B2182B", "RR10 Flight" = "#EF8A62", 
                 "Ground Young" = "#1B7837", "Return Young" = "#7FBF7B", 
                 "Ground Old" = "#2166AC", "Return Old" = "#92C5DE")
  # Create a table with fraction of cells in each type per subject, convert to df for ggplot
  proportions_df <- as.data.frame(prop.table(table(obj$subject, Idents(obj)), margin = 1))
  names(proportions_df) <- c("Subject", "CellType", "Fraction")
  
  # Assign cohorts based on the subject identifiers
  proportions_df <- proportions_df %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "SY") ~ "RR10 Flight",
      str_detect(Subject, "GY") ~ "Ground Young",
      str_detect(Subject, "FY") ~ "Return Young",
      str_detect(Subject, "GO") ~ "Ground Old",
      str_detect(Subject, "FO") ~ "Return Old"
    ))
  
  # Specify order for comparison
  proportions_df$Cohort <- factor(proportions_df$Cohort, levels = cohort_order)
  proportions_df$CellType <- factor(proportions_df$CellType, levels = celltype_order)
  
  # Plot with boxplot and points to show the data distribution
  ggplot(proportions_df, aes(x = Cohort, y = Fraction, fill = Cohort, color = Cohort)) +
    geom_point(size=1.5, position = position_dodge(width=0.75)) +
    geom_boxplot(alpha=0.75, outlier.shape = NA, position=position_dodge(width=0.75)) +
    facet_wrap(~CellType, scales = "free_x", ncol = 5) +  # Adjust the number of columns based on the number of cell types
    theme_minimal() + NoLegend() +
    labs(title = "Relative Abundance of Each Cell Type by Cohort",
         x = "Cohort",
         y = "Fraction of Cells") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, NA)) +
    scale_fill_manual(values = my_colors) +  # Apply custom colors
    scale_color_manual(values = my_colors)
}



#### Fraction of cells in each phase of the cell cycle for each cell type and each cohort

phase_by_celltype <- function(obj) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  
  for (cell_type in cell_types) {
    # Subset Seurat object by cell type
    sub_obj <- subset(obj, idents = cell_type)
    phases <- unique(sub_obj$Phase)
    ct_table <- table(sub_obj$subject, sub_obj$Phase)
    ct_proportions <- prop.table(ct_table, margin = 1)
    proportions_df <- as.data.frame(ct_proportions)
    names(proportions_df) <- c("Subject", "Phase", "Fraction")
    
    proportions_df <- proportions_df %>%
      mutate(Cohort = case_when(
          str_detect(Subject, "RR10GC") ~ "RR10 Ground",
          str_detect(Subject, "SY") ~ "RR10 Flight",
          str_detect(Subject, "GY") ~ "Ground Young",
          str_detect(Subject, "FY") ~ "Return Young",
          str_detect(Subject, "GO") ~ "Ground Old",
          str_detect(Subject, "FO") ~ "Return Old"
      ))
    proportions_df$Cohort <- factor(proportions_df$Cohort, levels = cohort_order)
    print(head(proportions_df, n=10))
    
    plot <- ggplot(proportions_df, aes(x = Cohort, y = Fraction, fill = Phase, color = Phase)) +
      geom_point(size=1.5, position = position_dodge(width=.75)) +
      geom_boxplot(alpha=0.5) +
      theme_minimal() +
      labs(title = paste0("Cell Phase Abundance in ",cell_type," by Cohort"),
           x = "Cohort",
           y = "Proportion of Phase") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(limits = c(0, NA))
    plots[[cell_type]] <- plot
  }
  plot_names <- sort(names(plots))
  sorted_plots <- plots[plot_names]
  
  # Combine all plots
  plot_grid(plotlist = sorted_plots, align = "v", ncol = 3)
}




#### Percent of cells in each cell type for each cohort that are expressing Gene

percent_celltype_expressing <- function(obj, thresh, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 Flight", "Ground Young", "Return Young", "Ground Old", "Return Old")
  
  # Compute percent expressing for the entire object
  expressiondf_all <- as.data.frame(Percent_Expressing(obj, layer="counts",assay="RNA", group_by = "subject", features = gene, threshold = thresh))
  long_data_all <- pivot_longer(expressiondf_all, cols = everything(), names_to = "Subject", values_to = "Percent")
  
  long_data_all <- long_data_all %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "SY") ~ "RR10 Flight",
      str_detect(Subject, "GY") ~ "Ground Young",
      str_detect(Subject, "FY") ~ "Return Young",
      str_detect(Subject, "GO") ~ "Ground Old",
      str_detect(Subject, "FO") ~ "Return Old"
    ))
  long_data_all$Cohort <- factor(long_data_all$Cohort, levels = cohort_order)
  
  plot_all <- ggplot(long_data_all, aes(x = Cohort, y = Percent, color = Cohort, fill = Cohort)) +
    geom_point(size=3.5, position = position_dodge(width=.75)) +
    geom_boxplot(alpha=0.5) +
    geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
    theme_minimal() + NoLegend() +
    labs(title = paste0("Percent of all Cells Expressing >",thresh," counts of ",gene), subtitle="by Cohort and Subject",
         y = paste("Percent Expressing ", gene)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
    scale_y_continuous(limits = c(0, NA))
  
  plots[["All Cell Types"]] <- plot_all
  
  # Loop over each cell type for specific plots
  for (cell_type in cell_types) {

    sub_obj <- subset(obj, idents = cell_type)
    
    expressiondf <- as.data.frame(Percent_Expressing(sub_obj, layer="counts",assay="RNA", group_by = "subject", features = gene, threshold = thresh))
    long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "Percent")
    
    long_data <- long_data %>%
      mutate(Cohort = case_when(
        str_detect(Subject, "RR10GC") ~ "RR10 Ground",
        str_detect(Subject, "SY") ~ "RR10 Flight",
        str_detect(Subject, "GY") ~ "Ground Young",
        str_detect(Subject, "FY") ~ "Return Young",
        str_detect(Subject, "GO") ~ "Ground Old",
        str_detect(Subject, "FO") ~ "Return Old"
      ))
    long_data$Cohort <- factor(long_data$Cohort, levels = cohort_order)
    
    plot <- ggplot(long_data, aes(x = Cohort, y = Percent, color = Cohort, fill = Cohort)) +
      geom_point(size=3.5, position = position_dodge(width=.75)) +
      geom_boxplot(alpha=0.5) +
      geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
      theme_minimal() + NoLegend() +
      labs(title = paste0("Percent of ", cell_type," Expressing >",thresh," counts of ",gene), subtitle="by Cohort and Subject",
           y = paste("Percent Expressing ", gene)) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
      scale_y_continuous(limits = c(0, NA))
    plots[[cell_type]] <- plot
  }
  plot_names <- sort(names(plots))
  sorted_plots <- plots[plot_names]
  
  # Combine all plots
  plot_grid(plotlist = sorted_plots, align = "v", ncol = 3)
}




#########################
######## Testing ########
#########################




##### Average expression of Gene across cohorts by subject

expression_by_spln_celltypes <- function(obj, layer_to_use, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("Young Ground", "Young Recovery", "Old Ground", "Old Recovery")
  
  # gather data for pooled dataset
  expressiondf_all <- as.data.frame(AverageExpression(obj, layer= layer_to_use, assay="RNA", group.by = "subject", features = c(gene)))
  long_data_all <- pivot_longer(expressiondf_all, cols = everything(), names_to = "Subject", values_to = "gene")
  
  # define cohort for the entire dataset
  long_data_all <- long_data_all %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.GY") ~ "Young Ground",
      str_detect(Subject, "\\.FY") ~ "Young Recovery",
      str_detect(Subject, "\\.GO") ~ "Old Ground",
      str_detect(Subject, "\\.FO") ~ "Old Recovery"
    ))
  long_data_all$Cohort <- factor(long_data_all$Cohort, levels = cohort_order)
  
  # generate plot for all cell types
  plot_all <- ggplot(long_data_all, aes(x = Cohort, y = gene, color = Cohort, fill = Cohort)) +
    geom_point(size=3.5, position = position_dodge(width=.75)) +
    geom_boxplot(alpha=0.5,) +
    geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
    theme_minimal() + NoLegend() +
    labs(title = paste(gene, "Expression Across All Cell Types"), subtitle="by Subject and Cohort",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
    scale_y_continuous(limits = c(0, NA))
  
  # add plot for all cell types to plots list
  plots[["All Cell Types"]] <- plot_all
  
  # loop over each cell type for specific plots
  for (cell_type in cell_types) {
    
    sub_obj <- subset(obj, idents = cell_type)
    
    # compute average expression for the subset
    expressiondf <- as.data.frame(AverageExpression(sub_obj, layer= layer_to_use, assay="RNA", group.by = "subject", features = c(gene)))
    long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
    
    # assign cohorts to subjects
    long_data <- long_data %>%
      mutate(Cohort = case_when(
        str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
        str_detect(Subject, "\\.SY") ~ "RR10 Flight",
        str_detect(Subject, "\\.GY") ~ "Ground Young",
        str_detect(Subject, "\\.FY") ~ "Return Young",
        str_detect(Subject, "\\.GO") ~ "Ground Old",
        str_detect(Subject, "\\.FO") ~ "Return Old"
      ))
    long_data$Cohort <- factor(long_data$Cohort, levels = cohort_order)
    
    # generate plot for current cell type
    plot <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort, fill = Cohort)) +
      geom_point(size=3.5, position = position_dodge(width=.75)) +
      geom_boxplot(alpha=0.5) +
      geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
      theme_minimal() + NoLegend() +
      labs(title = paste(gene, "Expression in", cell_type), subtitle="by Subject and Cohort",
           y = paste(gene, "Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
      scale_y_continuous(limits = c(0, NA))
    plots[[cell_type]] <- plot
  }
  plot_names <- sort(names(plots))
  sorted_plots <- plots[plot_names]
  
  # combine all plots
  plot_grid(plotlist = sorted_plots, align = "v", ncol = 3)
}




























expression_by_celltypes_shaded <- function(obj, gene, cohort) {
  # Subset the Seurat object for the specific cohort and calculate average expression
  sub_obj <- subset(obj, subset = Cohort == cohort)
  exp_data <- AverageExpression(sub_obj, features = gene, return.seurat = FALSE)
  
  # Prepare the data
  exp_data <- as.data.frame(exp_data$RNA)  # Change 'RNA' to the correct assay name if different
  colnames(exp_data)[1] <- "Subject"
  exp_data$GeneExpression <- exp_data[[gene]]
  exp_data <- exp_data[, c("Subject", "GeneExpression")]
  
  # Calculate the quartiles and IQR
  quartiles <- quantile(exp_data$GeneExpression, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- IQR(exp_data$GeneExpression, na.rm = TRUE)
  lower_bound <- quartiles[1] - 1.5 * iqr
  upper_bound <- quartiles[2] + 1.5 * iqr
  
  # Create the plot
  ggplot(exp_data, aes(x = Subject, y = GeneExpression)) +
    geom_point() +  # Scatter plot of expression
    geom_line(aes(group = 1), color = "blue", linewidth = 1) +  # Connect points with a line
    geom_ribbon(aes(ymin = lower_bound, ymax = upper_bound), fill = "blue", alpha = 0.2) +
    labs(title = paste("Expression of", gene, "in", cohort, "Cohort"),
         x = "Subject",
         y = "Gene Expression") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

































scatter_by_celltype <- function(obj, gene) {
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







scatter_by_celltype <- function(obj, gene) {
  # Compute average expression for each cell type
  expressiondf <- as.data.frame(AverageExpression(obj, layer= "counts", assay="SCT", group.by = c("ident", "subject"), features = c(gene)))
  long_data <- pivot_longer(expressiondf, cols = -c(ident), names_to = "Subject", values_to = "gene")
  
  # Define cohort
  long_data <- long_data %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "\\.SY") ~ "RR10 FLIGHT",
      str_detect(Subject, "\\.GY") ~ "RRRM2 Ground",
      str_detect(Subject, "\\.FY") ~ "RRRM2 Return",
      str_detect(Subject, "\\.GO") ~ "Ground Old",
      str_detect(Subject, "\\.FO") ~ "Flight Old"
    ))
  
  # List to store plots
  plots <- list()
  
  # Loop through each cell type
  for (cell_type in unique(long_data$ident)) {
    # Filter data for the current cell type
    type_data <- filter(long_data, ident == cell_type)
    
    # Generate plot for this cell type
    plot <- ggplot(type_data, aes(x = Cohort, y = gene, color = Cohort)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = Subject), hjust=-.3, size = 4) +
      theme_minimal() + 
      labs(title = paste(gene, "Expression in", cell_type, "by Subject and Cohort"),
           x = "Cohort",
           y = paste(gene, "Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(limits = c(0, NA))
    
    # Store plot in list
    plots[[cell_type]] <- plot
  }
  
  # Combine all plots
  plot_grid(plotlist = plots, align = "v", ncol = 1, labels = names(plots))
}











scatter_by_celltype <- function(obj, gene) {
  # Get unique cell types
  cell_types <- unique(Idents(obj))
  
  # Store plots in a list
  plots <- list()
  
  for (cell_type in cell_types) {
    # Subset Seurat object by cell type
    sub_obj <- subset(obj, idents = cell_type)
    
    # Compute average expression for the subset
    expressiondf <- as.data.frame(AverageExpression(sub_obj, layer= "counts", assay="RNA", group.by = "subject", features = c(gene)))
    long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
    
    # Define cohort
    long_data <- long_data %>%
      mutate(Cohort = case_when(
        str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
        str_detect(Subject, "\\.SY") ~ "RR10 FLIGHT",
        str_detect(Subject, "\\.GY") ~ "RRRM2 Ground",
        str_detect(Subject, "\\.FY") ~ "RRRM2 Return",
        str_detect(Subject, "\\.GO") ~ "Ground Old",
        str_detect(Subject, "\\.FO") ~ "Flight Old"
      ))
    
    # Generate plot for this cell type
    plot <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort)) +
      geom_point(size = 3) +
      geom_text_repel(aes(label = Subject), hjust=-.1, size = 3) +
      theme_minimal() + NoLegend() +
      labs(title = paste(gene, "Expression in", cell_type, "by Subject and Cohort"),
           x = cell_type,
           y = paste(gene, "Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      scale_y_continuous(limits = c(0, NA))
    
    # Store plot in list
    plots[[cell_type]] <- plot
  }
  
  # Combine all plots
  plot_grid(plotlist = plots, align = "v", ncol = 5)
}








#Boxed

scatter_by_celltype <- function(obj, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 FLIGHT", "RRRM2 Ground", "RRRM2 Return", "Ground Old", "Flight Old")
  
  for (cell_type in cell_types) {
    # Subset Seurat object by cell type
    sub_obj <- subset(obj, idents = cell_type)
    
    # Compute average expression for the subset
    expressiondf <- as.data.frame(AverageExpression(sub_obj, layer= "counts", assay="RNA", group.by = "subject", features = c(gene)))
    long_data <- pivot_longer(expressiondf, cols = everything(), names_to = "Subject", values_to = "gene")
    
    # Define cohort
    long_data <- long_data %>%
      mutate(Cohort = case_when(
        str_detect(Subject, "\\.RR10GC") ~ "RR10 Ground",
        str_detect(Subject, "\\.SY") ~ "RR10 FLIGHT",
        str_detect(Subject, "\\.GY") ~ "RRRM2 Ground",
        str_detect(Subject, "\\.FY") ~ "RRRM2 Return",
        str_detect(Subject, "\\.GO") ~ "Ground Old",
        str_detect(Subject, "\\.FO") ~ "Flight Old"
      ))
    long_data$Cohort <- factor(long_data$Cohort, levels = cohort_order)
    # Generate plot for this cell type
    plot <- ggplot(long_data, aes(x = Cohort, y = gene, color = Cohort)) +
      geom_point(size = 4) +
      geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3,max.overlaps = 12) +
      theme_minimal() + NoLegend()+
      labs(title = paste(gene, "Expression in", cell_type), subtitle="by Subject and Cohort",
           x = paste0(cell_type,"s by Cohort"),
           y = paste(gene, "Expression")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
      scale_y_continuous(limits = c(0, NA))
    plots[[cell_type]] <- plot
  }
  plot_grid(plotlist = plots, align = "v", ncol = 5)
}

VlnPlot(Femur, features=c("Mmut"), split.by = "Cohort", ncol=1,pt.size = 0,assay="RNA")


AverageExpression(Femur, features="Hbx",group.by = "Cohort")

AggregateExpression(Femur, features="Mmut",group.by = "Cohort", assays = "RNA")





#nice



#experiment with fraction of cells expressing gene mettric
scatter_pooled <- function(obj, gene) {
  cell_types <- unique(Idents(obj))
  plots <- list()
  cohort_order <- c("RR10 Ground", "RR10 FLIGHT", "RRRM2 Ground", "RRRM2 Return", "Ground Old", "Flight Old")
  
  # Compute average expression for the entire object
  expressiondf_all <- as.data.frame(Percent_Expressing(obj, layer="counts",assay="RNA", group_by = "subject", features = c(gene))*100)
  long_data_all <- pivot_longer(expressiondf_all, cols = everything(), names_to = "Subject", values_to = "Percent")
  
  # Define cohort for the entire dataset
  long_data_all <- long_data_all %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "SY") ~ "RR10 FLIGHT",
      str_detect(Subject, "GY") ~ "RRRM2 Ground",
      str_detect(Subject, "FY") ~ "RRRM2 Return",
      str_detect(Subject, "GO") ~ "Ground Old",
      str_detect(Subject, "FO") ~ "Flight Old"
    ))
  long_data_all$Cohort <- factor(long_data_all$Cohort, levels = cohort_order)
  
  # Generate plot for all cell types
  ggplot(long_data_all, aes(x = Cohort, y = Percent, color = Cohort)) +
    geom_point(size = 4) +
    geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 3, max.overlaps = 12) +
    theme_minimal() +
    labs(title = paste(gene, "Expression Across All Cell Types"), subtitle="by Subject and Cohort",
         x = "Cohorts",
         y = paste(gene, "Expression")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.background = element_rect(color = "black")) +
    scale_y_continuous(limits = c(0, NA))
  
}



















Femur$sample_number <- as.numeric(str_extract(Femur$subject, "\\d+$"))
Femur$sample_rank <- Femur$sample_number
#library(DescTools)
for(cohort in unique(Femur$Cohort)){
  Femur$sample_rank <- ifelse(Femur$Cohort == cohort, Rank(subset(Femur, subset = Cohort == cohort)$sample_number,ties.method = "dense"), Femur$sample_rank)}


Femur <- Femur %>% 
  group_by(Cohort) %>%
  mutate(rank_in_cohort = rank(sample_number)) 

ggplot(proportions_df, aes(x = Cohort, y = Fraction, color = Cohort)) +
  geom_point(size=3.5) +
  facet_wrap(~CellType, scales = "free_x",nrow = 1) +
  theme_minimal() +
  geom_text_repel(aes(label = str_sub(Subject, -2,-1)), hjust=-.2, size = 4,max.overlaps = 12) +
  labs(title = "Proportions of Cell Types by Subject and Cohort",
       x = "Cohort",
       y = "Proportion of Cell Types") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_brewer(palette = "Set1")










phase_proportions <- function(obj) {
  # Get unique cell types
  phases <- unique(marrow$Phase)
  
  # Create a table of counts of cell types per subject
  ct_table <- table(marrow$subject, marrow$Phase)
  
  # Compute proportions of cell types per subject
  ct_proportions <- prop.table(ct_table, margin = 1)
  
  # Convert to dataframe for ggplot
  proportions_df <- as.data.frame(ct_proportions)
  names(proportions_df) <- c("Subject", "Phase", "Fraction")
  
  # Define cohorts
  proportions_df <- proportions_df %>%
    mutate(Cohort = case_when(
      str_detect(Subject, "RR10GC") ~ "RR10 Ground",
      str_detect(Subject, "SY") ~ "RR10 FLIGHT",
      str_detect(Subject, "GY") ~ "RRRM2 Ground",
      str_detect(Subject, "FY") ~ "RRRM2 Return",
      str_detect(Subject, "GO") ~ "Ground Old",
      str_detect(Subject, "FO") ~ "Flight Old"
    ))
  
  # Plot proportions by cohort and cell type
  ggplot(proportions_df, aes(x = Cohort, y = Fraction, fill = Phase, color = Phase)) +
    geom_point(size=4) +
    geom_boxplot(alpha=0.5, position=position_dodge(width=0.75)) +
    theme_minimal() +
    labs(title = "Proportions of Cell Types by Subject and Cohort",
         x = "Cohort",
         y = "Proportion of Cell Types") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_y_continuous(limits = c(0, NA))
  #scale_fill_manual(values = c("red","red","red","red","orange","orange","orange","orange","blue", "blue","blue","blue","lightblue","lightblue","lightblue","lightblue"))
  
  
}




    

















as.data.frame(AverageExpression(Femur, layer= "counts", assay="RNA", group.by = "subject", features = c("Mthfr")))
as.data.frame(Percent_Expressing(Femur, layer="count",assay="RNA", group_by = "subject", features = "Mthfr"))



