##############################
# Simple 3D UMAP (Cube grid)
##############################

library(plotly)
# Load your seurat object and specify dims
# yourseuratobject <- Spleen
# n_dims <- 8

make_3d_umap <- function(yourseuratobject,n_dims){
  # Generate df
  yourseuratobject <- RunUMAP(yourseuratobject, dims = 1:n_dims, n.components = 3L)
  Embeddings(object = yourseuratobject, reduction = "umap")
  plot.data <- FetchData(object = yourseuratobject, vars = c("umap_1", "umap_2", "umap_3", "ident"))
  plot.data$label <- paste(rownames(plot.data))
  
  # Plot
  fig <- plot_ly(data = plot.data, 
          x = ~umap_1, y = ~umap_2, z = ~umap_3, 
          color = ~ident, 
          colors = c("lightseagreen",
                     "gray50",
                     "darkgreen",
                     "red4",
                     "red",
                     "turquoise4",
                     "black",
                     "yellow4",
                     "royalblue1",
                     "lightcyan3",
                     "peachpuff3",
                     "khaki3",
                     "gray20",
                     "orange2",
                     "royalblue4",
                     "yellow3",
                     "gray80",
                     "darkorchid1",
                     "lawngreen",
                     "plum2",
                     "darkmagenta"),
          type = "scatter3d", 
          mode = "markers", 
          opacity = .25,
  	  marker = list(size = 2, width=2),
          text=~label, 
          hoverinfo="text")
  
  axx <- list(
    nticks = 4,
    range = c(-10,10))
  
  axy <- list(
    nticks = 4,
    range = c(-10,10))
  
  axz <- list(
    nticks = 4,
    range = c(-10,10))
  
  
  fig <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
  fig}

fig_cube <- fig %>% layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz, aspectmode='cube')) # To maintain cubic aspect

#fig_cube
