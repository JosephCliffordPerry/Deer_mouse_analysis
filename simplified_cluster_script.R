
rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt", sep = "\t", header = TRUE)

Phenotype_data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt", sep = "\t", header = TRUE)
unified_format_image<-gsub("24_11_06", "slide", rawdata$Image)
Sample_number<- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)
data_with_sample<-cbind(rawdata,Sample_number)
data_with_phenotypes<-merge(data_with_sample,Phenotype_data,by = "Sample_number")

library(NMAcompanion)
ROIs<-NMAcompanion::get_regions_of_interest(data_with_phenotypes)
library(dbscan)
large_clustering<-list()
library(dbscan)

# Perform DBSCAN clustering on selected ROI features
large_clustering <- list()

large_clustering$`Angle_profile26:32` <- dbscan(ROIs$`Angle_profile26:32`, eps = 10)
large_clustering$`Angle_profile43:48` <- dbscan(ROIs$`Angle_profile43:48`, eps = 10)
large_clustering$Dataset_Length_of_hook_microns <- dbscan(ROIs$Dataset_Length_of_hook_microns, eps = 0.1)
large_clustering$Dataset_Width_of_body_microns <- dbscan(ROIs$Dataset_Width_of_body_microns, eps = 0.5)

# Combine clustering results back into ROIs
clusters <- list()

clusters$`Angle_profile26:32` <- cbind(ROIs$`Angle_profile26:32`,
                                       Clustering_file = large_clustering$`Angle_profile26:32`$cluster)

clusters$`Angle_profile43:48` <- cbind(ROIs$`Angle_profile43:48`,
                                       Clustering_file = large_clustering$`Angle_profile43:48`$cluster)

clusters$Dataset_Length_of_hook_microns <- cbind(ROIs$Dataset_Length_of_hook_microns,
                                                 Clustering_file = large_clustering$Dataset_Length_of_hook_microns$cluster)

clusters$Dataset_Width_of_body_microns <- cbind(ROIs$Dataset_Width_of_body_microns,
                                                Clustering_file = large_clustering$Dataset_Width_of_body_microns$cluster)
###################################
cont_clusters<-find_contious_clusters(data_with_sample)
clusters2<-append(clusters,cont_clusters)
library(ggplot2)
## test ROI vs umap
make_umap_plots <- function(roi_data, phenotype_data, cluster_data, roi_name, save_path = "umap_plots") {
  # Create output directory if it doesn’t exist
  if (!dir.exists(save_path)) dir.create(save_path)

  # Run UMAP
  umap_result <- UmapNMAdata(data = roi_data)

  # Create a list of ggplots
  plots <- list(
    Hypoxic_exposure = umap_result +
      aes(colour = phenotype_data$Hypoxic_exposure) +
      ggtitle(paste(roi_name, "- Hypoxic Exposure")),

    Population = umap_result +
      aes(colour = phenotype_data$Population) +
      ggtitle(paste(roi_name, "- Population")),

    Family = umap_result +
      aes(colour = phenotype_data$Family) +
      ggtitle(paste(roi_name, "- Family")),

    Clustering_file = umap_result +
      aes(colour = cluster_data) +
      ggtitle(paste(roi_name, "- Clustering"))
  )

  # Save each plot as a PNG
  for (plot_name in names(plots)) {
    file_name <- file.path(save_path, paste0(roi_name, "_", plot_name, ".png"))
    ggsave(filename = file_name, plot = plots[[plot_name]], width = 6, height = 5, dpi = 300)
  }

  return(plots)
}
#
# Angle_profile4348plots<-make_umap_plots(roi_data = ROIs$`Angle_profile43:48`,
#                                         phenotype_data = data_with_phenotypes,
#                                         cluster_data = clusters2$`Angle_profile43:48`$Clustering_file,
#                                         roi_name =  "Angle_profile43:48")
# Dataset_Length_of_hook_micronsplots<-make_umap_plots(roi_data = ROIs$Dataset_Length_of_hook_microns,
#                                         phenotype_data = data_with_phenotypes,
#                                         cluster_data = clusters2$Dataset_Length_of_hook_microns$Clustering_file,
#                                         roi_name =  "Dataset_Length_of_hook_microns")
# Angle_profile4348plots<-make_umap_plots(roi_data = ROIs$`Angle_profile43:48`,
#                                         phenotype_data = data_with_phenotypes,
#                                         cluster_data = clusters2$`Angle_profile43:48`$Clustering_file,
#                                         roi_name =  "Angle_profile43:48")
#
# amalgamated_features<-NMAcompanion::Amalgamate_morphological_features(data_with_sample,clusters2)
# List of ROIs you want to plot
roi_names <- c("Angle_profile26:32",
               "Angle_profile43:48",
               "Dataset_Length_of_hook_microns",
               "Dataset_Width_of_body_microns")

# Create a named list to store all plot sets
all_umap_plots <- list()

# Loop over each ROI
for (roi_name in roi_names) {

  message("Generating UMAP plots for: ", roi_name)

  all_umap_plots[[roi_name]] <- make_umap_plots(
    roi_data = ROIs[[roi_name]],
    phenotype_data = data_with_phenotypes,
    cluster_data = clusters2[[roi_name]]$Clustering_file,
    roi_name = roi_name,
    save_path = file.path("umap_plots", roi_name)  # Subfolder per ROI
  )
}






MakeNMAConsensus <- function(data, groups = NULL, highlight = NULL) {
  if (!hasName(data, "Outline_OrientedCoordinates_X_0")) {
    stop("Missing outlines check NMA export, consensus images will be unavailable.")
  }

  Outliney <- t(data %>% select(starts_with("Outline_OrientedCoordinates_Y")))
  Outlinex <- t(data %>% select(starts_with("Outline_OrientedCoordinates_X")))

  if (!is.null(groups)) {
    if (length(groups) != ncol(Outliney)) {
      stop("Length of groups must match the number of outlines.")
    }

    # Calculate mean outline per group
    plot_data <- purrr::map_dfr(unique(groups), function(g) {
      group_idx <- which(groups == g)
      group_y <- rowMeans(Outliney[, group_idx, drop = FALSE])
      group_x <- rowMeans(Outlinex[, group_idx, drop = FALSE])

      data.frame(
        x = group_x,
        y = group_y,
        group = g,
        point = seq_along(group_x)
      )
    })

    consensus <- ggplot(plot_data, aes(x = x, y = y, fill = group, group = group)) +
      geom_polygon(alpha = 0.5, color = "black", linewidth = 1) +
      facet_wrap(~ group) +
      theme(
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()
      ) +
      labs(title = "Stitched Consensus Polygons")

    # Add highlighted points if provided
    if (!is.null(highlight)) {
      consensus <- consensus +
        geom_point(
          data = subset(plot_data, point %in% highlight),
          aes(x = x, y = y),
          color = "red",
          size = 2
        )
    }

  } else {
    # Overall mean outline
    meanOutlineY <- rowMeans(Outliney)
    meanOutlineX <- rowMeans(Outlinex)

    polygon_df <- data.frame(
      x = meanOutlineX,
      y = meanOutlineY,
      point = seq_along(meanOutlineX)
    )

    consensus <- ggplot(polygon_df, aes(x = x, y = y)) +
      geom_polygon(alpha = 0.7) +
      theme_minimal()

    # Highlight points if provided
    if (!is.null(highlight)) {
      consensus <- consensus +
        geom_point(
          data = subset(polygon_df, point %in% highlight),
          aes(x = x, y = y),
          color = "red",
          size = 2
        )
    }
  }

  return(consensus)
}
