#potential sub-pipeline on varied cells
library(NMAcompanion)

rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt", sep = "\t", header = TRUE)

Phenotype_data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt", sep = "\t", header = TRUE)
unified_format_image<-gsub("24_11_06", "slide", rawdata$Image)
Sample_number<- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)
data_with_sample<-cbind(rawdata,Sample_number)
data_with_phenotypes<-merge(data_with_sample,Phenotype_data,by = "Sample_number")

cont_clusters<-find_contious_clusters(data_with_phenotypes)

# Sum all "Clustering_file" columns across dataframes
sum_clustering <- Reduce(`+`, lapply(cont_clusters, function(df) {
  # Ensure numeric: converts factors/characters safely
  as.numeric(as.character(df$Clustering_file))
}))
# Create a combined dataframe
tomerge <- data.frame(sum_clustering = sum_clustering,CellID = data_with_phenotypes$CellID)
combined_df <- merge(data_with_phenotypes,tomerge,by = "CellID")
# ---- Step 2: Filter based on thresholds ----

# <= 3
Standard_morph <- subset(combined_df, sum_clustering <= 3)

# 4 - 5
Small <- subset(combined_df, sum_clustering >= 4 & sum_clustering <= 5)

# = 6
Outlier <- subset(combined_df, sum_clustering == 6)

# ---- Step 3: (Optional) Print summary ----
cat("<=3 rows:", nrow(Standard_morph), "\n")
cat("4-5 rows:", nrow(Small), "\n")
cat("=6 rows:", nrow(Outlier), "\n")


#mini analysis of outliers
ROIs<-NMAcompanion::get_regions_of_interest(Outlier)
Clusters<- NMAcompanion::Cluster_ROI_list(ROIs)
morphological_features<-NMAcompanion::Amalgamate_morphological_features(Outlier,Clusters)

graph<-Graph_clustered_ROIs(Clusters,Outlier)

#mini analysis of Small
ROIs<-NMAcompanion::get_regions_of_interest(Small)
ROIs<-ROIs[1:11]
Clusters<- NMAcompanion::Cluster_ROI_list(ROIs)
morphological_features<-NMAcompanion::Amalgamate_morphological_features(Small,Clusters)

graph<-Graph_clustered_ROIs(Clusters,Small)

########################################
library(dplyr)
library(ggplot2)
library(viridis)

# --- Step 1: Rename column for consistency ---
names(data_with_phenotypes)[names(data_with_phenotypes) == "CellID"] <- "UUID"

# --- Step 2: Initialize list to store plots ---
plot_list <- list()

# --- Step 3: Loop through all morphological features ---
for (i in seq_along(morphological_features)) {

  # Merge phenotype and clustering data
  graphdata <- merge(
    data_with_phenotypes,
    morphological_features[[i]][["Cluster_dataframe"]],
    by = "UUID"
  )

  # Summarize number of cells per grouping and cluster
  graphdata_summary <- graphdata %>%
    group_by(Hypoxic_exposure, Population, Family, Clustering_file) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    # Calculate percentages per group (Hypoxia, Population, Family)
    group_by(Hypoxic_exposure) %>%
    mutate(percent_hypoxia = 100 * cell_count / sum(cell_count)) %>%
    ungroup() %>%
    group_by(Population) %>%
    mutate(percent_population = 100 * cell_count / sum(cell_count)) %>%
    ungroup() %>%
    group_by(Family) %>%
    mutate(percent_family = 100 * cell_count / sum(cell_count)) %>%
    ungroup()

  # --- Hypoxia plot ---
  p_hypoxia <- ggplot(graphdata_summary, aes(
    x = Hypoxic_exposure, y = percent_hypoxia, fill = as.factor(Clustering_file)
  )) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = "Hypoxia Level",
      y = "Percentage of Cells",
      title = paste("Feature", i, "- Hypoxia")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          panel.grid.major.x = element_blank())

  # --- Population plot ---
  p_population <- ggplot(graphdata_summary, aes(
    x = Population, y = percent_population, fill = as.factor(Clustering_file)
  )) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = "Population",
      y = "Percentage of Cells",
      title = paste("Feature", i, "- Population")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          panel.grid.major.x = element_blank())

  # --- Family plot ---
  p_family <- ggplot(graphdata_summary, aes(
    x = Family, y = percent_family, fill = as.factor(Clustering_file)
  )) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      x = "Family",
      y = "Percentage of Cells",
      title = paste("Feature", i, "- Family")
    ) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "top",
          panel.grid.major.x = element_blank())

  # Store plots in the list
  plot_list[[paste0("Feature_", i, "_Hypoxia")]] <- p_hypoxia
  plot_list[[paste0("Feature_", i, "_Population")]] <- p_population
  plot_list[[paste0("Feature_", i, "_Family")]] <- p_family
}




library(dplyr)
library(ggplot2)
library(viridis)

 # Summarize counts per cluster × Family × Hypoxia
graphdata_summary <- graphdata %>%
  group_by(Family, Hypoxic_exposure, Clustering_file) %>%
  summarise(cell_count = n(), .groups = "drop")

# Calculate percent within each cluster
graphdata_summary <- graphdata_summary %>%
  group_by(Clustering_file) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  ungroup() %>%
  # Create combined factor for X-axis
  mutate(Family_Hypoxia = paste(Family, Hypoxic_exposure, sep = "_"))

# Plot
ggplot(graphdata_summary, aes(
  x = Family_Hypoxia,
  y = percent,
  fill = as.factor(Clustering_file)
)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    x = "Family & Hypoxia",
    y = "Percentage of Cluster",
    title = "Distribution of Clusters Across Family & Hypoxia Levels"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )




