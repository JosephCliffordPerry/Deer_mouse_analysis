library(NMAcompanion)
rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Combined_relabelled_deer_mice_stats.txt", sep = "\t", header = TRUE)
ROI<-get_regions_of_interest(rawdata)
#install.packages("dbscan")
library(dbscan)
#big_clusterer
set.seed(234)
clusters<-list()
clusters2<-list()
error_tagged_angle_dataset <- Extreme_angle_detector(data = rawdata)
data <- filter(error_tagged_angle_dataset, suspected_detection_error == 1)
angle_data <- data %>% dplyr::select(starts_with("Angle_profile_"))
diameter_data <- data %>% dplyr::select(starts_with("Diameter_profile_"))
radius_data <- data %>% dplyr::select(starts_with("Radius_profile_"))
angle_outliers <- make_outlier_cluster(angle_data, "angle")
diameter_outliers <- make_outlier_cluster(diameter_data, "diameter")
radius_outliers <- make_outlier_cluster(radius_data, "radius")
continumouse<-list(angle_outliers=angle_outliers,radius_outliers=radius_outliers,diameter_outliers=diameter_outliers)


kNNdistplot(ROI$`Angle_profile26:32`,minPts = 2)
clusters$'Angle_profile26:32'<-dbscan(ROI$`Angle_profile26:32`,eps = 10)

#clusters$'Angle_profile27:32'<- cbind(continumouse[["angle_outliers"]][["angle outliers"]],clusters$'Angle_profile27:32'$cluster)

clusters2[["Angle_profile26:32"]] <- list(
    `Angle_profile26:32` =continumouse[["angle_outliers"]][["angle outliers"]],
    Clustering_file = clusters[["Angle_profile26:32"]][["cluster"]]
  )


kNNdistplot(ROI$`Angle_profile43:48`,minPts = 2)
clusters$'Angle_profile43:48'<-dbscan(ROI$`Angle_profile43:48`,eps = 10)
#clusters$'Angle_profile43:48'<- cbind(cellID = continumouse[["angle_outliers"]][["angle outliers"]],Clustering_file = clusters$'Angle_profile43:48'$cluster)
clusters2[["Angle_profile43:48"]] <- list(
  `Angle_profile43:48` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Angle_profile43:48"]][["cluster"]]
)
kNNdistplot(ROI$Dataset_Length_of_hook_microns,minPts = 2)
clusters$Dataset_Length_of_hook_microns<-dbscan(ROI$Dataset_Length_of_hook_microns,eps = 0.1)
#clusters$Dataset_Length_of_hook_microns<- cbind(continumouse[["angle_outliers"]][["angle outliers"]],clusters$Dataset_Length_of_hook_microns$cluster)
clusters2[["Dataset_Length_of_hook_microns"]] <- list(
  `Dataset_Length_of_hook_microns` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Dataset_Length_of_hook_microns"]][["cluster"]]
)
kNNdistplot(ROI$Dataset_Width_of_body_microns,minPts = 2)
clusters$Dataset_Width_of_body_microns<-dbscan(ROI$Dataset_Width_of_body_microns,eps = 0.5)
#clusters$Dataset_Width_of_body_microns<- cbind(continumouse[["angle_outliers"]][["angle outliers"]],clusters$Dataset_Width_of_body_microns$cluster)
clusters2[["Dataset_Width_of_body_microns"]] <- list(
  `Dataset_Width_of_body_microns` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Dataset_Width_of_body_microns"]][["cluster"]]
)
clusters2<-append(clusters2,continumouse)
clusties <- data$CellID
   for (i in names(clusters2)) {
     x <- as.numeric(clusters2[[i]][["Clustering_file"]])
      # Assign a unique name to x based on the cluster index (i)
        col_name <- paste0("Clustering_", i)
        # Create a new data frame with the same number of rows as data
          x_df <- as.data.frame(x = x)

            # Rename the column to col_name
             colnames(x_df) <- col_name
             # Combine x into clusties with the assigned column name
               clusties <- cbind(clusties, x_df)
             }
  ids<-paste0(clusties$`Clustering_Angle_profile26:32`,clusties$`Clustering_Angle_profile43:48`,clusties$Clustering_Dataset_Length_of_hook_microns,clusties$Clustering_Dataset_Width_of_body_microns,clusties$Clustering_angle_outliers,clusties$Clustering_radius_outliers,clusties$Clustering_diameter_outliers)
cluster_characterising <- function(data = data, ids = ids) {
  id_dataframe <- Make_cluster_id_df(data = data, ids = ids)
  seeds <- as.vector(Make_hamming_seeds(id_dataframe, data))
  output <- Do_hamming_amalgamation(df = id_dataframe, CharVectors = seeds)
  return(output)
}

hamming_cluster<-cluster_characterising(data = data, ids = ids)
hamming_cluster$Clustering_file <- hamming_cluster$Cluster_characterising_ids
outlinedata <- data %>% dplyr::select(starts_with("Outline_Oriented"))
hamming_consensus <- Make_hamming_consensus_images(outlinedata = outlinedata, cluster = hamming_cluster)
############
#make first NMA export
export<-cbind(hamming_cluster$UUID,hamming_cluster$Cluster_characterising_ids)
write.table(export,file = "D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/clustering1.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)

#################

#Secondary analysis round




rawdata2<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Group_1_Cluster_6_stats.txt", sep = "\t", header = TRUE)
ROI<-get_regions_of_interest(rawdata2)
#install.packages("dbscan")
library(dbscan)
#big_clusterer
set.seed(234)
clusters<-list()
clusters2<-list()
error_tagged_angle_dataset <- Extreme_angle_detector(data = rawdata2)
data <- filter(error_tagged_angle_dataset, suspected_detection_error == 1)
angle_data <- data %>% dplyr::select(starts_with("Angle_profile_"))
diameter_data <- data %>% dplyr::select(starts_with("Diameter_profile_"))
radius_data <- data %>% dplyr::select(starts_with("Radius_profile_"))
angle_outliers <- make_outlier_cluster(angle_data, "angle")
diameter_outliers <- make_outlier_cluster(diameter_data, "diameter")
radius_outliers <- make_outlier_cluster(radius_data, "radius")
continumouse<-list(angle_outliers=angle_outliers,radius_outliers=radius_outliers,diameter_outliers=diameter_outliers)


kNNdistplot(ROI$`Angle_profile27:32`,minPts = 2)
clusters$'Angle_profile26:32'<-dbscan(ROI$`Angle_profile27:32`,eps = 10)

#clusters$'Angle_profile27:32'<- cbind(continumouse[["angle_outliers"]][["angle outliers"]],clusters$'Angle_profile27:32'$cluster)

clusters2[["Angle_profile26:32"]] <- list(
  `Angle_profile26:32` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Angle_profile26:32"]][["cluster"]]
)


kNNdistplot(ROI$`Angle_profile43:48`,minPts = 2)
clusters$'Angle_profile43:48'<-dbscan(ROI$`Angle_profile43:48`,eps = 10)
#clusters$'Angle_profile43:48'<- cbind(cellID = continumouse[["angle_outliers"]][["angle outliers"]],Clustering_file = clusters$'Angle_profile43:48'$cluster)
clusters2[["Angle_profile43:48"]] <- list(
  `Angle_profile43:48` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Angle_profile43:48"]][["cluster"]]
)

clusters2<-append(clusters2,continumouse)
clusties <- data$CellID
for (i in names(clusters2)) {
  x <- as.numeric(clusters2[[i]][["Clustering_file"]])
  # Assign a unique name to x based on the cluster index (i)
  col_name <- paste0("Clustering_", i)
  # Create a new data frame with the same number of rows as data
  x_df <- as.data.frame(x = x)

  # Rename the column to col_name
  colnames(x_df) <- col_name
  # Combine x into clusties with the assigned column name
  clusties <- cbind(clusties, x_df)
}
ids<-paste0(clusties$`Clustering_Angle_profile27:32`,clusties$`Clustering_Angle_profile43:48`,clusties$Clustering_angle_outliers,clusties$Clustering_radius_outliers,clusties$Clustering_diameter_outliers)
cluster_characterising <- function(data = data, ids = ids) {
  id_dataframe <- Make_cluster_id_df(data = data, ids = ids)
  seeds <- as.vector(Make_hamming_seeds(id_dataframe, data))
  output <- Do_hamming_amalgamation(df = id_dataframe, CharVectors = seeds)
  return(output)
}

hamming_cluster<-cluster_characterising(data = data, ids = ids)
hamming_cluster$Clustering_file <- hamming_cluster$Cluster_characterising_ids
outlinedata <- data %>% dplyr::select(starts_with("Outline_Oriented"))
hamming_consensus <- Make_hamming_consensus_images(outlinedata = outlinedata, cluster = hamming_cluster)

#make second NMA export
export<-cbind(hamming_cluster$UUID,hamming_cluster$Cluster_characterising_ids)
write.table(export,file = "D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/clustering2.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)

#####################################


#clustering 3

rawdata3<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Group_1_Cluster_4_stats.txt", sep = "\t", header = TRUE)
ROI<-get_regions_of_interest(rawdata3)
#install.packages("dbscan")
library(dbscan)
#big_clusterer
set.seed(234)
clusters<-list()
clusters2<-list()
error_tagged_angle_dataset <- Extreme_angle_detector(data = rawdata3)
data <- filter(error_tagged_angle_dataset, suspected_detection_error == 1)
angle_data <- data %>% dplyr::select(starts_with("Angle_profile_"))
diameter_data <- data %>% dplyr::select(starts_with("Diameter_profile_"))
radius_data <- data %>% dplyr::select(starts_with("Radius_profile_"))
angle_outliers <- make_outlier_cluster(angle_data, "angle")
diameter_outliers <- make_outlier_cluster(diameter_data, "diameter")
radius_outliers <- make_outlier_cluster(radius_data, "radius")
continumouse<-list(angle_outliers=angle_outliers,radius_outliers=radius_outliers,diameter_outliers=diameter_outliers)



kNNdistplot(ROI$`Angle_profile43:47`,minPts = 2)
clusters$'Angle_profile43:47'<-dbscan(ROI$`Angle_profile43:47`,eps = 5)
#clusters$'Angle_profile43:48'<- cbind(cellID = continumouse[["angle_outliers"]][["angle outliers"]],Clustering_file = clusters$'Angle_profile43:48'$cluster)
clusters2[["Angle_profile43:47"]] <- list(
  `Angle_profile43:47` =continumouse[["angle_outliers"]][["angle outliers"]],
  Clustering_file = clusters[["Angle_profile43:47"]][["cluster"]]
)

clusters2<-append(clusters2,continumouse)
clusties <- data$CellID
for (i in names(clusters2)) {
  x <- as.numeric(clusters2[[i]][["Clustering_file"]])
  # Assign a unique name to x based on the cluster index (i)
  col_name <- paste0("Clustering_", i)
  # Create a new data frame with the same number of rows as data
  x_df <- as.data.frame(x = x)

  # Rename the column to col_name
  colnames(x_df) <- col_name
  # Combine x into clusties with the assigned column name
  clusties <- cbind(clusties, x_df)
}
ids<-paste0(clusties$`Clustering_Angle_profile43:47`,clusties$Clustering_angle_outliers,clusties$Clustering_radius_outliers,clusties$Clustering_diameter_outliers)
cluster_characterising <- function(data = data, ids = ids) {
  id_dataframe <- Make_cluster_id_df(data = data, ids = ids)
  seeds <- as.vector(Make_hamming_seeds(id_dataframe, data))
  output <- Do_hamming_amalgamation(df = id_dataframe, CharVectors = seeds)
  return(output)
}

hamming_cluster<-cluster_characterising(data = data, ids = ids)
hamming_cluster$Clustering_file <- hamming_cluster$Cluster_characterising_ids
outlinedata <- data %>% dplyr::select(starts_with("Outline_Oriented"))
hamming_consensus <- Make_hamming_consensus_images(outlinedata = outlinedata, cluster = hamming_cluster)

#make second NMA export
export<-cbind(hamming_cluster$UUID,hamming_cluster$Cluster_characterising_ids)
write.table(export,file = "D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/clustering3.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
hamming_consensus

