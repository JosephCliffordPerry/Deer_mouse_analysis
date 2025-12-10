library(umap)
library(ggplot2)
data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Combined_relabelled_deer_mice_stats.txt", sep = "\t", header = TRUE)
numdata <- dplyr::select_if(data, is.numeric)
print("umap1 started")
umapdata <- numdata[, !grepl("seg|outline", names(numdata), ignore.case = TRUE)]
umapo <- umap(umapdata, preserve.seed = TRUE)
#save(umapo,file = "D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/umapo.rda")




rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt", sep = "\t", header = TRUE)

Phenotype_data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt", sep = "\t", header = TRUE)
unified_format_image<-gsub("24_11_06", "slide", rawdata$Image)
Sample_number<- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)
data_with_sample<-cbind(rawdata,Sample_number)
data_with_hypoxia<-merge(data_with_sample,Phenotype_data,by = "Sample_number")


graphdata<-as.data.frame(cbind(umapo$layout,data_with_hypoxia$Hypoxic_exposure))


ggplot(data = graphdata)+
  geom_point(aes(x = V1,y = V2,colour = factor(V3)))+
  theme(axis.text.x = element_blank(),axis.text.y = element_blank()) +
  labs(title = "General UMAP coloured by hypoxic exposure", colour = NULL)

graphdata

graphdata<-as.data.frame(cbind(umapo$layout,data_with_hypoxia$))


library(dbscan)
kNNdistplot(umapo$layout,minPts = 2)

umap_cluster<-dbscan(umapo$layout,eps = 0.05)
export<-cbind(rawdata$CellID,umap_cluster$cluster)
write.table(export,file = "D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/umap_clusters.txt",sep = "\t",row.names = FALSE, col.names = TRUE,quote = FALSE)
