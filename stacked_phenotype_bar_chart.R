library(dplyr)
library(ggplot2)

rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt", sep = "\t", header = TRUE)

Phenotype_data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt", sep = "\t", header = TRUE)
rawdata$Image[1]
unified_format_image<-gsub("24_11_06", "slide", rawdata$Image)
Sample_number<- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)
data_with_sample<-cbind(rawdata,Sample_number)
data<-merge(data_with_sample,Phenotype_data,by = "Sample_number")

data$Dataset<-gsub("Group_1_Cluster_8","small_cells",data$Dataset)
data$Dataset<-gsub("Group_1_Cluster_7","normal_cells",data$Dataset)
data$Dataset <- ifelse(data$Dataset %in% c("small_cells", "normal_cells"),
                            data$Dataset,
                            "other_cells")

graphdata<-as.data.frame(cbind(data$Sample_number,data$Hypoxic_exposure,data$Dataset))





sample_ratios<-as.data.frame(table(graphdata$V1,graphdata$V3))
hypoxic_samples<-as.data.frame(table(graphdata$V1,graphdata$V2))

graphdata2<-merge(hypoxic_samples,sample_ratios,by = "Var1")

graphdata3 <- graphdata2[graphdata2$Freq.x != 0, ]
graphdata3$percent<-(graphdata3$Freq.y/graphdata3$Freq.x)*100


sorted_graphdata <- graphdata3 %>%
  group_by(Var1) %>%
  mutate(normal_cell_percent = ifelse(Var2.y == "normal_cells", percent, 0)) %>%
  summarise(total_normal_cell_percent = sum(normal_cell_percent, na.rm = TRUE)) %>%
  right_join(sorted_graphdata, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Order descending by "normal_cells"


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by hypoxia"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend

###############
#populations

graphdata<-as.data.frame(cbind(data$Sample_number,data$Population,data$Dataset))





sample_ratios<-as.data.frame(table(graphdata$V1,graphdata$V3))
Population_samples<-as.data.frame(table(graphdata$V1,graphdata$V2))

graphdata2<-merge(Population_samples,sample_ratios,by = "Var1")

graphdata3 <- graphdata2[graphdata2$Freq.x != 0, ]
graphdata3$percent<-(graphdata3$Freq.y/graphdata3$Freq.x)*100


# Step 1: Calculate the total "normal_cells" percentage for each Var1
normal_cell_order <- graphdata3 %>%
  group_by(Var1) %>%
  summarise(total_normal_cell_percent = sum(ifelse(Var2.y == "normal_cells", percent, 0), na.rm = TRUE))

# Step 2: Join this summary back with the original dataset
sorted_graphdata <- graphdata3 %>%
  left_join(normal_cell_order, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Reorder Var1 based on total_normal_cell_percent


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by Population"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend


#populations

graphdata<-as.data.frame(cbind(data$Sample_number,data$Population,data$Dataset))





sample_ratios<-as.data.frame(table(graphdata$V1,graphdata$V3))
Population_samples<-as.data.frame(table(graphdata$V1,graphdata$V2))

graphdata2<-merge(Population_samples,sample_ratios,by = "Var1")

graphdata3 <- graphdata2[graphdata2$Freq.x != 0, ]
graphdata3$percent<-(graphdata3$Freq.y/graphdata3$Freq.x)*100


# Step 1: Calculate the total "normal_cells" percentage for each Var1
normal_cell_order <- graphdata3 %>%
  group_by(Var1) %>%
  summarise(total_normal_cell_percent = sum(ifelse(Var2.y == "normal_cells", percent, 0), na.rm = TRUE))

# Step 2: Join this summary back with the original dataset
sorted_graphdata <- graphdata3 %>%
  left_join(normal_cell_order, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Reorder Var1 based on total_normal_cell_percent


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by Population"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend




#populations

graphdata<-as.data.frame(cbind(data$Sample_number,data$Population,data$Dataset))





sample_ratios<-as.data.frame(table(graphdata$V1,graphdata$V3))
Population_samples<-as.data.frame(table(graphdata$V1,graphdata$V2))

graphdata2<-merge(Population_samples,sample_ratios,by = "Var1")

graphdata3 <- graphdata2[graphdata2$Freq.x != 0, ]
graphdata3$percent<-(graphdata3$Freq.y/graphdata3$Freq.x)*100


# Step 1: Calculate the total "normal_cells" percentage for each Var1
normal_cell_order <- graphdata3 %>%
  group_by(Var1) %>%
  summarise(total_normal_cell_percent = sum(ifelse(Var2.y == "normal_cells", percent, 0), na.rm = TRUE))

# Step 2: Join this summary back with the original dataset
sorted_graphdata <- graphdata3 %>%
  left_join(normal_cell_order, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Reorder Var1 based on total_normal_cell_percent


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by Population"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend


####################

graphdata<-as.data.frame(cbind(data$Sample_number,data$Family,data$Dataset))

sample_ratios<-as.data.frame(table(graphdata$V1,graphdata$V3))
Family_samples<-as.data.frame(table(graphdata$V1,graphdata$V2))

graphdata2<-merge(Family_samples,sample_ratios,by = "Var1")

graphdata3 <- graphdata2[graphdata2$Freq.x != 0, ]
graphdata3$percent<-(graphdata3$Freq.y/graphdata3$Freq.x)*100


# Step 1: Calculate the total "normal_cells" percentage for each Var1
normal_cell_order <- graphdata3 %>%
  group_by(Var1) %>%
  summarise(total_normal_cell_percent = sum(ifelse(Var2.y == "normal_cells", percent, 0), na.rm = TRUE))

# Step 2: Join this summary back with the original dataset
sorted_graphdata <- graphdata3 %>%
  left_join(normal_cell_order, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Reorder Var1 based on total_normal_cell_percent


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "identity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by Population"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend

# Step 1: Calculate the total "normal_cells" percentage for each Var1
normal_cell_order <- graphdata3 %>%
  group_by(Var1) %>%
  summarise(total_normal_cell_percent = sum(ifelse(Var2.y == "normal_cells", percent, 0), na.rm = TRUE))

# Step 2: Join this summary back with the original dataset
sorted_graphdata <- graphdata3 %>%
  left_join(normal_cell_order, by = "Var1") %>%
  mutate(Var1 = reorder(Var1, -total_normal_cell_percent)) # Reorder Var1 based on total_normal_cell_percent


ggplot(sorted_graphdata, aes(x = Var1, y = percent, fill = Var2.y))+
  geom_bar(stat = "id entity") +
  facet_wrap(~ Var2.x,nrow = 1,scale = "free") +
  labs(
    x = "Category",
    y = "Percentage",
    fill = "Subgroup",
    title = "Morphology grouped by Family"
  )  +
  theme(axis.text.x = element_blank(),legend.title = element_blank(),
        legend.key = element_rect(color = "black", size = 0.1),  # Set legend key outline color
        legend.key.size = unit(1, "lines"),
        strip.text = element_text(size = 12)) +
  guides(alpha = "none")   # Remove alpha legend














