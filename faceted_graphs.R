library(dplyr)
library(ggplot2)

load("umapo.rda")



rawdata<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt", sep = "\t", header = TRUE)

Phenotype_data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt", sep = "\t", header = TRUE)
Phenotype_data$Hypoxic_exposure <- dplyr::recode(
  Phenotype_data$Hypoxic_exposure,
  "2w" = "2wk",
  "6w" = "6wk"
)

unified_format_image<-gsub("24_11_06", "slide", rawdata$Image)
Sample_number<- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)
data_with_sample<-cbind(rawdata,Sample_number)
all_combo_colors <- c(
  KN_Control = "#F48367",
  KN_2day    = "#EF552E",
  KN_2wk     = "#D13610",
  KN_6wk     = "#8E260B",
  ME_Control = "#689BCA",
  ME_2day    = "#407CB5",
  ME_2wk     = "#356897",
  ME_6wk     = "#254869"
)
data_with_hypoxia <- dplyr::left_join(
  data_with_sample,
  Phenotype_data,
  by = "Sample_number"
)
graphdata <- data.frame(
  UMAP1 = umapo$layout[,1],
  UMAP2 = umapo$layout[,2],
  Population = data_with_hypoxia$Population,
  Hypoxia_label = data_with_hypoxia$Hypoxic_exposure
) %>%
  dplyr::mutate(
    combo = paste0(Population, "_", Hypoxia_label)
  )
ggplot(graphdata,
  aes(UMAP1, UMAP2, colour = combo)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = all_combo_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  )+
  facet_wrap(~ Population)

