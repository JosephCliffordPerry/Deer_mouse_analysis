library(Rtsne)
library(ggplot2)
library(dplyr)

data<-read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Combined_relabelled_deer_mice_stats.txt", sep = "\t", header = TRUE)
numdata <- dplyr::select_if(data, is.numeric)

TSNEdata <- numdata[, grepl("_profile_", names(numdata))]

tsneo <- Rtsne(
  as.matrix(scale(TSNEdata)),
  dims = 2,
  perplexity = 30,
  max_iter = 1000,
  theta = 0.5,
  check_duplicates = FALSE
)

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
  KN_2day    = "green",
  KN_2wk     = "seagreen",
  KN_6wk     = "#8E260B",
  ME_Control = "#689BCA",
  ME_2day    = "yellow",
  ME_2wk     = "purple",
  ME_6wk     = "#254869"
)
data_with_hypoxia <- dplyr::left_join(
  data_with_sample,
  Phenotype_data,
  by = "Sample_number"
)
graphdata <- data.frame(
  tsne1 = tsneo$Y[,1],
  tsne2 = tsneo$Y[,2],
  Population = data_with_hypoxia$Population,
  Hypoxia_label = data_with_hypoxia$Hypoxic_exposure
) %>%
  dplyr::mutate(
    combo = paste0(Population, "_", Hypoxia_label)
  )
ggplot(graphdata,
       aes(tsne1, tsne2, colour = combo)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = all_combo_colors) +
  theme_bw()



ggplot(graphdata,
       aes(tsne1, tsne2, colour = combo)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = all_combo_colors) +
  facet_wrap(~ combo, ncol = 4) +
  theme_bw() +
  theme(
    legend.position = "none"
  )+
  facet_wrap(~ combo, scales = "free")








pop_colors <- c(
  KN = "#D13610",
  ME = "#356897"
)

graphdata <- data.frame(
  tSNE1 = tsneo$Y[,1],
  tSNE2 = tsneo$Y[,2],
  Population = data_with_hypoxia$Population
) %>%
  filter(!is.na(Population))

ggplot(graphdata,
       aes(tSNE1, tSNE2, colour = Population)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = pop_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )


ggplot(graphdata,
       aes(tSNE1, tSNE2, colour = Population)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_colour_manual(values = pop_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  facet_wrap(~ Population) +
  theme_bw()
