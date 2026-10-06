############################
# Libraries
############################
library(dplyr)
library(ggplot2)
library(ggh4x)
library(patchwork)

############################
# Load and prepare data
############################
rawdata <- read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Clustered_deer_mice_stats.txt",
                      sep = "\t", header = TRUE)

Phenotype_data <- read.table("D:/2024-2025_sperm_images/Deer_Mouse/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt",
                             sep = "\t", header = TRUE)

unified_format_image <- gsub("24_11_06", "slide", rawdata$Image)
Sample_number <- sub("^.*?_(.*?)_.*$", "\\1", unified_format_image)

data_with_sample <- cbind(rawdata, Sample_number)
data <- merge(data_with_sample, Phenotype_data, by = "Sample_number")

data$Dataset <- gsub("Group_1_Cluster_8","small_cells", data$Dataset)
data$Dataset <- gsub("Group_1_Cluster_7","normal_cells", data$Dataset)

data$Dataset <- ifelse(data$Dataset %in% c("small_cells","normal_cells"),
                       data$Dataset,
                       "other_cells")

plot_df <- data.frame(
  Sample     = data$Sample_number,
  Hypoxia    = data$Hypoxic_exposure,
  Morph      = data$Dataset,
  Population = data$Population
)

plot_df <- plot_df %>%
  mutate(Hypoxia_days = case_when(
    Hypoxia == "Control" ~ 0,
    Hypoxia == "2day"    ~ 2,
    Hypoxia == "2wk"     ~ 14,
    Hypoxia == "6wk"     ~ 42,
    TRUE ~ NA_real_
  ))
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
############################
# Base percentage table
############################
BasePercentageTable <- plot_df %>%
  count(Sample, Morph, Population, Hypoxia_days) %>%
  group_by(Sample) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    Morph = factor(
      Morph,
      levels = c("normal_cells","small_cells","other_cells"),
      labels = c("Normal nuclei","Intermediate","Abnormal nuclei")
    )
  )


############################################################
# -------- GRAPH ORDERING --------
# Ordered by abnormal % inside Hypoxia × Population
############################################################

order_b <- df %>%
  filter(Morph == "Abnormal nuclei") %>%
  group_by(Hypoxia_days, Population) %>%
  arrange(desc(percent)) %>%
  mutate(order_id = row_number()) %>%
  ungroup() %>%
  select(Sample, Hypoxia_days, Population, order_id)

df_b <- df %>%
  left_join(order_b, by = c("Sample","Hypoxia_days","Population")) %>%
  mutate(x_id = interaction(Hypoxia_days, Population, Sample, drop = TRUE))

df <- df %>%
  mutate(
    Hypoxia_label = case_when(
      Hypoxia_days == 0  ~ "Control",
      Hypoxia_days == 2  ~ "2day",
      Hypoxia_days == 14 ~ "2wk",
      Hypoxia_days == 42 ~ "6wk"
    ),
    combo = paste0(Population, "_", Hypoxia_label)
  )

legend_df <- data.frame(
  Morph = "Abnormal nuclei",
  xmin  = c(3.20, 3.35, 3.50, 3.65,
            3.20, 3.35, 3.50, 3.65),
  xmax  = c(3.35, 3.50, 3.65, 3.80,
            3.35, 3.50, 3.65, 3.80),
  ymin  = c(95,95,95,95,85,85,85,85),
  ymax  = c(100,100,100,100,90,90,90,90),
  fill  = c(
    "#F48367","#EF552E","#D13610","#8E260B",
    "#689BCA","#407CB5","#356897","#254869"
  )
)

label_df <- data.frame(
  Morph = "Abnormal nuclei",
  x = c(3.9, 3.9),
  y = c(97.5, 87.5),
  label = c("KN", "MBS")
)

legend_df$Morph <- factor(
  legend_df$Morph,
  levels = levels(df$Morph)
)

label_df$Morph <- factor(
  label_df$Morph,
  levels = levels(df$Morph)
)
base_size<-8
b <- ggplot(
  df,
  aes(
    x = factor(Hypoxia_days),
    y = percent,
    colour = combo,
  )
) +

  geom_point(
    position = position_jitter(width = 0.15, height = 0),
    size = 2,
    alpha = 0.7
  ) +

  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.5,
    fatten = 0,
    colour = "black"
  ) +

  facet_nested(
    . ~ Morph,
    scales = "free_x",
    space = "free_x"
  ) +

  scale_colour_manual(
    values = all_combo_colors,
    labels = c(
      KN_Control = "KN Control",
      KN_2day    = "KN 2 days",
      KN_2wk     = "KN 2 weeks",
      KN_6wk     = "KN 6 weeks",
      ME_Control = "MBS Control",
      ME_2day    = "MBS 2 days",
      ME_2wk     = "MBS 2 weeks",
      ME_6wk     = "MBS 6 weeks"
    ),
    name = NULL
  ) +

  theme_bw() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = base_size),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      margin = margin(t = 8),
      size = base_size
    ),
    axis.text.y = element_text(size = base_size),  # increase y-axis labels too
    strip.text.x = element_text(size = base_size+2),
    plot.margin = unit(c(0.4, 0.07, 0, 0.07), "cm"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.length.x = unit(0, "cm"),
    axis.ticks.length.y = unit(0, "cm"),
  ) +

  scale_x_discrete(
    labels = c(
      "0"  = "Control",
      "2"  = "2 days",
      "14" = "2 weeks",
      "42" = "6 weeks"
    )
  ) +

  ylab("Percentage") +
  xlab("Duration of hypoxic exposure")+
  geom_rect(
    data = legend_df,
    aes(xmin = xmin, xmax = xmax,
        ymin = ymin, ymax = ymax,
        fill = fill),
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2
  )+coord_cartesian(clip = "off")

b

# #1205 X 340 px
# width_in  <- 1205 / 300
# height_in <- 502 / 300
#
# width_in
# # 4.02
#
# height_in
# # 1.13
# ggsave(
#   "Morph_percentages.pdf",
#   b,
#   device = cairo_pdf,
#   width = width_in,
#   height = height_in,
#   units = "in",
#   bg = "white"
# )
