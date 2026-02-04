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

############################
# Base percentage table
############################
df <- plot_df %>%
  count(Sample, Morph, Population, Hypoxia_days) %>%
  group_by(Sample) %>%
  mutate(percent = n / sum(n) * 100) %>%
  ungroup() %>%
  mutate(
    Morph = factor(
      Morph,
      levels = c("normal_cells","small_cells","other_cells"),
      labels = c("Normal nuclei","Small nuclei","Abnormal nuclei")
    )
  )


############################################################
# -------- GRAPH A ORDERING --------
# Order by GENERAL % of each Morph within Population
# while visually grouping by Hypoxia
############################################################

# 1. Calculate general morphology ordering
order_a <- df %>%
  group_by(Morph, Population, Sample) %>%
  summarise(morph_percent = sum(percent), .groups = "drop") %>%
  group_by(Morph, Population) %>%
  arrange(desc(morph_percent)) %>%
  mutate(order_id = row_number()) %>%
  ungroup()

# 2. Join ordering back
df_a <- df %>%
  left_join(order_a, by = c("Morph","Population","Sample"))

# 3. Create X axis grouped by Hypoxia but ordered by percentage
df_a <- df_a %>%
  arrange(Morph, Population, Hypoxia_days, order_id) %>%
  mutate(
    x_id = interaction(Morph, Population, Hypoxia_days, Sample, drop = TRUE),
    x_id = factor(x_id, levels = unique(x_id))
  )

############################
# GRAPH A
############################
a <- ggplot(
  df_a,
  aes(x = x_id, y = percent, fill = factor(Hypoxia_days))
) +
  geom_bar(stat = "identity") +
  facet_nested(. ~ Morph + Population,
               scales = "free_x",
               space  = "free_x") +
  scale_x_discrete(labels = function(x) sub("^.*\\.", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(
    title = "Samples grouped by hypoxia and ordered by morphology %",
    x = "Sample",
    y = "Percentage",
    fill = "Hypoxia (days)"
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 11)
  )


############################################################
# -------- GRAPH B ORDERING --------
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

############################
# GRAPH B
############################
b <- ggplot(
  df_b,
  aes(x = reorder(x_id, -order_id), y = percent, fill = Morph)
) +
  geom_bar(stat = "identity") +
  facet_nested(. ~ Hypoxia_days + Population,
               scales = "free_x",
               space  = "free_x") +
  scale_fill_manual(values = c(
    "Abnormal nuclei" = "grey80",
    "Normal nuclei"   = "grey20",
    "Small nuclei"    = "grey50"
  )) +
  scale_x_discrete(labels = function(x) sub("^.*\\.", "", x)) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  labs(
    title = "Samples with morphology % faceted by hypoxia and population",
    x = "Sample",
    y = "Percentage"
  ) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x = element_text(size = 11),
    legend.title = element_blank()
  )

############################
# Combine plots
############################
a / b
