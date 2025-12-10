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

# Prepare plotting data
plot_df <- data.frame(
  Sample  = data$Sample_number,
  Hypoxia = data$Hypoxic_exposure,
  Morph   = data$Dataset,
  Population = data$Population
)

# Convert hypoxia to days
plot_df <- plot_df %>%
  mutate(Hypoxia_days = case_when(
    Hypoxia == "Control" ~ 0,
    Hypoxia == "2day"    ~ 2,
    Hypoxia == "2wk"     ~ 14,
    Hypoxia == "6wk"     ~ 42,
    TRUE ~ NA_real_
  ))

# Counts per sample and morphology
cell_counts <- plot_df %>%
  count(Sample, Morph)

sample_totals <- cell_counts %>%
  group_by(Sample) %>%
  summarise(Total = sum(n), .groups = "drop")

df <- left_join(cell_counts, sample_totals, by = "Sample") %>%
  mutate(percent = n / Total * 100)

# Order samples by "other_cells"
other_cells_order <- df %>%
  filter(Morph == "other_cells") %>%
  select(Sample, other_percent = percent)

df <- df %>%
  left_join(other_cells_order, by = "Sample") %>%
  mutate(Sample = reorder(Sample, -other_percent))


# Rename Morph levels for plotting
df <- df %>%
  mutate(Morph = factor(Morph,
                        levels = c("normal_cells", "small_cells","other_cells" ),
                        labels = c( "Normal nuclei", "Small nuclei","Abnormal nuclei")))

# Prepare hypoxia vertical lines (top 10% of each bar)
hypoxia_line <- plot_df %>%
  distinct(Sample, Hypoxia_days) %>%
  left_join(other_cells_order, by = "Sample") %>%   # ensures same sample order
  mutate(
    x = as.numeric(reorder(Sample, -other_percent)), # numeric position of each bar
    y_min = 90,  # start at 90% of bar
    y_max = 100  # end at top
  )

# hypoxia_abnormality_graph -----------------------------------------------
############################
# Libraries
############################
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(umap)

set.seed(123)  # for reproducibility

############################
# Prepare bar chart (top panel)
############################
plot_df <- data %>%
  mutate(Morph = case_when(
    Dataset == "normal_cells" ~ "Normal nuclei",
    Dataset == "small_cells"  ~ "Small nuclei",
    TRUE                      ~ "Abnormal nuclei"
  )) %>%
  select(Sample_number, Population, Morph)

# Counts per sample and morphology
cell_counts <- plot_df %>%
  count(Sample_number, Morph, Population)

sample_totals <- cell_counts %>%
  group_by(Sample_number) %>%
  summarise(Total = sum(n), .groups = "drop")

df_pop <- left_join(cell_counts, sample_totals, by = "Sample_number") %>%
  mutate(percent = n / Total * 100)

# Order samples by "Abnormal nuclei"
abnormal_order <- df_pop %>%
  filter(Morph == "Abnormal nuclei") %>%
  select(Sample_number, abnormal_percent = percent)

df_pop <- df_pop %>%
  left_join(abnormal_order, by = "Sample_number") %>%
  mutate(Sample_number = reorder(Sample_number, -abnormal_percent),
         Morph = factor(Morph, levels = c("Normal nuclei","Small nuclei","Abnormal nuclei")))

top_plot <- ggplot(df_pop, aes(x = Sample_number, y = percent, fill = Morph)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Population, nrow = 1, scales = "free_x") +
  scale_fill_manual(values = c(
    "Abnormal nuclei" = "grey80",
    "Normal nuclei" = "grey20",
    "Small nuclei" = "grey50"
  )) +
  labs(
    title = "Morphology composition per sample",
    x = "Sample",
    y = "Percentage"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_blank(),
    strip.text = element_text(size = 12)
  )

############################
# Prepare UMAP (bottom panel)
############################

# Select 1000 random rows per population
umap_df_sample <- data %>%
  group_by(Population) %>%
  slice_sample(n = 5000) %>%
  ungroup()

# Use Diameter_profile_0-99 for UMAP
diam_cols <- paste0("Diameter_profile_", 0:99)
umap_matrix <- umap_df_sample %>%
  select(all_of(diam_cols)) %>%
  as.matrix()

# Run UMAP
set.seed(123)
umap_res <- umap(umap_matrix)

umap_plot_df <- as.data.frame(umap_res$layout)
colnames(umap_plot_df) <- c("UMAP1","UMAP2")
umap_plot_df$Population <- umap_df_sample$Population

bottom_plot <- ggplot(umap_plot_df, aes(x = UMAP1, y = UMAP2, color = Population)) +
  geom_point(size = 1) +
  labs(title = "UMAP of randomly sampled cells (5000 per Population)") +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

############################
# Combine panels
############################
top_plot / bottom_plot + plot_layout(heights = c(1,2))

