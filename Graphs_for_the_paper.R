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
  Morph   = data$Dataset
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


# Plot
ggplot(df, aes(x = Sample, y = percent, fill = Morph)) +
  geom_bar(stat = "identity") +
  geom_segment(
    data = hypoxia_line,
    aes(x = x, xend = x,
        y = y_min, yend = y_max,
        color = factor(Hypoxia_days)),  # default ggplot2 colors
    inherit.aes = FALSE,
    size = 0.8
  ) +
  scale_fill_manual(values = c("Abnormal nuclei" = "grey80",
                               "Normal nuclei" = "grey20",
                               "Small nuclei" = "grey50")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Samples ordered by abnormal nuclei with hypoxia overlaid",
    x = "Sample",
    y = "Percentage",
    color = "Hypoxia (days)"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9)
  )

