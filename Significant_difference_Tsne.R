library(dplyr)
library(ggplot2)
library(umap)
library(Rtsne)

# Extract significant profile columns
sig_cols <- sig_results$Variable

# Keep only significant profile variables
sig_profiledata <- data %>%
  select(all_of(sig_cols))

# Remove rows with missing values
keep <- complete.cases(sig_profiledata)

sig_profiledata <- sig_profiledata[keep, ]
plot_meta <- data[keep, ]

set.seed(123)

um <- umap(
  scale(sig_profiledata)
)

umap_df <- data.frame(
  UMAP1 = um$layout[,1],
  UMAP2 = um$layout[,2],
  Population = plot_meta$Population
)

ggplot(
  umap_df,
  aes(UMAP1, UMAP2, colour = Population)
) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_bw()

ggplot(
  umap_df,
  aes(UMAP1, UMAP2, colour = Population)
) +
  geom_point(size = 0.5, alpha = 0.5) +
  facet_wrap(~Population) +
  theme_bw()

set.seed(123)

tsne <- Rtsne(
  scale(sig_profiledata),
  perplexity = 30,
  check_duplicates = FALSE
)

tsne_df <- data.frame(
  TSNE1 = tsne$Y[,1],
  TSNE2 = tsne$Y[,2],
  Population = plot_meta$Population
)
tsne_df$combo <- paste(tsne_df$Population,
                       data$Hypoxic_exposure,
                       sep = "_")
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

a<-ggplot(
  tsne_df,
  aes(TSNE1, TSNE2, colour = combo)
) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_colour_manual(values = all_combo_colors) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
a
#1205 X 502 px
width_in  <- 803 / 300
height_in <- 502 / 300

width_in
# 2.7

height_in
# 1.67
ggsave(
  "population_TSNE.pdf",
  a,
  device = cairo_pdf,
  width = width_in,
  height = height_in,
  units = "in",
  bg = "white"
)
ggplot(
  tsne_df,
  aes(TSNE1, TSNE2, colour = combo)
) +
  geom_point(size = 0.5, alpha = 0.5)+
  scale_colour_manual(values = all_combo_colors) +
  facet_wrap(~Population) +
  theme_bw()

###########
#TSNE values cluster export
###########


cluster_df <- data.frame(
  CellID = plot_meta$CellID,
  Population = plot_meta$Population,
  TSNE1 = tsne_df$TSNE1,
  TSNE2 = tsne_df$TSNE2,
  Cluster = 0
)

#  cluster 1
cluster_df$Cluster[
    cluster_df$TSNE1 > -2 &
    cluster_df$TSNE1 < 8 &
    cluster_df$TSNE2 > 5 &
    cluster_df$TSNE2 < 15
] <- 1

# cluster 2
cluster_df$Cluster[
    cluster_df$TSNE1 > -20 &
    cluster_df$TSNE1 < -10 &
    cluster_df$TSNE2 > -10 &
    cluster_df$TSNE2 < 0
] <- 2

# Check cluster composition
table(
  cluster_df$Cluster,
  cluster_df$Population
)
export<-cbind(cluster_df$CellID,cluster_df$Cluster)
# Export
write.table(
  export,
  file = "Manual_TSNE_clusters.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

ggplot(
  cluster_df,
  aes(TSNE1, TSNE2, colour = Population)
) +
  geom_point(size = 0.5, alpha = 0.5) +
  annotate(
    "rect",
    xmin = -2,
    xmax = 8,
    ymin = 5,
    ymax = 15,
    fill = NA,
    colour = "red",
    linewidth = 1
  ) +
  annotate(
    "rect",
    xmin = -20,
    xmax = -10,
    ymin = -10,
    ymax = 0,
    fill = NA,
    colour = "blue",
    linewidth = 1
  ) +
  theme_bw()
