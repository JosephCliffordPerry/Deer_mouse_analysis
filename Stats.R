morph_matrix <- df %>%
  select(Sample, Population, Morph, percent) %>%
  tidyr::pivot_wider(
    names_from = Morph,
    values_from = percent
  )

library(vegan)

adonis2(
  morph_matrix[,c(
    "Normal nuclei",
    "Intermediate",
    "Abnormal nuclei"
  )] ~ Population,
  data = morph_matrix,
  permutations = 9999,
  method = "euclidean"
)


dist_mat <- dist(
  morph_matrix[,c(
    "Normal nuclei",
    "Intermediate",
    "Abnormal nuclei"
  )]
)

disp <- betadisper(
  dist_mat,
  morph_matrix$Population
)

anova(disp)


lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(df, Morph == "Abnormal nuclei"))

lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(df, Morph == "Intermediate"))


lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(df, Morph == "Normal nuclei"))


wilcox.test(
  percent ~ Population,
  data = subset(df, Morph == "Abnormal nuclei")
)

wilcox.test(
  percent ~ Population,
  data = subset(df, Morph == "Intermediate")
)

wilcox.test(
  percent ~ Population,
  data = subset(df, Morph == "Normal nuclei")
)

