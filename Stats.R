library(dplyr)
morph_matrix <- BasePercentageTable %>%
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


summary(lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(BasePercentageTable, Morph == "Abnormal nuclei")))

summary(lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(BasePercentageTable, Morph == "Intermediate")))


summary(lm(percent ~ Population + factor(Hypoxia_days),
   data = subset(BasePercentageTable, Morph == "Normal nuclei")))


wilcox.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Abnormal nuclei")
)

wilcox.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Intermediate")
)

wilcox.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Normal nuclei")
)

fligner.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Abnormal nuclei")
)

fligner.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Intermediate")
)

fligner.test(
  percent ~ Population,
  data = subset(BasePercentageTable, Morph == "Normal nuclei")
)

library(dplyr)

variability <- BasePercentageTable %>%
  group_by(Morph, Population) %>%
  summarise(
    SD = sd(percent),
    IQR = IQR(percent),
    .groups = "drop"
  )

variability_wide <- variability %>%
  tidyr::pivot_wider(
    names_from = Population,
    values_from = c(SD, IQR)
  ) %>%
  mutate(
    SD_ratio = SD_KN / SD_ME,
    IQR_ratio = IQR_KN / IQR_ME
  )

variability_wide
