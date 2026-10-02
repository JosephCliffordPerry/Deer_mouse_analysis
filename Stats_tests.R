library(dplyr)

rawdata <- read.table(
  "E:/Deer_mouse_analysis/Clustered_deer_mice_stats.txt",
  sep = "\t",
  header = TRUE
)

Phenotype_data <- read.table(
  "E:/Deer_mouse_analysis/Deer_mouse_phenotype_data.txt",
  sep = "\t",
  header = TRUE
)



unified_format_image <- gsub(
  "24_11_06",
  "slide",
  rawdata$Image
)

Sample_number <- sub(
  "^.*?_(.*?)_.*$",
  "\\1",
  unified_format_image
)

data <- cbind(rawdata, Sample_number)

data <- merge(
  data,
  Phenotype_data,
  by = "Sample_number"
)

library(dplyr)

profiledata <- data %>%
  select(contains("_profile_"))

fit <- manova(
  as.matrix(profiledata) ~ Population,
  data = data
)

summary(fit, test = "Pillai")


fit <- manova(
  as.matrix(profiledata) ~ Population + factor(Hypoxic_exposure),
  data = data
)

summary(fit, test = "Pillai")


profile_cols <- names(profiledata)

results <- lapply(profile_cols, function(col){

  fit <- lm(
    reformulate("Population", response = col),
    data = data
  )

  data.frame(
    Variable = col,
    P = anova(fit)[["Pr(>F)"]][1]
  )

})

results <- do.call(rbind, results)

results$FDR <- p.adjust(
  results$P,
  method = "BH"
)

results <- results[order(results$FDR), ]

head(results, 100)
# export <- cbind(data_with_hypoxia$CellID[!is.na(data_with_hypoxia$Population)], as.numeric(as.factor(data_with_hypoxia$Population[!is.na(data_with_hypoxia$Population)])))
# export <- na.omit(export)
# write.table(export,file = "Populations.txt",sep = "\t",row.names = FALSE, col.names = FALSE,quote = FALSE)
results$Bonferroni <- p.adjust(
  results$P,
  method = "bonferroni"
)
table(results$Bonferroni)


sig_results <- results %>%
  filter(Bonferroni < 0.05)
