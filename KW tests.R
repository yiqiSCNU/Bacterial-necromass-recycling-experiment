rm(list=ls())


library(tidyverse)
library(reshape2)



## Kruskal-Wallis test, 
## 1. identify metabolites with significant abundances across 12 single species necromass medium

metab_data <- read.csv("annotated metabolites necromass composition table.csv", header = TRUE, row.names = 1)

metadata <- read.csv("Sample information.csv",header = T)
meta1 <- subset.data.frame(metadata,Experiment=="composition")

groups <- meta1$Treatment
unique_groups <- unique(groups)
group_counts <- table(groups)  


results <- data.frame(
  Metabolite = rownames(metab_data),
  KW_statistic = numeric(nrow(metab_data)),
  KW_pvalue = numeric(nrow(metab_data)),
  stringsAsFactors = FALSE
)


for (i in 1:nrow(metab_data)) {
  values <- as.numeric(metab_data[i, ])
  test_result <- kruskal.test(values ~ groups)
  
  results$KW_statistic[i] <- test_result$statistic
  results$KW_pvalue[i] <- test_result$p.value
}

results$FDR <- p.adjust(results$KW_pvalue, method = "BH")


group_means <- metab_data %>%
  t() %>%
  as.data.frame() %>%
  mutate(Group = groups) %>%
  group_by(Group) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  column_to_rownames("Group") %>%
  t() %>%
  as.data.frame()



results$Max_group_mean <- apply(group_means, 1, function(x) max(x, na.rm = TRUE))
results$Min_group_mean <- apply(group_means, 1, function(x) min(x, na.rm = TRUE))


results$Min_nonzero_group_mean <- apply(group_means, 1, function(x) {

  nonzero_vals <- x[x > 0]

  if (length(nonzero_vals) == 0) return(NA)

  min(nonzero_vals, na.rm = TRUE)
})


results$Max_FC <- apply(group_means, 1, function(x) {

  max_val <- max(x, na.rm = TRUE)

  min_nonzero <- min(x[x > 0], na.rm = TRUE)

  if (is.infinite(min_nonzero) || min_nonzero <= 0 || length(min_nonzero) == 0) {
    return(NA)
  }

  max_val / min_nonzero
})


final_results <- cbind(results, group_means)


write.csv(final_results, "KW test for necromass composition data.csv", row.names = FALSE)



## 2. identify metabolites utilized by 12 species with significantly different efficiency

metab_data <- read.csv("depleted metabolites abundance table.csv", header = TRUE, row.names = 1) [,-(1:5)]

 
metadata <- read.csv("Sample information.csv", header = TRUE)
metadata <- metadata[66:125,]

groups <- metadata$Treatment
unique_groups <- unique(groups)
group_counts <- table(groups)  

results <- data.frame(
  Metabolite = rownames(metab_data),
  KW_statistic = numeric(nrow(metab_data)),
  KW_pvalue = numeric(nrow(metab_data)),
  stringsAsFactors = FALSE
)

for (i in 1:nrow(metab_data)) {
  values <- as.numeric(metab_data[i, ])
  test_result <- kruskal.test(values ~ groups)
  
  results$KW_statistic[i] <- test_result$statistic
  results$KW_pvalue[i] <- test_result$p.value
}


results$FDR <- p.adjust(results$KW_pvalue, method = "BH")

write.csv(results,"sig consummed KW test results.csv")
