rm(list=ls())



## independent t-tests, OPLDA analyses,
## identify metabolites that were significantly depleted by each of the 12 species


library(tidyverse)
library(ropls)

abundance_data <- read.csv("annotated metabolites necromass consumption data.csv", row.names = 1, check.names = FALSE)
abundance_data <- as.matrix(abundance_data)


metadata <- read.csv("Sample information.csv",header=T)
meta2 <- metadata[61:125,]


meta2$Treatment <- factor(meta2$Treatment)
control_group <- "Mix"
treatment_groups <- setdiff(levels(meta2$Treatment), control_group)




for (i in 1:12) {
  treatment <- treatment_groups[i]
 

  treatment_samples <- meta2$Sample[meta2$Treatment == treatment]
  control_samples <- meta2$Sample[meta2$Treatment == control_group]
  

  subset_data <- abundance_data[, c(control_samples, treatment_samples)]
  

  mean_treatment <- rowMeans(subset_data[, treatment_samples])
  mean_control <- rowMeans(subset_data[, control_samples])

  fold_change <- mean_treatment / mean_control
  log10_fc <- ifelse(fold_change > 0, log10(fold_change), NA)
  

  group_factor <- factor(c(
    rep(control_group, length(control_samples)),
    rep(treatment, length(treatment_samples))
  ))
  

  oplsda_result <- opls(
    x = t(subset_data),
    y = group_factor,
    log10L=TRUE,
    predI = 1,        
    orthoI = 1       

  )
  

  vip_values <- oplsda_result@vipVn
  

  p_values <- apply(subset_data, 1, function(x) {
    tryCatch(
      t.test(x ~ group_factor,paired=F,var.equal=T)$p.value,
      error = function(e) NA
    )
  })
  

  fdr <- p.adjust(p_values, method = "fdr")
  

  result_df <- data.frame(
    id = rownames(subset_data),              # 1. 代谢物ID
    Mean_Treatment = mean_treatment,         # 2. 处理组均值
    Mean.mix = mean_control,                 # 3. 对照组均值
    VIP.OPLDA = vip_values,                  # 4. VIP值
    p.value = p_values,                      # 5. p值
    FDR = fdr,                               # 6. FDR
    FoldChange = fold_change,                # 7. FoldChange
    log10FC = log10_fc,                      # 8. log10FC
    stringsAsFactors = FALSE
  )
  

  colnames(result_df)[2] <- paste0("Mean_", treatment)
  

  output_file <- paste0("Spent_", treatment, "_vs_", control_group, ".csv")
  write.csv(result_df, output_file, row.names = FALSE)

}
