library(tidyverse)
library(psych)
library(viridis)
library(ComplexHeatmap)

# AR vs. MED expression
# TCGA PRAD

med_ar <- read_tsv("denseDataOnlyDownload.tsv") %>%
    select(-samples) %>%
    drop_na

med_ar_corr <- corr.test(med_ar[, -1], method = "spearman")

my_pal <- viridis_pal()(256)

Heatmap(med_ar_corr$r, heatmap_legend_param = list(title = "Spearman"), 
        col = my_pal, show_row_dend = F, show_column_dend = F)

write.table(med_ar_corr$p, "PRAD_AR_MED_corr_pval.txt", sep="\t", 
            col.names = NA, row.names = T, quote = F)
