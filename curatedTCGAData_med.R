library(curatedTCGAData)
library(tidyverse)
library(psych)
library(ComplexHeatmap)
library(viridis)
library(glue)
library(ggpointdensity)

# Get PRAD RNA-seq data
data <- curatedTCGAData(diseaseCode = "PRAD", assays = "RNA*", dry.run = FALSE)

# Subset RNA-seq for AR and Med subunits
data_med <- assay(data[c("AR", "MED1", "MED4", "MED6", "MED7", "MED8", "MED9",
                         "MED10", "MED11", "MED12", "MED12L", "MED13", "MED13L",
                         "CDK8", "CDK19", "CCNC", "MED14", "MED15", "MED16", "MED17",
                         "MED18", "MED19", "MED20", "MED21", "MED22", "MED23",
                         "MED24", "MED25", "MED26", "MED27", "MED28", "MED29",
                         "MED30", "MED31")]) %>% 
  as.data.frame %>% 
  rownames_to_column("gene") %>%
  pivot_longer(!gene, names_to = "sample", values_to = "exp") %>%
  mutate(exp = log2(exp)) %>%
  pivot_wider(names_from = gene, values_from = exp)

cor <- corr.test(data_med[, -1], method = "spearman", adjust = "holm", use = "complete")

# Generate heatmap
pal <- viridis_pal()(256)

Heatmap(cor$r, col = pal, name = "Spearman's rho", 
        show_row_dend = FALSE, 
        show_column_dend = FALSE,
        border = TRUE, 
        heatmap_legend_param = list(at = c(-1, 0, 1), 
                                    direction = "horizontal"))

rho <- round(cor$r["AR", "MED12"], digits = 3)

data_med %>%
  ggplot(aes(x = AR, y = MED12)) +
  geom_pointdensity(size = 0.5) + 
  theme_bw() +
  scale_color_viridis_c() +
  geom_label(label = glue("Spearman's \u03c1 = {rho}"),
             x = 4, y = 11.4, vjust = 0, hjust = 0)
