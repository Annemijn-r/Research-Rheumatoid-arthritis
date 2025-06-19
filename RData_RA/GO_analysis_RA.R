library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
setwd(dir = "Transcriptoomanalyse/Casus/Research Rheumtaoid Arhtritis")
getwd()
  #Loading file in sigs2.
sigs2 <- read.csv("Results analysis dds/Results_case_RA.csv", sep = " ", , header = TRUE, row.names = 1)

  #Creating two files of up- and down-regulated.
up_genes2 <- rownames(sigs2[sigs2$log2FoldChange > 1, ])
down_genes2 <- rownames(sigs2[sigs2$log2FoldChange < -1, ])

  #GO analysis.
GO_up2 <- enrichGO(gene = up_genes2, 
                  OrgDb = 'org.Hs.eg.db', 
                  keyType = 'SYMBOL', #Symbols are used as gene names
                  ont = 'BP') #BP = biological processes

GO_down2 <- enrichGO(gene = down_genes2,
                    OrgDb = 'org.Hs.eg.db', 
                    keyType = 'SYMBOL', #Symbols are used as gene names
                    ont = 'BP') #BP = biological processes

GO_up_df2 <- as.data.frame(GO_up2)
GO_down_df2 <- as.data.frame(GO_down2)
  
#Saving GO-file.
write.csv2(GO_up_df2, file = "GO-analysis/Excel_results_GO/GO_enrichment_results_up.csv", row.names = FALSE)
write.csv2(GO_down_df2, file = "GO-analysis/Excel_results_GO/GO_enrichment_results_down.csv", row.names = FALSE)

  #Creating barplots.
fit_up <- plot(barplot(GO_up2, showCategory = 10) + #Top 10 will be shown
                   scale_fill_gradientn(
                     colours = c("#f497b6", "#90e0c4", "#b8a2d6"),
                     name = "p value"
                   ))
fit_up
png("GO-analysis/Barplot_GO/GO_Upregulated_results_plot.png", res = 250, width = 1800, height = 3200)
print(fit_up)  
dev.off()

fit_down <- plot(barplot(GO_down2, showCategory = 10) + #Top 10 will be shown
                   scale_fill_gradientn(
                     colours = c("#f497b6", "#90e0c4", "#b8a2d6"),
                     name = "p value"
                   )) 
fit_down
png("GO-analysis/Barplot_GO/GO_Downregulated_results_plot.png", res = 250, width = 2000, height = 3200)
print(fit_down)  
dev.off()


  #Analysing the results.
top_up_processes <- head(GO_up_df2[order(GO_up_df2$p.adjust), ], 10)
top_down_processes <- head(GO_down_df2[order(GO_down_df2$p.adjust), ], 10)

  #Showing all the genes included in the specific pathway.
top_up_processes$geneID[3]

top_down_processes$geneID[9]

  #Loading packages.
library(dplyr)
library(ggplot2)

  #Creating point diagram with counts and hits(%).
point_up_GO <- GO_up2@result %>%
  top_n(10, wt = -pvalue) %>%
  mutate(hitsPerc = Count * 100 / as.numeric(gsub("/.*", "", BgRatio))) %>%
  ggplot(aes(x = hitsPerc, 
             y = Description, 
             colour = pvalue, 
             size = Count)) +
  geom_point() +
  expand_limits(x = 0) +
  scale_colour_gradientn(
    colours = c("#90e0c4", "#b8a2d6", "#f497b6"),  
    name = "p value") +
  labs(x = "Hits (%)", y = "GO term upregulated", colour = "p value", size = "Count")

point_up_GO
png("GO-analysis/Dotplot_GO/GO_Upregulated_results_dotplot.png", width = 10, height = 6, units = "in", res = 300)
print(point_up_GO)  
dev.off()

point_down_GO <- GO_down2@result %>%
  top_n(10, wt = -pvalue) %>%
  mutate(hitsPerc = Count * 100 / as.numeric(gsub("/.*", "", BgRatio))) %>%
  ggplot(aes(x = hitsPerc, 
             y = Description, 
             colour = pvalue,
             size = Count)) +
  geom_point() +
  expand_limits(x = 0) +
  scale_colour_gradientn(
    colours = c("#90e0c4", "#b8a2d6", "#f497b6"),
    name = "p value") +
  labs(x = "Hits (%)", y = "GO term downregulated", colour = "p value", size = "Count")

point_down_GO
png("GO-analysis/Dotplot_GO/GO_Downregulated_results_dotplot.png", width = 8, height = 6, units = "in", res = 300)
print(point_down_GO)  
dev.off()

