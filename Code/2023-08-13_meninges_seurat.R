### Collaboration with Marina and Brant 
## This code is to perform all analysis using the mesenchyme seurat object and generate some plots for Figure 7 of Marina's paper.

##Load libraries
library(Seurat)
library(Matrix)
library(ggplot2)
library(dplyr)

############### FIGURE 7A ####################
##Load mama
message(paste0(Sys.time(), ": Loading mama"))
mama <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_ds_slim_seurat_25Jan2022.rds")

# Load cell annotations to introduce for downsample purposes
cell.annot.path <- "~/Box/zfext/annotations_celltype_curated_newMama/celltype_annotations_df.tsv"
message(Sys.time(), ": Loading ", cell.annot.path)
cell.annot <- read.table(cell.annot.path, sep = " ", header = T, stringsAsFactors = F)
cell.annot$tissue <- unlist(lapply(strsplit(x = cell.annot$clust, split = "\\."), function(x) x[1]))
cell.annot$clust.sg <- paste0(cell.annot$clust, "_", gsub(" ", "", cell.annot$stage.group))

mama@meta.data <- cbind(
  mama@meta.data[,setdiff(colnames(mama@meta.data), colnames(cell.annot))],
  cell.annot[rownames(mama@meta.data),]
)

##Figure 7A - highlight in mama where the meninges cells are
### Set idents on mama
Idents(mama) <- mama@meta.data$clust
cells.meninges <- WhichCells(mama, idents = c("mese.21", "mese.29", "mese.32"))

##Make a slot in mama showing the three clusters colored in mama

mama@meta.data$meninges.clusters <- NA
mama@meta.data[WhichCells(mama, idents = c("mese.21")), "meninges.clusters"] <- "meninges precursors"
mama@meta.data[WhichCells(mama, idents = c("mese.29")), "meninges.clusters"] <- "leptomeninges"
mama@meta.data[WhichCells(mama, idents = c("mese.32")), "meninges.clusters"] <- "meningeal fibroblasts"
mama@meta.data[setdiff(WhichCells(mama), WhichCells(mama, idents = c("mese.21", "mese.29", "mese.32"))), "meninges.clusters"] <- "rest"

colors.meninges <- setNames(c("#17A589", "#FF7F00", "#0000FF"), c("meninges precursors", "leptomeninges", "meningeal fibroblasts"))

# Assuming umap_df is your dataframe and cell_group is the column defining groups
umap_df <- as.data.frame(mama@reductions$umap@cell.embeddings)
umap_df$meninges.clusters <- mama@meta.data$meninges.clusters
#for(i in rownames(umap_df)){
#  umap_df[i, "stage.nice"] <- mama@meta.data[i, "stage.nice"]
#}

# Plot all cells in gray
p <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(color = "#f5f7f6", size = 1)

# Overlay specific groups with different colors
overlay_colors <- c("grey75", "#17A589", "#FF7F00", "#0000FF")  # Define colors for different groups
groups_to_overlay <- unique(mama@meta.data$meninges.clusters) # Define groups to overlay

for (i in seq_along(groups_to_overlay)) {
  group <- groups_to_overlay[i]
  color <- overlay_colors[i]
  
  overlay_points <- umap_df[umap_df$meninges.clusters == group, ]
  p <- p + geom_point(data = overlay_points, aes(x = UMAP_1, y = UMAP_2), color = color, size = 1.5) + theme_minimal() +  # Choose a minimal theme
    theme(panel.grid = element_blank(),  # Remove grid lines
          panel.background = element_rect(fill = "white"))  # Set white background
}

# Print the plot
dpi <- 300
png("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/mama_meninges_cells_highlighted_v4.png", height = 5*dpi, width = 5*dpi)
p
dev.off()

dpi <- 300
png("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/mama_meninges_cells_highlighted_v2.png", height = 5*dpi, width = 5*dpi)
DimPlot(mama, group.by = "meninges.clusters", cols = colors.meninges, pt.size = 2, raster = F, na.value = "grey75", order = T) + NoAxes() + NoLegend()
dev.off()


########### FIGURE 7B #################
##Now load mesenchyme object
sample <- "mesenchyme"
obj <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))

##Highlight the meninges cells in the mesenchyme object
##Make a slot in the meninges object for these cells - including the rest of non-meningeal cells in the mesenchyme
obj@meta.data$meninges.clusters <- NA
obj@meta.data[WhichCells(obj, idents = c("21")), "meninges.clusters"] <- "meninges precursors"
obj@meta.data[WhichCells(obj, idents = c("29")), "meninges.clusters"] <- "leptomeninges"
obj@meta.data[WhichCells(obj, idents = c("32")), "meninges.clusters"] <- "meningeal fibroblasts"


dpi <- 300
png("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/mesenchyme_obj_meninges_cells_highlighted.png", height = 5*dpi, width = 5*dpi)
DimPlot(obj, group.by = "meninges.clusters", cols = colors.meninges, pt.size = 2, raster = F, na.value = "grey75") + NoAxes() + NoLegend()
dev.off()

##Find markers for the 3 categories above - meninges clusters compared to the rest of the meninges
######### FIGURE 7D ########
##Set Idents
Idents(obj) <- obj@meta.data$meninges.clusters

##Now find differentially expressed markers between the 2 cell types
markers <- FindAllMarkers(obj, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

##Get markers which are specific to the 2 populations
markers.lm <- markers[markers$cluster == "leptomeninges" & markers$pct.1 - markers$pct.2 >= 0.8, "gene"]
markers.mf <- markers[markers$cluster == "meningeal fibroblasts" & markers$pct.1 - markers$pct.2 >= 0.65, "gene"]
markers.shared <- markers[markers$cluster == "meningeal fibroblasts" & markers$pct.1 - markers$pct.2 <= 0.05 & markers$avg_log2FC >= 0.9, "gene"]
markers.rest <- markers[markers$cluster == "mesenchyme_rest" & markers$pct.1 - markers$pct.2 >= 0.2, "gene"]

##Plot select genes that are shared
genes.shared <- c("jdp2b", "btg2", "ppdpfb", "igfbp2a", "msx1b", "hexb", "zic2a", "fzd7b")

##Make a list of all genes to plot
genes.to.plot <- c("ctsla", "ppdpfb", "fxyd1", "fabp11a", "zic2a", "epd", "ggctb", "soul5", "slc38a3a", "atp1b4",
                   "slc4a7", "slc47a2.1", "alpl", "slc7a3a", "slc5a6a", "prrx1b", "col6a1", "col6a2", "mab21l2", "meis1b")

genes.to.plot <- c("prrx1b", "col6a1", "col6a2", "mab21l2", "meis1b", "slc4a7", "slc47a2.1", "alpl", "slc7a3a", "slc5a6a",
                   "soul5", "slc38a3a", "atp1b4", "epd", "ggctb", "fxyd1", "fabp11a", "zic2a", "ctsla", "ppdpfb")

Idents(obj) <- factor(x = Idents(obj), levels = c("meninges precursors", "leptomeninges", "meningeal fibroblasts", "mesenchyme_rest"))

##Plot Dotplot
pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/meninges_comparisons_dotplot_all_v2.pdf", width = 12, height = 18)
DotPlot(obj, features = unique(genes.to.plot), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/meninges_comparisons_dotplot_v1.pdf", width = 12, height = 18)
DotPlot(obj, features = c("igfbp2a", "fabp11a", "fxyd1", "tfr1a", "slc16a9a", "slc16a5b", "slc5a6a",
                          "slc4a7", "slc47a2.1", "vtnb", "cp", "c4", "slc13a4", "ggctb", "epd", "ctsla", "ppdpfb"), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

##Now plot a dotplot for the most differentially expressed genes
##Make smaller object for just the meninges
cells.meninges <- WhichCells(obj, idents = c("21", "29", "32"))
obj.meninges <- subset(obj, cells = cells.meninges)

##Add the meninges annotations to the smaller object
obj.meninges@meta.data$meninges.clusters <- NA
obj.meninges@meta.data[WhichCells(obj.meninges, idents = c("21")), "meninges.clusters"] <- "meninges precursors"
obj.meninges@meta.data[WhichCells(obj.meninges, idents = c("29")), "meninges.clusters"] <- "leptomeninges"
obj.meninges@meta.data[WhichCells(obj.meninges, idents = c("32")), "meninges.clusters"] <- "meningeal fibroblasts"

##Set Idents
Idents(obj.meninges) <- obj.meninges@meta.data$meninges.clusters

##Now find differentially expressed markers between the 2 cell types
markers <- FindAllMarkers(obj.meninges, assay = "RNA", logfc.threshold = 0.25, only.pos = T, min.pct = 0.25)
top20 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

##Get markers which are specific to the 2 populations
markers.lm <- markers[markers$cluster == "leptomeninges" & markers$pct.1 - markers$pct.2 >= 0.4, ]
markers.mf <- markers[markers$cluster == "meningeal fibroblasts" & markers$pct.1 - markers$pct.2 >= 0.4, ]
markers.shared <- markers[markers$cluster == "meningeal fibroblasts" & markers$pct.1 - markers$pct.2 <= 0.1, ]

##Save the marker lists
write.csv(markers.lm, "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/gene_lists/EPD_markers_df_pct_diff_0.4_seurat.csv")
write.csv(markers.mf, "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/gene_lists/fibroblasts_markers_df_pct_diff_0.4_seurat.csv")
write.csv(markers.shared, "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/gene_lists/shared_markers_df_pct_diff_0.1.csv")

genes.shared <- c("sparc", "ctsla", "ppdpfb", "igfbp2a", "msx1b", "hexb", "zic2a", "fzd7b")

##Make a list of all genes to plot
genes.to.plot <- unlist(unique(list(markers.lm, markers.mf, genes.shared)))

##Plot Dotplot
pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/meninges_comparisons_dotplot_v5_new_colors.pdf", width = 12, height = 14)
DotPlot(obj.meninges, features = c("tfr1a", "slc16a9a", "slc16a5b", "slc5a6a",
                                    "slc4a7", "slc47a2.1", "vtnb", "cp", "c4", "slc13a4", "ggctb", "epd", 
                                   "igfbp2a", "zic2a", "fabp11a", "fxyd1", "ctsla", "ppdpfb"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_gradient2(low = "#416DF2", mid = "#C1ABF4", high = "#E43847", midpoint = 4) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/meninges_comparisons_dotplot_v6_new_colors.pdf", width = 12, height = 14)
DotPlot(obj.meninges, features = c("tfr1a", "slc16a9a", "slc16a5b", "slc5a6a",
                                   "slc4a7", "slc47a2.1", "vtnb", "cp", "c4", "slc13a4", "ggctb", "epd", 
                                   "igfbp2a", "zic2a", "fabp11a", "fxyd1", "ctsla", "ppdpfb"), dot.scale = 18, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_gradient2(low = "#416DF2", mid = "#FFFFFF", high = "#E43847", midpoint = 4) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##Plot another dotplot for the adult genes that Marina sent
pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/meninges_leptomeninges_fibroblasts_comparisons_dotplot_v3_new_colors.pdf", width = 12, height = 24)
DotPlot(obj.meninges, features = c("cygb1", "col1a1a", "col1a1b", "zgc:158423", "cldn11a", "wu:fj16a03", "rbp4", "prodha", "slc5a6a", "vtnb", "XLOC-018373", "slc13a4", "slc22a7b.1", "nlgn2b", "si:dkey-166k12.1", 
                                   "slc15a2", "b3gnt7", "col18a1a", "col18a1b", "sat1a.1", "clu", "slc4a4a", "si:ch211-195b13.1", 
                                   "slc6a9", "soul5", "ppdpfa", "f3a", "slc38a4", "slc6a22.1", "ca4a", "ggctb", "epd", "aldh1a2", "ngfrb", "cldn11b", "prrx1b", "plat", 
                                   "alpl", "slc4a2b", "slc7a3a", "kcnq5a", "slc4a7", "slc47a2.1", "slc16a9a", "slc16a5b", "csf1a", "si:ch211-236l14.4", 
                                   "crabp2a", "crabp2b", "crhbp", "sost", "csf1b", "si:ch211-105c13.3", "coch", "ptgdsb.1", "ptgdsb.2"),
                                   dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_gradient2(low = "#416DF2", mid = "#C1ABF4", high = "#E43847", midpoint = 4) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/meninges_leptomeninges_fibroblasts_comparisons_dotplot_v4_new_colors.pdf", width = 12, height = 24)
DotPlot(obj.meninges, features = c("cygb1", "col1a1a", "col1a1b", "zgc:158423", "cldn11a", "wu:fj16a03", "rbp4", "prodha", "slc5a6a", "vtnb", "XLOC-018373", "slc13a4", "slc22a7b.1", "nlgn2b", "si:dkey-166k12.1", 
                                   "slc15a2", "b3gnt7", "col18a1a", "col18a1b", "sat1a.1", "clu", "slc4a4a", "si:ch211-195b13.1", 
                                   "slc6a9", "soul5", "ppdpfa", "f3a", "slc38a4", "slc6a22.1", "ca4a", "ggctb", "epd", "aldh1a2", "ngfrb", "cldn11b", "prrx1b", "plat", 
                                   "alpl", "slc4a2b", "slc7a3a", "kcnq5a", "slc4a7", "slc47a2.1", "slc16a9a", "slc16a5b", "csf1a", "si:ch211-236l14.4", 
                                   "crabp2a", "crabp2b", "crhbp", "sost", "csf1b", "si:ch211-105c13.3", "coch", "ptgdsb.1", "ptgdsb.2"),
        dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_colour_gradient2(low = "#416DF2", mid = "#FFFFFF", high = "#E43847", midpoint = 4) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))
dev.off()


##Are there any markers specific to the progenitors?
markers.prog <- rownames(markers)[markers$cluster == "meninges precursors"]

##Get the progenitor genes that are transcription factors
markers.prog.tfs <- markers.prog[markers.prog %in% tf.list$Symbol]


########### FIGURE 7C ###########
##Draw feature-plots
##Ok, now plot the featureplots with the new URD color scheme
# Define the genes to plot
colors_feature_plot <- c("#CECECE80", "#B2B2B280", "#7D9FD180", "#5A90E0", "#307DF0", "#0065FF", "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF",
                                    "#00DEFF", "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", "#FFC900", "#FFB700", "#FFA500", "#FF9200",
                                    "#FF8000", "#FF6D00", "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000")
genes.plot <- c("slc4a7")

# Create a list of plots
fplot_peri <- lapply(genes.plot, function(gene) {
  FeaturePlot(object = obj, features = gene, pt.size = 1.8)
})

# Set the gradient color scale
scale <- scale_color_gradientn(colors = colors_feature_plot, na.value = "#CECECE") 

# Customize the plot theme
theme <- theme(panel.background = element_rect(fill = "white"),
               panel.grid = element_blank(),
               axis.line = element_line(colour = "black"),
               axis.text = element_text(colour = "black"),
               axis.title = element_text(colour = "black"),
               legend.text = element_text(colour = "black"),
               legend.title = element_text(colour = "black"))

# Modify the plots with the custom theme and color scale
fplots <- lapply(fplot_peri, function(plot) {
  plot + scale + theme 
})

dpi <- 300
png(paste0("FeaturePlot_", genes.plot, "_fig.png"), width = 2*dpi, height = 2*dpi)
fplots
dev.off()

DotPlot(obj.meninges, features = c("vtnb", "XLOC-018373", "slc13a4", "slc22a7b.1", "nlgn2b", "si:dkey-166k12.1", 
                                   "slc15a2", "b3gnt7", "col18a1a", "sat1a.1", "clu", "slc4a4a", "si:ch211-195b13.1",
                                   "slc6a9", "soul5", "ppdpfa", "f3a", "slc38a4", "slc6a22.1", "ca4a", "ggctb", "epd"), dot.scale = 16, scale = F) + coord_flip() +
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) +
  scale_color_gradientn(colours = defaultURDContinuousColors(with.grey = TRUE, evenly.spaced = TRUE)) +
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white")))


##Supplementary Figure 6 stage and cell cycle distribution
cells.progenitors <- WhichCells(obj, idents = c("21"))
cells.epd <- WhichCells(obj, idents = c("29"))
cells.fibroblasts <- WhichCells(obj, idents = c("32"))

##Get stage distribution dataframe
plotCellSizeBar <- function(seurat.object, clust, clustering = "clust", stage.meta = "stage.group") {
  cells.in.clust <- rownames(seurat.object@meta.data)[seurat.object@meta.data[,clustering] == clust]
  cells.per.stage <- table(seurat.object@meta.data[, stage.meta])
  cells.clust.stage <- as.data.frame(table(seurat.object@meta.data[cells.in.clust, stage.meta]), stringsAsFactors = F)
  colnames(cells.clust.stage) <- c("Stage", "n")
  stage.groups.to.add <- setdiff(unique(seurat.object@meta.data$stage.group), cells.clust.stage$Stage)
  df <- as.data.frame(stage.groups.to.add)
  df$n <- 0
  colnames(df) <- c("Stage", "n")
  cells.clust.stage <- rbind(df, cells.clust.stage)
  cells.clust.stage$total <- cells.per.stage
  cells.clust.stage$n.norm <- cells.clust.stage$n / cells.clust.stage$total
  cells.clust.stage$n.norm.percent <- round(cells.clust.stage$n.norm / sum(cells.clust.stage$n.norm) * 100, digits = 1)
  this.plot <- ggplot(data = cells.clust.stage, aes(x = Stage, y = n.norm.percent)) + geom_bar(stat = "identity") + scale_y_continuous()
  this.plot <- this.plot + theme_bw() + labs(x = "Hours post-fertilization", y = "% of cluster cells (normalized)")
  this.plot <- this.plot + theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
  return(this.plot)
}

plot(plotCellSizeBar(mama, "mese.21", "clust", "stage.group"))
plot(plotCellSizeBar(mama, "mese.29", "clust", "stage.group"))
plot(plotCellSizeBar(mama, "mese.32", "clust", "stage.group"))

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/mese-21_normalized_cellcounts_per_cluster.pdf", width = 4, height = 3)
plot(plotCellSizeBar(mama, "mese.21", "clust", "stage.group"))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/mese-29_normalized_cellcounts_per_cluster.pdf", width = 4, height = 3)
plot(plotCellSizeBar(mama, "mese.29", "clust", "stage.group"))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/mese-32_normalized_cellcounts_per_cluster.pdf", width = 4, height = 3)
plot(plotCellSizeBar(mama, "mese.32", "clust", "stage.group"))
dev.off()


plotCellInfoBar <- function(seurat.object, plot.meta, clust, clustering = "clust", stage.meta = "stage.group", normalized = T, legend = T, colors = NULL, xlab = "Hours post-fertilization", ylab = "% of cells", labels = waiver()) {
  cells.in.clust <- rownames(seurat.object@meta.data)[seurat.object@meta.data[,clustering] == clust]
  data.plot <- seurat.object@meta.data[cells.in.clust, c(clustering, stage.meta, plot.meta)]
  xt <- table(data.plot[,stage.meta], data.plot[,plot.meta])
  stage.groups.to.add <- setdiff(unique(seurat.object@meta.data$stage.group), rownames(xt))
  df <- as.data.frame(stage.groups.to.add)
  df$cycling <- 0
  df$non_cycling <- 0
  colnames(df) <- c("stage", "cycling", "non-cycling")
  rownames(df) <- df$stage
  if (normalized) xt <- round(sweep(xt, 1, rowSums(xt), "/") * 100, digits = 1)
  xtm <- reshape2::melt(xt)
  df2 <- reshape2::melt(df)
  colnames(df2) <- c("Var1", "Var2", "value") 
  xtm <- rbind(df2, xtm)
  this.plot <- ggplot(data = xtm, aes(x = Var1, fill = Var2, y = value)) + geom_bar(stat="identity") + theme_bw()
  this.plot <- this.plot + labs(x = xlab, y = ylab, fill = "")
  this.plot <- this.plot + theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
  if (!legend) this.plot <- this.plot + theme(legend.position = "none")
  if (!is.null(colors)) this.plot <- this.plot + scale_fill_manual(values = colors, labels = labels)
  return(this.plot)
}

plotCellInfoBar <- function(seurat.object, plot.meta, clust, clustering = "clust", stage.meta = "stage.group", normalized = T, legend = T, colors = NULL, xlab = "Hours post-fertilization", ylab = "% of cells", labels = NULL) {
  cells.in.clust <- rownames(seurat.object@meta.data)[seurat.object@meta.data[,clustering] == clust]
  data.plot <- seurat.object@meta.data[cells.in.clust, c(clustering, stage.meta, plot.meta)]
  xt <- table(data.plot[,stage.meta], data.plot[,plot.meta])
  
  if (normalized) xt <- round(sweep(xt, 1, rowSums(xt), "/") * 100, digits = 1)
  xtm <- reshape2::melt(xt)
  this.plot <- ggplot(data = xtm, aes(x = Var1, fill = Var2, y = value)) + geom_bar(stat="identity") + theme_bw()
  this.plot <- this.plot + labs(x = xlab, y = ylab, fill = "")
  this.plot <- this.plot + theme(axis.text.x = element_text(angle = 45, hjust = 0.95))
  if (!legend) this.plot <- this.plot + theme(legend.position = "none")
  if (!is.null(colors)) {
    if (is.null(labels)) {
      this.plot <- this.plot + scale_fill_manual(values = colors)
    } else {
      if (is.null(names(labels))) names(labels) <- labels
      this.plot <- this.plot + scale_fill_manual(values = colors, labels = labels,  breaks = names(labels))
    }
  }
  return(this.plot)
}

##Load cell cycle scores
cc.score <- readRDS("~/Box/zfext/02-Clustering/mama+ds_integrated/merged_mama_cellCycle_score.rds")
mama@meta.data <- cbind(mama@meta.data, cc.score)

##Find which cells belong to G1/S phase and which belong to G2/M phase
cells.s.phase <- rownames(mama@meta.data)[mama@meta.data$s.score > 0]
cells.g2m.phase <- rownames(mama@meta.data)[mama@meta.data$g2m.score > 0]

##Define a slot in the global metadata assigning cells to the different cell.cycle phases
mama@meta.data$cc.phase <- NA
mama@meta.data[cells.s.phase, "cc.phase"] <- "G1/S-phase"
mama@meta.data[cells.g2m.phase, "cc.phase"] <- "G2M-phase"

colors.cellcycle <- c("#E96E3A", "#1D6D8B")
names(colors.cellcycle) <- c("cycling", "non-cycling")

cycling.labels <- c("Cycling", "Not Cycling")
names(cycling.labels) <- c("cycling", "non-cycling")

##Use the above information to classify cells into cycling and non-cycling groups
##Classify cells into either cycling or non-cycling based on the above cell cycle scores
cells.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G1/S-phase")]
cells.also.cycling <- rownames(mama@meta.data)[which(mama@meta.data$cc.phase == "G2M-phase")]

##Total cycling  cells
cells.cycling.total <- unlist(unique(list(cells.cycling, cells.also.cycling)))

##Non cycling  cells
#cells.non.cycling <- rownames(mama@meta.data)[which(is.na(mama@meta.data$cc.phase))]
cells.non.cycling <- setdiff(WhichCells(mama), cells.cycling.total)

##Add these cycling and non-cycling groups to the global dataset metadata
mama@meta.data$cc.status <- NA
mama@meta.data[cells.non.cycling, "cc.status"] <- "non-cycling"
mama@meta.data[cells.cycling.total, "cc.status"] <- "cycling"

##Add non-cycling cells to the "cc.phase" column
mama@meta.data[cells.non.cycling, "cc.phase"] <- "non-cycling"

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/figure_panels/2024-09_panel_edits/mese-32_cellcycle_per_cluster.pdf", width = 4, height = 3)
plot(plotCellInfoBar(seurat.object = mama, plot.meta = "cc.status", clust = "mese.32", legend = F, colors = colors.cellcycle, normalized = F, ylab = "Number of cells"))
dev.off()

