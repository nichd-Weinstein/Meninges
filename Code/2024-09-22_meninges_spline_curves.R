##This code is to calculate combined spline curves for the progenitors to the EPD and mFB branches

##Load libraries
library(Seurat)
library(URD)
library(Matrix)
library(ggplot2)

##Calculate varying genes
#Consider all genes expressed in 1% of meningeal cells
frac.exp <- rowSums(as.matrix(obj.tree@logupx.data > 0))/ncol(obj.tree@logupx.data > 0)
meningeal.genes <- names(frac.exp)

##Calculate smoothed spline fits for expression
##Branch 2 - Leptomeninges
lm.spline <- geneSmoothFit(obj.tree, method = "spline", pseudotime = "pseudotime.2", 
                           cells = cellsInCluster(obj.tree, "segment", c("3", "2")), 
                           genes = meningeal.genes, moving.window = 1, cells.per.window = 5,
                           spar = 0.9)

##Branch 1 - fibroblasts
fb.spline <- geneSmoothFit(obj.tree, method = "spline", pseudotime = "pseudotime.2", 
                           cells = cellsInCluster(obj.tree, "segment", c("3", "1")), 
                           genes = meningeal.genes, moving.window = 4, cells.per.window = 8,
                           spar = 0.9)

##Find genes that change their actial mean expression value by 75%
lm.change.real <- apply(lm.spline$mean.smooth, 1, function(x) diff(range(x)))
fb.change.real <- apply(fb.spline$mean.smooth, 1, function(x) diff(range(x)))
lm.genes.mean <- names(which(lm.change.real >= 0.75))
fb.genes.mean <- names(which(fb.change.real >= 0.75))

##Which genes are well fit by the spline curve? Remove noise that is not well fit
lm.spline.fit <- apply(lm.spline$scaled.smooth - lm.spline$scaled.expression.red, 1, 
                       function(i) sum(i^2))
fb.spline.fit <- apply(fb.spline$scaled.smooth - fb.spline$scaled.expression.red, 1, 
                       function(i) sum(i^2))

lm.spline.fit.norm <- lm.spline.fit/ncol(lm.spline$scaled.smooth)
fb.spline.fit.norm <- fb.spline.fit/ncol(fb.spline$scaled.smooth)
lm.spline.wellfit <- names(which(lm.spline.fit.norm <= 0.02))
fb.spline.wellfit <- names(which(fb.spline.fit.norm <= 0.02))

##Which genes change in their scaled log2 mean expression value by atleast 33%? 
lm.change.scale <- apply(lm.spline$scaled.smooth, 1, function(x) diff(range(x)))
fb.change.scale <- apply(fb.spline$scaled.smooth, 1, function(x) diff(range(x)))

lm.genes.scale <- names(which(lm.change.scale >= lm.spline.fit.norm * 0.15/0.02 + 0.33))
fb.genes.scale <- names(which(fb.change.scale >= fb.spline.fit.norm * 0.15/0.02 + 0.33))

##Ensure that genes are fit by the spline curve significantly better than a flat line with slope 0
lm.w <- ncol(lm.spline$scaled.smooth) - 1
lm.weight <- diff(as.numeric(colnames(lm.spline$scaled.smooth))) * 1000
lm.spline.fit.weighted <- apply(lm.spline$scaled.smooth[, 1:lm.w] - lm.spline$scaled.expression.red[, 1:lm.w],
                                1, function(i) sum(lm.weight * i^2))

fb.w <- ncol(fb.spline$scaled.smooth) - 1
fb.weight <- diff(as.numeric(colnames(fb.spline$scaled.smooth))) * 1000
fb.spline.fit.weighted <- apply(fb.spline$scaled.smooth[, 1:fb.w] - fb.spline$scaled.expression.red[, 1:fb.w],
                                1, function(i) sum(fb.weight * i^2))

lm.flat.fit.weighted <- apply(lm.spline$scaled.expression.red[, 1:lm.w], 1, function(x) sum(lm.weight * (x - mean(x)) ^ 2))
lm.spline.fit.ratio <- log2(lm.flat.fit.weighted/lm.spline.fit.weighted)
lm.spline.fit.betterthanflat <- names(which(lm.spline.fit.ratio > 0.25))

fb.flat.fit.weighted <- apply(fb.spline$scaled.expression.red[, 1:fb.w], 1, function(x) sum(fb.weight * (x - mean(x)) ^ 2))
fb.spline.fit.ratio <- log2(fb.flat.fit.weighted/fb.spline.fit.weighted)
fb.spline.fit.betterthanflat <- names(which(fb.spline.fit.ratio > 0.25))

#Take intersection of these genes and use them as the "varying genes" in the heatmap and analysis
lm.genes <- intersect(intersect(intersect(lm.genes.scale, lm.genes.mean), lm.spline.wellfit), 
                      lm.spline.fit.betterthanflat)
fb.genes <- intersect(intersect(intersect(fb.genes.scale, fb.genes.mean), fb.spline.wellfit), 
                      fb.spline.fit.betterthanflat)

##Combine pieces of spline fits into a single list of multi-plotting identify the pesudotime of the branchpoint
pt.crop <- as.numeric(unlist(obj.tree@tree$segment.pseudotime.limits)[1])
#Crop according to the pseudotime of the branchpoint
prog.only.spline <- cropSmoothFit(c(lm.spline, fb.spline), pt.max = pt.crop)
lm.only.spline <- cropSmoothFit(lm.spline, pt.min = pt.crop)
fb.only.spline <- cropSmoothFit(fb.spline, pt.min = pt.crop)

##Combine the list
splines <- list(prog.only.spline, lm.only.spline, fb.only.spline)
names(splines) <- c("Progenitors", "LM", "mFB")

##Plot expression of genes and fit splines
genes.plot <- c("tgfbi", "ctsla", "zbtb16a", "col11a1a")

plotSmoothFitMultiCascade(smoothed.fit = splines, genes = genes.plot, scaled = T, 
                          alpha.data = 0.2, alpha.smooth = 1, 
                          colors = c("Progenitors" = "#17A589", "LM" = "#FF0000", "mFB" = "#0065FF"))

saveRDS(splines, file = "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/cascades/splines_EPD_mFB_Progenitors_combined.rds")
splines <- readRDS("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/cascades/splines_EPD_mFB_Progenitors_combined.rds")

##Now find TFs that are commonly expressed between the EPDs and mFBs or in Progenitors
##Load TF list
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")

##TFs expressed in EPDs
##Which LM specific markers are TFs?
markers.lm.tfs <- lm.genes[lm.genes %in% tf.list$Symbol]
markers.fb.tfs <- fb.genes[fb.genes %in% tf.list$Symbol]

tfs.lm <- lm.genes.mean[lm.genes.mean %in% tf.list$Symbol]
tfs.fb <- fb.genes.mean[fb.genes.mean %in% tf.list$Symbol]

plotSmoothFitMultiCascade(smoothed.fit = splines, genes = markers.lm.tfs, scaled = T, 
                          alpha.data = 0.2, alpha.smooth = 1, 
                          colors = c("Progenitors" = "#17A589", "LM" = "#FF0000", "mFB" = "#0065FF"))

plotSmoothFitMultiCascade(smoothed.fit = splines, genes = c("crebrf", "foxd1", "foxl2a", "pias4a", "sox5", "zeb2b", "zic2b", "alx4a", "epd"), scaled = T, 
                          alpha.data = 0.2, alpha.smooth = 1, 
                          colors = c("Progenitors" = "#17A589", "LM" = "#FF0000", "mFB" = "#0065FF"))

epd.tf.list <- c("foxc1a", "foxc1b", "foxd1", "foxf2a", "foxl2a", "klf2a", "msx1b", "osr1", "zeb2b", "hopx",
                 "six1a", "six2a", "sox9a", "tbx15", "twist1a", "twist1b", "zic1", "zic2b", "zic3", "zic5")

fb.tf.list <- c("bhlhe41", "bhlhe40", "cebpd", "pbx1b", "stat5a", "tfe3a", "ebf2", "klf15", "tbx15", "zic2a")

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/plots/meninges_LM_mFB_splines_selectTFs.pdf", width = 8, height = 8)
plotSmoothFitMultiCascade(smoothed.fit = splines, genes = c("foxl2a", "zeb2b", "six1a", "klf2a", "klf15", "bhlhe40", "sox9a"), scaled = T, 
                          alpha.data = 0.2, alpha.smooth = 1, 
                          colors = c("Progenitors" = "black", "mFB" = "#fb75b1", "LM" = "#17A589"))
dev.off()

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/plots/meninges_LM_mFB_impulse_progenitor_plots.pdf", width = 20, height = 20)
plotSmoothFitMultiCascade(smoothed.fit = splines, genes = unlist(unique(list(fb.tf.list, epd.tf.list))), scaled = T, 
                          alpha.data = 0.2, alpha.smooth = 2, 
                          colors = c("Progenitors" = "#17A589", "LM" = "#FF0000", "mFB" = "#0065FF"), ncol = 5)
dev.off()
