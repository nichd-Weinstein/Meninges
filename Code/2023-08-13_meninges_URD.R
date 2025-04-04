##### Collaboration with Marina and Brant 
## This code is to perform all analysis using the mesenchyme URD object and generate some plots for Figure 7 of Marina's paper.

##Load libraries
library(Seurat)
library(URD)

sample = "meninges"

save.path <- "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/obj/"
plot.path <- "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/plots/"

##Load endoderm object
obj1 <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_seurat.rds"))
markers <- readRDS(paste0(file = "~/Box/zfext/annotations_celltype_curated_newMama/", sample, "/obj_seurat/", sample, "_markers.rds"))
DimPlot(obj1, label = T)    
DimPlot(obj1, group.by = "stage.nice")

##Remove the hybrid cells between the 2 cell types for the trajectory analysis
###1. Get the cells for the EPD cluster and subset them into a sub object to get the UMAP coordinates
cells.epd <- WhichCells(obj1, idents = c("29")). ##The putative hybrid cells are part of the EPD cluster
obj.epd <- subset(obj1, cells = cells.epd)

##Make a dataframe of the UMAP coordinates of the subsetted object
df <- as.data.frame(obj.epd@reductions$umap@cell.embeddings)
df2 <- as.data.frame(obj.fb@reductions$umap@cell.embeddings)

##Plot UMAP and get the hybrid cells using the UMAP coordinates
DimPlot(obj.epd)
cells.hybrid <- rownames(df)[df$UMAP_1 >= 1.5]
cells.fb.hybrid <- rownames(df2)[df2$UMAP_1 <= 2.5]

##Check in the mesenchyme UMAP to make sure you have captured the hybrid cells
DimPlot(obj1, cells.highlight = cells.hybrid)

##Cells to subset for trajectory analysis
cells.epd.sub <- setdiff(cells.epd, cells.hybrid)
cells.fb <- WhichCells(obj1, idents = c("32"))
cells.progenitors <- WhichCells(obj1, idents = c("21"))

##Union of all these above cell ids
cells.meninges <- unlist(unique(list(cells.epd.sub, cells.fb, cells.progenitors)))
obj.meninges <- subset(obj1, cells = cells.meninges) ##Subset object

#Create function to transfer Seurat information to URD
seuratToURD2 <- function(seurat.object) {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    # Create an empty URD object
    ds <- new("URD")
    
    # Copy over data
    ds@logupx.data <- as(as.matrix(seurat.object@assays$RNA@data), "dgCMatrix")
    if(!any(dim(seurat.object@assays$RNA@counts) == 0)) ds@count.data <- as(as.matrix(seurat.object@assays$RNA@counts[rownames(seurat.object@assays$RNA@data), colnames(seurat.object@assays$RNA@data)]), "dgCMatrix")
    
    # Copy over metadata
    ## TO DO - grab info
    get.data <- NULL
    if (.hasSlot(seurat.object, "data.info")) { 
      get.data <- as.data.frame(seurat.object@assays$RNA@data.info)
    } else if (.hasSlot(seurat.object, "meta.data")) { 
      get.data <- as.data.frame(seurat.object@meta.data) 
    }
    if(!is.null(get.data)) {
      di <- colnames(get.data)
      m <- grep("res|cluster|Res|Cluster", di, value=T, invert = T) # Put as metadata if it's not the result of a clustering.
      discrete <- apply(get.data, 2, function(x) length(unique(x)) / length(x))
      gi <- di[which(discrete <= 0.015)]
      ds@meta <- get.data[,m,drop=F]
      ds@group.ids <- get.data[,gi,drop=F]
    }
    
    # Copy over var.genes
    if(length(seurat.object@assays$RNA@var.features > 0)) ds@var.genes <- seurat.object@assays$RNA@var.features
    
    # Move over tSNE projection
    if (.hasSlot(seurat.object, "tsne.rot")) {
      if(!any(dim(seurat.object@tsne.rot) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@tsne.rot)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("tsne" %in% names(seurat.object@reductions)) && !any(dim(seurat.object@reductions$tsne) == 0)) {
        ds@tsne.y <- as.data.frame(seurat.object@reductions$tsne@cell.embeddings)
        colnames(ds@tsne.y) <- c("tSNE1", "tSNE2")
      }
    }
    
    # Move over PCA results
    if (.hasSlot(seurat.object, "pca.x")) {
      if(!any(dim(seurat.object@pca.x) == 0)) {
        ds@pca.load <- seurat.object@pca.x
        ds@pca.scores <- seurat.object@pca.rot
        warning("Need to set which PCs are significant in @pca.sig")
      }
      ## TO DO: Convert SVD to sdev
    } else if (.hasSlot(seurat.object, "reductions")) {
      if(("pca" %in% names(seurat.object@reductions)) && !any(dim(Loadings(seurat.object, reduction = "pca")) == 0)) {
        ds@pca.load <- as.data.frame(Loadings(seurat.object, reduction = "pca"))
        ds@pca.scores <- as.data.frame(seurat.object@reductions$pca@cell.embeddings)
        ds@pca.sdev <- seurat.object@reductions$pca@stdev
        ds@pca.sig <- pcaMarchenkoPastur(M=dim(ds@pca.scores)[1], N=dim(ds@pca.load)[1], pca.sdev=ds@pca.sdev)
      }
    }
    return(ds)
  } else {
    stop("Package Seurat is required for this function. To install: install.packages('Seurat')\n")
  }
}

#Creating URD object
obj <- seuratToURD2(obj.meninges)

saveRDS(obj, paste0(save.path, sample, "_URD.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD.rds"))

stages <- unique(obj@meta$stage.group)
cells.each.stage <- lapply(stages, function(stage) rownames(obj@meta)[which(obj@meta$stage.group == stage)])

#Calculate TSNE and graph clusterings
obj <- calcPCA(obj)
obj <- calcTsne(obj, perplexity = 30, theta = 0.5)
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(15, 20, 30, 40, 50), do.jaccard = T, method = "Louvain")
obj <- graphClustering(obj, dim.use = "pca", num.nn = c(15, 20, 25, 30, 35, 40, 45), do.jaccard = T, method = "Infomap")

umap.coord <- obj.meninges@reductions$umap@cell.embeddings
head(umap.coord)
head(obj@tsne.y)
colnames(umap.coord) <- c("tSNE1", "tSNE2")
class(umap.coord)
class(obj@tsne.y)

obj@tsne.y <- as.data.frame(umap.coord)
plotDim(obj, "cdx1b")

sample = "meninges"
saveRDS(obj, paste0(save.path, sample, "_URD.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD.rds"))

##Further subset the data - get rid of all the outlier cells
df <- as.data.frame(umap.coord)
cells.to.remove <- rownames(df)[df$tSNE2 > -2]
cells.remove.again <- rownames(df)[df$tSNE1 < -1]
cells.remove.all <- unlist(unique(list(cells.to.remove, cells.remove.again, cells.fb.hybrid)))
cells.keep <- setdiff(colnames(obj@logupx.data), cells.remove.all)

##Subset URD object
obj <- urdSubset(obj, cells.keep = cells.keep)

#Plotting the clustering and assessing tSNE plot
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDim(obj, "Louvain-50", legend = T, plot.title = "Louvain_Jaccard Graph-based Clustering (15 NNs)", alpha = 1)

#Removing outliers
#Calculating a k-nn graph
message(paste0(Sys.time(), ": Calculating a k-nearest neighbour graph"))
obj <- calcKNN(obj)

#Plot cells according to their distance to their nearest and 20th nearest neighbours, and identify those with unusually large distances.
outliers <- knnOutliers(obj, nn.1 = 1, nn.2 = 12, x.max = 24, slope.r = 2, int.r = 5, slope.b = 0.8, int.b = 7.5, title = "Identifying outliers by k-NN distance")       

gridExtra::grid.arrange(grobs=list(#Plot some apoptotic markers
  plotDim(obj, "isg15", alpha = 0.4, point.size = 0.5), 
  plotDim(obj, "foxo3b", alpha = 0.4, point.size = 0.5),
  plotDim(obj, "gadd45aa", alpha = 0.4, point.size = 0.5),
  #Figure out which clusters correspond to these cells
  plotDimHighlight(obj, clustering = "Louvain-40", cluster = "12", legend = F)))

#Subset object to eliminate outliers
cells.keep <- setdiff(colnames(obj@logupx.data), c(outliers))
obj <- urdSubset(obj, cells.keep = cells.keep)

##Saving trimmed object
saveRDS(obj, paste0(save.path, sample, "_URD_trimmed.rds"))

#Calculate diffusion map
obj <- calcDM(obj, knn = 100, sigma.use = 7.5) ##sigma used 7.5

dm.8 <- obj@dm
stage.colors <- c("darkblue", "#FFCCCC", "#99CC00", "#33CC00", "cyan3",
                            "gold", "goldenrod", "darkorange", "indianred1", "indianred2", "plum", "deepskyblue2",
                            "lightgrey")


plotDimArray(obj, reduction.use = "dm", dims.to.plot = 1:18, outer.title = "Diffusion Map (Sigma 10, 100 NNs): Stage", label="stage.group", alpha = 0.4, plot.title="", legend=T, discrete.colors = stage.colors)
plotDim(obj, "stage.group", transitions.plot = 10000, plot.title="Developmental stage (with transitions)")

saveRDS(obj, paste0(save.path, sample, "_URD_withDM.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD_withDM.rds"))

#Calculate pseudotime
message(paste0(Sys.time(), ": Defining earliest stage as the root.cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == " 24-34"]

plotDim(obj, "epd", plot.title="epd")
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDimHighlight(obj, "stage.group", " 24-34", plot.title = "tip-cells")
markersAUCPR(obj, "1", clustering = "Louvain-10", exp.thresh = 1.0, frac.must.express = 0.6)

#Do the flood
message(paste0(Sys.time(), ": Do the pseudotime_flood"))
flood.result <- floodPseudotime(obj, root.cells = root.cells, n = 100, minimum.cells.flooded = 2, verbose = T)

#Save the result
message(paste0(Sys.time(), ": Saving the flood object"))
saveRDS(flood.result, paste0(save.path, "meninges_URD_flood.rds"))
flood.result <- readRDS(paste0(save.path, sample, "_URD_flood.rds"))

#Process pseudotime floods
message(paste0(Sys.time(), ": Processing pseudotime_flood"))
obj <- floodPseudotimeProcess(obj, flood.result, floods.name = "pseudotime")
pseudotimePlotStabilityOverall(obj)

#Inspect pseudotime using tSNE plot
message(paste0(Sys.time(), ": Plotting Pseudotime on tSNE plot"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDim(obj, "stage.group", legend = T, plot.title = "Developmental Stage", alpha = 0.5)
plotDim(obj, "pseudotime.2", plot.title = "Pseudotime")
plotDim(obj, "epd")
dev.off()

#Plotting distances
message(paste0(Sys.time(), ": Plotting Distances"))
pdf(file=paste0(plot.path, sample, "_pseudotime.pdf"), width = 8, height = 8)
plotDists(obj, "pseudotime", "stage.group", plot.title="Pseudotime by stage")
dev.off()

gg.data <- cbind(obj@pseudotime, obj@meta[rownames(obj@pseudotime), ])
#Plot
ggplot(gg.data, aes(x=pseudotime, color = stage.group, fill = stage.group)) + geom_density(alpha = 0.4)

#Save object
message(paste0(Sys.time(), ": Saving object"))
saveRDS(obj, file=paste0(save.path, sample, "_URD_withDMandPT.rds"))
obj <- readRDS(paste0(save.path, sample, "_URD_withDMandPT.rds"))

##Get the pseudotime ranges for mFB and EPDs
##Make 2 dataframes containing the pseudotimes for each cell type
df1 <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.epd.sub, ]
df2 <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.fb, ]
cells.progenitors <- setdiff(colnames(obj@logupx.data), unlist(unique(list(cells.epd.sub, cells.fb))))
df3 <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.progenitors, ]

##Calculate range of pseudotime values for each dataframe
epd.range <- max(df1) - min(df1)
fb.range <- max(df2) - min(df2)
##Difference in range
diff.range <- epd.range/fb.range

##Create a separate column called pseudotime.2
obj@pseudotime$pseudotime.2 <- NA
obj@pseudotime[rownames(obj@pseudotime) %in% cells.epd.sub, "pseudotime.2"] <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.epd.sub, "pseudotime"] 
obj@pseudotime[rownames(obj@pseudotime) %in% cells.progenitors, "pseudotime.2"] <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.progenitors, "pseudotime"]
obj@pseudotime[rownames(obj@pseudotime) %in% cells.fb, "pseudotime.2"] <- obj@pseudotime[rownames(obj@pseudotime) %in% cells.fb, "pseudotime"] * diff.range

##Save the object with the new normalized pseudotime
saveRDS(obj, file=paste0(save.path, sample, "_URD_withDMandPT_2.rds"))

#Part 3 - Determining Tips
message(paste0(Sys.time(), ": Cropping the cells from the final stage_group"))
#cells_120h <- grep("120", colnames(obj@logupx.data), value = T)
cells_120h <- rownames(obj@meta)[obj@meta$stage.group == "120"]
##Subset URD object to just 120 hpf cells
obj_120h <- urdSubset(obj, cells.keep = cells_120h)

#Perform PCA/tSNE on final stage cells
message(paste0(Sys.time(), ": Load the variable genes specific to this stage"))
var.genes.120h <- scan(var.path, sample, "120_var.txt", what = "character")
obj_120h@var.genes <- obj.meninges@assays$RNA@var.features

#Calculate PCA
message(paste0(Sys.time(), ": Calculating PCA"))
obj_120h <- calcPCA(obj_120h)

#Calculate tSNE
message(paste0(Sys.time(), ": Calculating tSNE"))
set.seed(18)
obj_120h <- calcTsne(obj_120h, perplexity = 30, theta = 0.5)

obj_120h <- graphClustering(obj_120h, num.nn = 40, do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(5, 8, 10, 15, 20, 30), do.jaccard = T, method = "Louvain")
obj_120h <- graphClustering(obj_120h, num.nn = c(10, 15, 20, 30, 40, 50), do.jaccard = T, method = "Infomap")

clusterings <- c(paste0("Infomap-", c(10, 15, 20, 30, 40, 50)), paste0("Infomap-", c(10, 15, 20, 30, 40, 50)))
clusterings <- c(paste0("Louvain-", c(10, 15, 20, 30, 40)))

for (c in clusterings) {
  plot(plotDim(obj_120h, c, legend = T))
}

#Looking at batch information
pdf(file=paste0(plot.path, sample, "_tSNE_batch.pdf"), width = 8, height = 8)
plotDim(obj_120h, "Louvain-40", plot.title = "Louvain-40_graph", legend = T, label.clusters = T)
dev.off()

##Calculate differential markers for clusters
clusters <- sort(unique(obj_120h@group.ids$`Louvain-40`))
pr.markers <- lapply(clusters, function(c) markersAUCPR(obj_120h, clust.1 = c, clustering = "Louvain-40", genes.use = obj_120h@var.genes))
names(pr.markers) <- clusters

plotDim(obj_120h, "epd", plot.title="EPD (Leptomeninges marker)")

#Make a set of data.frames to keep track during cluster assignment
I30.n <- length(unique(obj_120h@group.ids$`Louvain-40`))
I30.cluster.assignments <- data.frame(cluster = 1:I30.n, name = rep(NA, I30.n, name = rep(NA, I30.n), tip = rep(NA, I30.n)), row.names = 1:I30.n)

plotDot(obj_120h, "slc4a7", clustering = "Louvain-40")

##Assign cluster names
I30.cluster.assignments["1", "name"] <- "mFB" #Use as tip
I30.cluster.assignments["2", "name"] <- "EPD" #Use as tips
#I30.cluster.assignments["3", "name"] <- "meningeal progeniors" #Don't use as tip

#Generate final clusterings
message(paste0(Sys.time(), ": Combine clustering assignments from two clusterings"))
I30.cluster.assignments$clustering <- "Louvain-40"
#I20.cluster.assignments$clustering <- "Louvain-40"
cluster.assignments <- rbind(I30.cluster.assignments)

#Remove any clusters that weren't assigned an identity
message(paste0(Sys.time(), ": Removing clusters without an identity"))
cluster.assignments <- cluster.assignments[!is.na(cluster.assignments$name), ]

#Renumber clusters
cluster.assignments$cluster.new <- 1:nrow(cluster.assignments)

#Create blank clusterings for obj_120h
obj_120h@group.ids$clusters.120h.name <- NA
obj_120h@group.ids$clusters.120h.num <- NA

#Copy cell identities over for each cluster
for (i in 1:nrow(cluster.assignments)) {
  cells <- cellsInCluster(obj_120h, clustering = cluster.assignments[i, "clustering"],
                          cluster = cluster.assignments[i, "cluster"])
  obj_120h@group.ids[cells, "clusters.120h.name"] <- cluster.assignments[i, "name"]
  obj_120h@group.ids[cells, "clusters.120h.num"] <- as.character(cluster.assignments[i, "cluster.new"])
}

#Transfer clusterings to main object
obj@group.ids$`Louvain-40` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "tip.clusters"] <- obj_120h@group.ids$`Louvain-40`

obj@group.ids$`Cluster` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster"] <- obj_120h@group.ids$clusters.120h.name

obj@group.ids$`Cluster-Num` <- NA
obj@group.ids[rownames(obj_120h@group.ids), "Cluster-Num"] <- obj_120h@group.ids$clusters.120h.num

#Save objects
message(paste0(Sys.time(), ": Saving the 120h_seurat object"))
saveRDS(obj_120h, file = paste0(save.path, "meninges_URD_120h_2.rds"))
obj_120h <- readRDS(file = paste0(save.path, "meninges_URD_120h_2.rds"))

message(paste0(Sys.time(), ": Saving the full object with 120h clustering added"))
saveRDS(obj, file = paste0(save.path, "meninges_URD_withTips_2.rds"))

##Read in 120 hpf URD object with tips
obj <- readRDS(file= paste0(save.path, "meninges_URD_withTips.rds"))

plotDim(obj, "tip.clusters", label.clusters = T)

message(paste0(Sys.time(), ": Saving the data.frame with tips"))
write.csv(cluster.assignments, file = paste0(save.path, "meninges_tips-use_Louvain-40.csv"))

#Plot tips in diffusion map
obj@group.ids$pop <- NA


#PART-4
#Biased random walks
message(paste0(Sys.time(), ": Load previous saved object"))
obj <- readRDS(file= paste0(save.path, "meninges_URD_withTips.rds"))
#Define parameters of logistic function to bias transition probablities
diffusion.logistic <- pseudotimeDetermineLogistic(obj, "pseudotime.2", optimal.cells.forward = 10, max.cells.back = 20, pseudotime.direction = "<",
                                                  do.plot = T, print.values = T)

#Create biased transition matrix
message(paste0(Sys.time(), ": Creating biased transition matrix"))
biased.tm <- pseudotimeWeightTransitionMatrix(obj, pseudotime = "pseudotime.2",
                                              logistic.params = diffusion.logistic, pseudotime.direction = "<")

#Define the root cells
message(paste0(Sys.time(), ": Defining root-cells"))
root.cells <- rownames(obj@meta)[obj@meta$stage.group == " 24-34"]

#Define tip cells
message(paste0(Sys.time(), ": Defining tip-cells"))
clustering <- "Louvain-40"
tips <- setdiff(unique(obj@group.ids[, clustering]), NA)
this.tip <- tips[tip.to.walk]

#Simulate the biased random walks from each tip
message(paste0(Sys.time(), ": Simulating random walks from each tip"))
tip.walks <- simulateRandomWalksFromTips(obj, tip.group.id = "tip.clusters", root.cells = root.cells,
                                         transition.matrix = biased.tm, n.per.tip = 25000, root.visits = 1,
                                         max.steps = 5000, verbose = T)

walks <- lapply(rownames(table(obj@group.ids$`Cluster-Num`)), function(c) {
  # Exclude any tip cells that for whatever reason didn't end up in the
  # biased TM (e.g. maybe not assigned a pseudotime).
  tip.cells <- rownames(obj@group.ids)[which(obj@group.ids$tip.clusters == c)]
  tip.cells.good <- intersect(tip.cells, rownames(biased.tm))
  # Perform the random walk simulation
  this.walk <- simulateRandomWalk(start.cells = tip.cells.good, transition.matrix = biased.tm,
                                  end.cells = root.cells, n = 50000, end.visits = 1, verbose.freq = 1000,
                                  max.steps = 5000)
  return(this.walk)
})
names(walks) <- rownames(table(obj@group.ids$`Cluster-Num`))

saveRDS(walks, file = paste0(save.path, sample, "_walks.rds"))
saveRDS(tip.walks, file = paste0(save.path, sample, "_tipWalks_2.rds"))
#Process the biased random walks into visitation frequencies
message(paste0(Sys.time(), ": Processing the biased random walks"))
obj <- processRandomWalksFromTips(obj, tip.walks, verbose = T)

#Visualize visitation of cells from each tip
plotDim(obj, "tip.clusters", plot.title = "Cells in each tip")

plotDim(obj, "visitfreq.log.1", plot.title = "Visitation frequency from tip-1 (log10)", transitions.plot = 10000)
plotDim(obj, "visitfreq.log.2", plot.title = "Visitation frequency from tip-2 (log10)", transitions.plot = 10000)

saveRDS(obj, file=paste0(save.path, sample, "_URD_withWalks_2.rds"))
saveRDS(obj, file=paste0(save.path, "/", folder_name, "/", dataset, "_URD_withWalks_more_var_genes.rds"))
obj <- readRDS(file=paste0(save.path, sample, "_URD_withWalks.rds"))
obj <- readRDS(file=paste0(save.path, sample, "_URD_withWalks_withTuftCells.rds"))

#PART 5 - Building URD Tree
library(URD)
library(rgl)

#Set up knitr to capture rgl output
rgl::setupKnitr()

#Load tip cells
message(paste0(Sys.time(), ": Loading tip cells"))
obj <- loadTipCells(obj, tips = "Cluster-Num")

#Build the actual tree
message(paste0(Sys.time(), ": Decide on tips to use for tree construction"))
tips.to.exclude <- c("3")
tips.to.use <- setdiff(as.character(1:2), tips.to.exclude)
tips.to.use <- c("1", "2")

#Build the tree
obj.built <- buildTree(object = obj, pseudotime = "pseudotime", divergence.method = "ks",
                       tips.use = tips.to.use, weighted.fusion = T, use.only.original.tips = T,
                       cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 5, minimum.visits = 1,
                       visit.threshold = 0.7, p.thresh = 0.001, save.breakpoint.plots = NULL, dendro.node.size = 100,
                       min.cells.per.segment = 10, min.pseudotime.per.segment = 0.1, verbose = F)

obj.tree <- buildTree(obj, pseudotime = "pseudotime", tips.use = tips.to.use, divergence.method = "preference", cells.per.pseudotime.bin = 25,
                      bins.per.pseudotime.window = 5, save.all.breakpoint.info = T, p.thresh = 0.1, verbose = F)

#Name the tips
tip.names <- unique(obj@group.ids[, c("Cluster", "Cluster-Num")])
tip.names <- tip.names[complete.cases(tip.names), ]
obj.tree <- nameSegments(obj.tree, segments = tip.names$`Cluster-Num`, segment.names = tip.names$`Cluster`)

##Finally, plot the tree
plotTree(obj.tree, "stage.group", label.segments = T, cell.alpha = 0.8, cell.size = 2, title = "Developmental_Stage_tree")

##Plot some genes
genes.plot <- c("epd", "ggctb", "ctsla", "sox4b")
for (gene in genes.plot) {
  plot(plotTree(obj.tree, gene))
}

pdf(paste0(plot.path, sample, "meninges_tree_genes.pdf"), width = 48, height = 48)
gridExtra::grid.arrange(grobs = lapply(genes.to.plot, plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T), ncol = 6)
dev.off()


##Check where all the cells from the Seurat object lie in the URD tree
#cells.15 <- WhichCells(obj1, idents = c("15"))

##Save the tree object
saveRDS(obj.tree, file=paste0(save.path, sample, "_TREE_v2.rds"))
obj.tree <- readRDS(file=paste0(save.path, sample, "_TREE.rds"))

#Manual refinement
message(paste0(Sys.time(), ": Descriptive names to be used on the dendrogram"))
new.seg.names <- c("mFB", "LM")
segs.to.name <- c("1", "2")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)

plotTree(obj.tree, "stage.group", title = "Developmental_stage_tree", label.segments = T)

diff.exp <- function(segment.1, segment.2) {
  markers.comp <- markersAUCPR(obj.tree, cells.1 = obj.tree@tree[["cells.in.segment"]][[segment.1]], cells.2 = obj.tree@tree[["cells.in.segment"]][[segment.2]], effect.size = 0.5, frac.min.diff = 0.1)
  
  write.csv(markers.comp, file = paste0(save.path, sample, "_markers_", segment.1, "vs", segment.2, ".csv"))
}
diff.exp("2", "1")

##Draw FDL or another layout to visualize trajectory
np.layout <- branchpointPreferenceLayout(obj.tree, pseudotime = "pseudotime", lineages.1 = "1", lineages.2 = "2", parent.of.lineages = "3",  min.visit = 2)
saveRDS(np.layout, file=paste0(save.path, sample, "_URD_branchpoint_preference.rds"))
np.layout <- readRDS(file=paste0(save.path, sample, "_URD_branchpoint_preference.rds"))
tree.colors <- c("#0000FF", "#2b8cbe", "#17A589", "#00FF00", "#f4be1d", "#8b5500", "#FF7F00", "#FF1493", "#FF00FF", "#8B008B")
plotBranchpoint(obj.tree, np.layout, label = "stage.group", point.alpha = 0.5, populations = c("mFB", "LM"),
                pt.lim = c(0.7, 0.1), xlab = "", ylab = "", legend = T, axis.lines = F,
                fade.low = 0, title = "Segment", visited.size = T, discrete.colors = tree.colors)
genes.plot <- c("foxc1a", "foxc1b", "foxf2b", "foxd2", "foxd1", "sox9a")
branch.plots <- lapply(genes.plot, function(gene) plotBranchpoint(obj.tree, np.layout, label = gene, point.alpha = 1, populations = c("mFB", "LM"),
                                                                  pt.lim = c(0.7, 0.11), xlab = "", ylab = "", title = gene, legend = F,
                                                                  axis.lines = F, fade.low = 0.66))
gridExtra::grid.arrange(grobs = branch.plots, ncol = 3)

gridExtra::grid.arrange(grobs = lapply(c("tgfbr2b", "fzd7a", "tfdp2",  "lin28a", "nfixb", "mn1b"), plotTree,
                                       object = obj.tree, label.x = F, plot.cells = T, color.limits = c(0,1)), ncol = 3)



##################### CASCADES ALONG THE TWO MENINGES SUBTYPES #############################
#Determine tips to run DE for
tips.to.run <- as.character(obj.tree@tree$segment.names)
genes.use <- NULL #Calculate for all genes

#Calculate the markers for each of the populations
gene.markers <- list()
for(tipn in 1:length(tips.to.run)) {
  tip <- tips.to.run[[tipn]]
  print(paste0(Sys.time(), ":", tip))
  markers <- aucprTestAlongTree(obj.tree, pseudotime = "pseudotime", tips = tip, log.effect.size = 0.4, auc.factor = 0.6, max.auc.threshold = 0.85, frac.must.express = 0.1, frac.min.diff = 0, genes.use = genes.use, root = "3", only.return.global = F, must.beat.sibs = 0.3, report.debug = T)
  saveRDS(markers, paste0(save.path, "/aucpr_markers/", tip, ".rds"))
  gene.markers[[tip]] <- markers
}

saveRDS(gene.markers, paste0(save.path, "/aucpr_markers/", "meninges_gene_markers_allTips.rds"))
gene.markers <- readRDS(paste0(save.path, "/aucpr_markers/", "meninges_gene_markers_allTips.rds"))

# Separate actual marker lists from the stats lists
gene.markers.de <- lapply(gene.markers, function(x) x[[1]])
gene.markers.stats <- lapply(gene.markers[1:8], function(x) x[[2]])
names(gene.markers.de) <- names(gene.markers)
names(gene.markers.stats) <- names(gene.markers)

## Examine the relationsip between DE genes and library complexity
# Compile all comparison stats into a single table
all.de.stats <- do.call("rbind", gene.markers.stats)
all.de.stats$tip <- substr(rownames(all.de.stats),1,nchar(rownames(all.de.stats))-2)
# Do a few plots
p1 <- ggplot(all.de.stats, aes(x=pt.1.mean, y=pt.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Pseudotime (Group 1)", y="Mean Pseudotime (Group 2)")
p2 <- ggplot(all.de.stats, aes(x=genes.1.mean, y=genes.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Detected Genes (Group 1)", y="Mean Detected Genes (Group 2)")
p3 <- ggplot(all.de.stats, aes(x=trans.1.mean, y=trans.2.mean)) + geom_point() + theme_bw() + geom_abline(slope = 1, intercept=0, col='red', lty=2) + labs(x="Mean Transcripts (Group 1)", y="Mean Transcripts (Group 2)")
cowplot::plot_grid(p1,p2,p3, ncol = 3)


# Determine temporal gene expression with impulse fitting
#Impulse fits
gene.cascades <- lapply(tips.to.run, function(tip) {
  print(paste0(Sys.time(), ": Impulse Fit ", tip))
  seg.cells <- cellsAlongLineage(obj.tree, tip, remove.root=F)
  casc <- geneCascadeProcess(object = obj.tree, pseudotime='pseudotime', cells = seg.cells, 
                             genes= rownames(gene.markers.de[[tip]]), 
                             moving.window=5, cells.per.window=18, 
                             pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)
  return(casc)
})
names(gene.cascades) <- tips.to.run
saveRDS(gene.cascades, file = "./data/allGene.cascades.rds")


##Load the calcuted cascade object
cascade.lm <- readRDS("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/cascades/meninges_LM_cascade.rds")
cascade.fb <- readRDS("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/cascades/meninges_mFB_cascade.rds")

##Calculate cascades on commonly shared genes between the two populations
casc.lm <- geneCascadeProcess(object = obj.tree, pseudotime='pseudotime', cells = cellsAlongLineage(obj.tree, c("3", "2"), remove.root=F), 
                   genes = unlist(unique(list(rownames(gene.markers.de$leptomeninges), markers.shared$gene, rownames(gene.markers.de$`meningeal fibroblasts`)))), 
                   moving.window=5, cells.per.window=18, 
                   pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)



##Plot the spline plots for differentially expressed transcription factors along the 2 lineages

##load the transcription factor list 
tf.list <- read.delim(file="~/Box/zfext/02-Clustering/2021-03 Iterative Clustering/gene_info/tfs/2021-06-25_zebrafish_LTA_TFs.txt", header = T, sep = "")
##Get differentially expressed TFs for the EPDs
genes <- rownames(cascade.lm$scaled.expression)
epd.tfs <- genes[which(genes %in% tf.list$Symbol)]
gene.num <- nrow(cascade.lm$scaled.expression)

##Get the TFs for fibroblasts
genes <- rownames(cascade.fb$scaled.expression)
fb.tfs <- genes[which(genes %in% tf.list$Symbol)]

##Plot branchpoint plots for some of these TFs
genes.plot <- fb.tfs[61:67]
pdf(paste0(plot.path, "meninges_mFB_specific_TFs_tree_61-67.pdf"), width = 48, height = 48)
branch.plots <- lapply(genes.plot, function(gene) plotBranchpoint(obj.tree, np.layout, label = gene, point.alpha = 1, populations = c("mFB", "LM"),
                                                                  pt.lim = c(0.7, 0.11), xlab = "", ylab = "", title = gene, legend = F,
                                                                  axis.lines = F, fade.low = 0.66))
gridExtra::grid.arrange(grobs = branch.plots, ncol = 4)

dev.off()

##plot spline curves with select genes 
genes.to.plot <- c("zeb2b", "zic2b", "sox5", "epd", "ggctb", "pias4a", "crebrf", "foxl2a", "foxd1")
##Plot Impulse plots
plotSmoothFit(smoothed.fit = cascade.lm, genes = genes.to.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1, plot.data = F)

##Plot similar spline curves for the fibroblast populations
genes.to.plot <- c("slc4a7", "slc47a2.1", "bhlhe41", "six2a", "smad7", "pbx1b", "tfe3a", "stat5a", "cebpd")
genes.to.plot <- c("alpl", "cebpb", "tfe3a", "zbtb4", "cebpg", "foxo3b", "tefa")
##Plot Impulse plots
plotSmoothFit(smoothed.fit = cascade.fb, genes = genes.to.plot, scaled = T, multiplot = F, alpha.data = 0.2, alpha.smooth = 1, plot.data = F)

##For each transcription factor, in each window plot two spline curves for EPDs and mFBs
cascades <- list(cascade.fb, cascade.lm)
names(cascades) <- c("mFB", "EPDs")
plotSmoothFitMultiCascade(smoothed.fits = cascades, genes = genes.to.plot, scaled = T, alpha.data = 0.2, alpha.smooth = 1)

##List all the classic markers of EPDs
epd.ph <- c("epd", "ggctb", "atp1b4", "slc38a3a", "soul5", "clu", "f3b", "fxyd1", 
            "vtnb", "slc13a4", "c4", "cl25b", "cp", "rbp4")

##Add bars  to annotate which genes are transcription factors
anno <- list(green=intersect(epd.tfs, genes),
             blue=intersect(epd.ph, genes))

##Define color scheme
pond.with.grey <- c("#CECECE", "#CBDAC2", RColorBrewer::brewer.pal(9, "YlGnBu")[3:9])

##Rerun the Impulse fitting on the transcription factors expressed in the progenitors and how they change expression along each trajectory
# Determine temporal gene expression with impulse fitting

##Change tip names to shorter names
message(paste0(Sys.time(), ": Shorter names to be used on the dendrogram"))
new.seg.names <- c("mFB", "LM")
segs.to.name <- c("1", "2")
obj.tree <- nameSegments(obj.tree, segments = segs.to.name, segment.names = new.seg.names)

##Which LM specific markers are TFs?
markers.lm.tfs <- rownames(cascade.lm$mean.smooth)[which(rownames(cascade.lm$mean.smooth) %in% tf.list$Symbol)]
markers.fb.tfs <- rownames(cascade.fb$mean.smooth)[which(rownames(cascade.fb$mean.smooth) %in% tf.list$Symbol)]

genes.plot <- unlist(unique(list("tsc22d3", "epd", "ggctb", "crebrf", "foxd1", "foxl2a", "pias4a", "sox5", "zeb2b", "zic2b", 
                "bhlhe41", "cebpd", "cebpd", "pbx1b", "six2a", "slc47a2.1", "slc4a7", "smad7", "stat5a", "tfe3a", 
                "id2b", markers.lm.tfs, markers.fb.tfs)))

tips.to.run <- c("mFB", "LM")
gene.cascades <- lapply(tips.to.run, function(tip) {
  print(paste0(Sys.time(), ": Impulse Fit ", tip))
  seg.cells <- cellsAlongLineage(obj.tree, tip, remove.root=F)
  casc <- geneCascadeProcess(object = obj.tree, pseudotime='pseudotime', cells = seg.cells, 
                             genes=genes.plot, 
                             moving.window=5, cells.per.window=18, 
                             pseudotime.per.window = 0.01, verbose = T, verbose.genes = T)
  return(casc)
})
names(gene.cascades) <- tips.to.run
saveRDS(gene.cascades, file = "~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/cascades/cascades_EPD_mFB_TFs.rds")

plots_per_page <- 9
pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/plots/meninges_LM_mFB_impulse_TFs_plots.pdf", width = 12, height = 12)
# Loop through the genes and create plots, 9 per page
for (i in seq(1, length(genes.plot), by = plots_per_page)) {
  
  # Select the genes for the current page
  genes_to_plot <- genes.plot[i:min(i + plots_per_page - 1, length(genes.plot))]
  
  # Create list to store individual plots
  plot_list <- list()
  
  # Generate individual plots for the current set of genes
  for (gene in genes_to_plot) {
    # Capture the plot generated by plotSmoothFitMultiCascade() for each gene
    p <- plotSmoothFitMultiCascade(smoothed.fits = gene.cascades, 
                                   genes = gene, 
                                   scaled = TRUE, 
                                   alpha.data = 0.2, 
                                   alpha.smooth = 1, 
                                   ncol = 1)
    # Add plot to list
    plot_list[[gene]] <- p
  }
  
  # Arrange the plots into a grid for the current page (ncol x nrow)
  grid.arrange(grobs = plot_list, ncol = 3, nrow = 3)
}
dev.off()

impulse.plots <- lapply(markers.prog.tfs, function(gene) plotSmoothFitMultiCascade(smoothed.fits = gene.cascades, 
                                                                                   genes = gene, scaled = T, 
                                                                                   alpha.data = 0.2, alpha.smooth = 1, ncol = 3))

pdf("~/Box/Farrell Lab/Manuscripts/2023 Brant and Marina_meninges/URD_trajectory/plots/meninges_LM_mFB_impulse_progenitor_plots.pdf", width = 20, height = 20)
gridExtra::grid.arrange(grobs = impulse.plots, ncol = 3)
dev.off()

