#' ## Batch effect correction and reclustering Endothelial
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
endo = normctxp5[,normctxp5$celltype=="Endothelial"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(endo, 
                        method = "CV", 
                        threshold = 2, 
                        minCells = (0.01*ncol(endo)),  
                        cutoff = 0.9)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
endo = dim_reduce(endo, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.87)

# Define initial clusters by density peak
endo = cluster_sc(endo, dimension = "2d", method = "density", num_clust = 2)
plot_tsne_metadata(endo, color_by = "Cluster")
plot_tsne_metadata(endo, color_by = "beads")
plot_tsne_metadata(endo, color_by = "lib_size")
plot_tsne_gene(endo,"Ncan")
plot_tsne_gene(endo,"Slc1a3")
plot_tsne_gene(endo,"Fabp7")
plot_tsne_gene(normctxp5,"Ncan", log_scale = T)
plot_tsne_gene(normctxp5,"Slc1a3", log_scale = T)
plot_tsne_gene(normctxp5,"Fabp7", log_scale = T)

#' ## results up

rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
colnames(pData(endo)) = gsub("^","endo_",colnames(pData(endo)))
merged = merge(pData(rawctxp5), pData(endo), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
plot_tsne_metadata(rawctxp5, 
                   color_by = "endo_Cluster", 
                   xcol = "endo_x", ycol = "endo_y",
                   facet_by = "endo_Cluster")
plot_tsne_gene(normctxp5,"Ncan")

### Label subtypes
rawctxp5$subtype[rawctxp5$endo_Cluster=="Cluster1"] = "doublets"
rawctxp5$subtype[rawctxp5$endo_Cluster=="Cluster2"] = "Endothelial"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "endo_x", ycol = "endo_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")