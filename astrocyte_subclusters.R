#' ## Batch effect correction and reclustering astro
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
astro = normctxp5[,normctxp5$celltype=="Astrocytes"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(astro, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(astro)),  
                        cutoff = 0.75)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
astro = dim_reduce(astro, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.93)

# Define initial clusters by density peak
astro = cluster_sc(astro, dimension = "2d", method = "density", s = 1)
plot_tsne_metadata(astro, color_by = "Cluster")
plot_tsne_metadata(astro, color_by = "beads")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(astro[,astro$beads=="V2"])),
                         log1p(exprs(astro[,astro$beads=="V3"])),
                         subset.row = vargenes,
                         d=14)

mnnastro = astro
pData(mnnastro) = pData(mnnastro)[,grep("iPC",colnames(pData(mnnastro)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(astro))[grep("iPC", colnames(pData(astro)))]
pData(mnnastro) = cbind(pData(mnnastro),pcs)
tSNE_result = Rtsne::Rtsne(pcs, dims = 2, perplexity = 30, theta = 0.5, check_duplicates = F, pca = F, max_iter = 500, verbose = T)
mnnastro$x = tSNE_result$Y[,1]
mnnastro$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnastro = cluster_sc(mnnastro, dimension = "2d", method = "density", s = 2)
plot_tsne_metadata(mnnastro, color_by = "Cluster")
plot_tsne_metadata(mnnastro, color_by = "beads")
mnnastro = calc_agg_bulk(mnnastro, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

cellcut = min(table(mnnastro$Cluster)/ncol(mnnastro))
rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            pd = pData(mnnastro),
            DEgroup="Cluster",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads",
            minCells = cellcut,
            outdir="results/DEgenes_markers/Astro/")

colnames(pData(mnnastro)) = gsub("^","astro_",colnames(pData(mnnastro)))
merged = merge(pData(rawctxp5), pData(mnnastro), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")

rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "astro_Cluster")
markers = scan("results/DEgenes_markers/Astro/results_markers", what = character())
cluster_order = c("Cluster11",
                  "Cluster3",
                  "Cluster5",
                  "Cluster1",
                  "Cluster8",
                  "Cluster7",
                  "Cluster6",
                  "Cluster10",
                  "Cluster4",
                  "Cluster9",
                  "Cluster2")
cluster_order = gsub("$","_bulk",cluster_order)
colnames(fData(rawctxp5)) = gsub("_num_genes_.*","_bulk",colnames(fData(rawctxp5)))
fData(rawctxp5) = fData(rawctxp5)[,cluster_order]
plot_heatmap(rawctxp5, genes = markers, 
             type = "bulk", cluster_by = "row", 
             cluster_type = "kmeans",
             k=4)
ggsave("results/DEgenes_markers/Astro/markers.pdf", h=10, w=5)
plot_tsne_metadata(rawctxp5, 
                   color_by = "astro_Cluster", 
                   xcol = "astro_x", ycol = "astro_y",
                   facet_by = "astro_Cluster")
rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "celltype")
plot_heatmap(rawctxp5, genes = markers, 
             type = "bulk", cluster_by = "both")
### Label subtypes
rawctxp5$subtype[rawctxp5$astro_Cluster%in%c("Cluster11","Cluster3")] = "Astrocytes_2"
rawctxp5$subtype[rawctxp5$astro_Cluster=="Cluster2"] = "doublets"
rawctxp5$subtype[rawctxp5$astro_Cluster%in%c("Cluster1","Cluster4","Cluster5","Cluster6","Cluster7","Cluster8","Cluster9","Cluster10")] = "Astrocytes_1"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "astro_x", ycol = "astro_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")