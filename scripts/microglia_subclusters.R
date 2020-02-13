#' ## Batch effect correction and reclustering microglia
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
microglia = normctxp5[,normctxp5$celltype=="Microglia"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(microglia, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(microglia)),  
                        cutoff = 0.9)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
microglia = dim_reduce(microglia, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.88)

# Define initial clusters by density peak
microglia = cluster_sc(microglia, dimension = "2d", method = "density", s=2)
plot_tsne_metadata(microglia, color_by = "Cluster")
plot_tsne_metadata(microglia, color_by = "beads")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(microglia[,microglia$beads=="V2"])),
                         log1p(exprs(microglia[,microglia$beads=="V3"])),
                         subset.row = vargenes,
                         d=length(grep("iPC",colnames(pData(microglia)))))

mnnmicroglia = microglia
pData(mnnmicroglia) = pData(mnnmicroglia)[,grep("iPC",colnames(pData(mnnmicroglia)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(microglia))[grep("iPC", colnames(pData(microglia)))]
pData(mnnmicroglia) = cbind(pData(mnnmicroglia),pcs)
tSNE_result = Rtsne::Rtsne(pcs, 
                           dims = 2, 
                           perplexity = 15, 
                           theta = 0, 
                           check_duplicates = F, 
                           pca = F, 
                           max_iter = 500, 
                           verbose = T)
mnnmicroglia$x = tSNE_result$Y[,1]
mnnmicroglia$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnmicroglia = cluster_sc(mnnmicroglia, dimension = "2d", method = "density", num_clust = 11)
plot_tsne_metadata(mnnmicroglia, color_by = "Cluster")
plot_tsne_metadata(mnnmicroglia, color_by = "beads")
mnnmicroglia = calc_agg_bulk(mnnmicroglia, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

cellcut = min(table(mnnmicroglia$Cluster)/ncol(mnnmicroglia))
rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            pd = pData(mnnmicroglia),
            DEgroup="Cluster",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            minCells = cellcut,
            outdir="results/DEgenes_markers/Microglia/")

colnames(pData(mnnmicroglia)) = gsub("^","microglia_",colnames(pData(mnnmicroglia)))
merged = merge(pData(rawctxp5), pData(mnnmicroglia), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "microglia_Cluster")
markers = scan("results/DEgenes_markers/Microglia/results_markers", what = character())
cluster_order = c("Cluster1",
                  "Cluster11",
                  "Cluster2",
                  "Cluster6",
                  "Cluster4",
                  "Cluster5",
                  "Cluster10",
                  "Cluster8",
                  "Cluster3",
                  "Cluster7",
                  "Cluster9")
cluster_order = gsub("$","_bulk",cluster_order)
colnames(fData(rawctxp5)) = gsub("_num_genes_.*","_bulk",colnames(fData(rawctxp5)))
fData(rawctxp5) = fData(rawctxp5)[,cluster_order]
microglia_cluster = plot_heatmap(rawctxp5, 
             genes = markers, 
             type = "bulk", 
             cluster_type = "kmeans",
             k=4,
             cluster_by = "row",
             show_k = T)
ggsave("results/DEgenes_markers/Microglia/markers.pdf", h=10, w=5)
Ms4a_markers = names(microglia_cluster[[2]]$cluster[microglia_cluster[[2]]$cluster==3])
write.table(Ms4a_markers, "genelist/Ms4a_cluster_markers.txt", quote = F, row.names = F, col.names = F)
plot_tsne_metadata(rawctxp5, 
                   color_by = "microglia_Cluster", 
                   xcol = "microglia_x", ycol = "microglia_y",
                   facet_by = "microglia_Cluster")
### Label subtypes
rawctxp5$subtype[rawctxp5$microglia_Cluster=="Cluster1"] = "Ms4a7_high"
rawctxp5$subtype[rawctxp5$microglia_Cluster=="Cluster9"] = "Border_Macrophages"
rawctxp5$subtype[rawctxp5$microglia_Cluster=="Cluster7"] = "Tcells"
rawctxp5$subtype[rawctxp5$microglia_Cluster%in%c("Cluster2","Cluster3","Cluster4","Cluster5","Cluster6","Cluster8","Cluster10","Cluster11")] = "Microglia"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "microglia_x", ycol = "microglia_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")
