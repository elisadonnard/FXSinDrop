#' ## Batch effect correction and reclustering pericytes
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
peri = normctxp5[,normctxp5$celltype=="Pericytes"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(peri, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(peri)),  
                        cutoff = 0.85)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
peri = dim_reduce(peri, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.8)

# Define initial clusters by density peak
peri = cluster_sc(peri, dimension = "2d", method = "density", num_clust = 2)
plot_tsne_metadata(peri, color_by = "Cluster")
plot_tsne_metadata(peri, color_by = "beads")
plot_tsne_gene(peri,"Col1a1")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(peri[,peri$beads=="V2"])),
                         log1p(exprs(peri[,peri$beads=="V3"])),
                         subset.row = vargenes,
                         d=length(grep("iPC",colnames(pData(peri)))))

mnnperi = peri
pData(mnnperi) = pData(mnnperi)[,grep("iPC",colnames(pData(mnnperi)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(peri))[grep("iPC", colnames(pData(peri)))]
pData(mnnperi) = cbind(pData(mnnperi),pcs)
tSNE_result = Rtsne::Rtsne(pcs, dims = 2, 
                           perplexity = 50, 
                           theta = 0.5, 
                           check_duplicates = F, 
                           pca = F, max_iter = 500, 
                           verbose = T)
mnnperi$x = tSNE_result$Y[,1]
mnnperi$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnperi = cluster_sc(mnnperi, dimension = "2d", method = "density", num_clust = 3)
plot_tsne_metadata(mnnperi, color_by = "Cluster")
plot_tsne_metadata(mnnperi, color_by = "beads")
mnnperi = calc_agg_bulk(mnnperi, aggregate_by = "Cluster")
plot_tsne_gene(mnnperi,"Col1a1")
plot_tsne_gene(mnnperi,"Sfrp2")
plot_tsne_gene(mnnperi,"Kif1a")

#' ## Identify celltypes and results up

cellcut = min(table(mnnperi$Cluster)/ncol(mnnperi))
rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            pd = pData(mnnperi),
            DEgroup="Cluster",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            minCells = cellcut,
            outdir="results/DEgenes_markers/Pericytes/")

colnames(pData(mnnperi)) = gsub("^","peri_",colnames(pData(mnnperi)))
merged = merge(pData(rawctxp5), pData(mnnperi), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "peri_Cluster")
markers = scan("results/DEgenes_markers/Pericytes/results_markers", what = character())
cluster_order = c("Cluster1",
                  "Cluster3",
                  "Cluster2")
cluster_order = gsub("$","_bulk",cluster_order)
colnames(fData(rawctxp5)) = gsub("_num_genes_.*","_bulk",colnames(fData(rawctxp5)))
fData(rawctxp5) = fData(rawctxp5)[,cluster_order]
plot_heatmap(rawctxp5, genes = markers, 
             type = "bulk", cluster_by = "row", 
             ceiling = 3, cluster_type = "kmeans",
             k=6)
ggsave("results/DEgenes_markers/Pericytes/markers.pdf",h=8,w=4)
plot_tsne_metadata(rawctxp5, 
                   color_by = "peri_Cluster", 
                   xcol = "peri_x", ycol = "peri_y",
                   facet_by = "peri_Cluster")
plot_tsne_gene(rawctxp5,"Col1a1")

### Label subtypes
rawctxp5$subtype[rawctxp5$peri_Cluster=="Cluster1"] = "Fibroblasts"
rawctxp5$subtype[rawctxp5$peri_Cluster%in%c("Cluster2","Cluster3")] = "Pericytes"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "peri_x", ycol = "peri_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")
