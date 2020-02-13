#' ## Batch effect correction and reclustering oligo
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
oligo = normctxp5[,normctxp5$celltype=="Oligodendrocytes"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(oligo, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(oligo)),  
                        cutoff = 0.85)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
oligo = dim_reduce(oligo, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.87)

# Define initial clusters by density peak
oligo = cluster_sc(oligo, dimension = "2d", method = "density", s = 2)
plot_tsne_metadata(oligo, color_by = "Cluster")
plot_tsne_metadata(oligo, color_by = "beads")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(oligo[,oligo$beads=="V2"])),
                         log1p(exprs(oligo[,oligo$beads=="V3"])),
                         subset.row = vargenes,
                         d=length(grep("iPC",colnames(pData(oligo)))))

mnnoligo = oligo
pData(mnnoligo) = pData(mnnoligo)[,grep("iPC",colnames(pData(mnnoligo)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(oligo))[grep("iPC", colnames(pData(oligo)))]
pData(mnnoligo) = cbind(pData(mnnoligo),pcs)
tSNE_result = Rtsne::Rtsne(pcs, dims = 2, perplexity = 30, theta = 0.5, check_duplicates = F, pca = F, max_iter = 500, verbose = T)
mnnoligo$x = tSNE_result$Y[,1]
mnnoligo$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnoligo = cluster_sc(mnnoligo, dimension = "2d", method = "density", s = 2)
plot_tsne_metadata(mnnoligo, color_by = "Cluster")
plot_tsne_metadata(mnnoligo, color_by = "beads")
mnnoligo = calc_agg_bulk(mnnoligo, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

cellcut = min(table(mnnoligo$Cluster)/ncol(mnnoligo))
rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            pd = pData(mnnoligo),
            DEgroup="Cluster",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            minCells = cellcut,
            outdir="results/DEgenes_markers/Oligo/")

colnames(pData(mnnoligo)) = gsub("^","oligo_",colnames(pData(mnnoligo)))
merged = merge(pData(rawctxp5), pData(mnnoligo), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "oligo_Cluster")
markers = scan("results/DEgenes_markers/Oligo/results_markers", what = character())
markers = c(markers,"Btg2","Dll1","Dbx2")
cluster_order = c("Cluster6",
                  "Cluster5",
                  "Cluster11",
                  "Cluster9",
                  "Cluster1",
                  "Cluster4",
                  "Cluster10",
                  "Cluster12",
                  "Cluster2",
                  "Cluster7",
                  "Cluster8",
                  "Cluster3",
                  "Cluster13")
cluster_order = gsub("$","_bulk",cluster_order)
colnames(fData(rawctxp5)) = gsub("_num_genes_.*","_bulk",colnames(fData(rawctxp5)))
fData(rawctxp5) = fData(rawctxp5)[,cluster_order]
heat = plot_heatmap(rawctxp5, genes = markers, 
             type = "bulk", 
             cluster_by = "row", 
             ceiling = 3, 
             show_k = T,
             cluster_type = "kmeans",
             k=6)
reordered = c(names(heat[[2]]$cluster)[heat[[2]]$cluster==6],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==2],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==4],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==1],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==3],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==5])
plot_heatmap(rawctxp5,
             genes = reordered, 
             type = "bulk", 
             gene_labels = name, 
             gene_labels_nudge = 1,
             cluster_by = "F", 
             ceiling = 3)
ggsave("results/DEgenes_markers/Oligo/markers.pdf",h=12,w=6)
plot_tsne_metadata(rawctxp5, 
                   color_by = "oligo_Cluster", 
                   xcol = "oligo_x", ycol = "oligo_y",
                   facet_by = "oligo_Cluster")
### Label subtypes
rawctxp5$subtype[rawctxp5$oligo_Cluster%in%c("Cluster13","Cluster3","Cluster8")] = "Oligodendrocytes"
rawctxp5$subtype[rawctxp5$oligo_Cluster=="Cluster6"] = "NSC"
rawctxp5$subtype[rawctxp5$oligo_Cluster%in%c("Cluster1","Cluster2","Cluster4","Cluster5","Cluster7","Cluster9","Cluster10","Cluster11","Cluster12")] = "OPC_1"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "oligo_x", ycol = "oligo_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")