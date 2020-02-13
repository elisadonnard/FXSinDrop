#' ## Batch effect correction and reclustering Neurons
library(SignallingSingleCell)
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
neurons = normctxp5[,normctxp5$celltype=="Neurons"]
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(neurons, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(neurons)),  
                        cutoff = 0.75)


# Reduce dimentions with PCA and tSNE in the top PCs that explain 90% of the variance
neurons = dim_reduce(neurons, 
                       genelist = vargenes, 
                       pre_reduce = "vPCA", 
                       iterations = 500, 
                       nComp = 50,
                       nVar = 0.93)

# Define initial clusters by density peak
neurons = cluster_sc(neurons, dimension = "2d", method = "density", s = 1)
plot_tsne_metadata(neurons, color_by = "Cluster")
plot_tsne_metadata(neurons, color_by = "beads")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(neurons[,neurons$beads=="V2"])),
                         log1p(exprs(neurons[,neurons$beads=="V3"])),
                         subset.row = vargenes,
                         d=20)

mnnneurons = neurons
pData(mnnneurons) = pData(mnnneurons)[,grep("iPC",colnames(pData(mnnneurons)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(neurons))[grep("iPC", colnames(pData(neurons)))]
pData(mnnneurons) = cbind(pData(mnnneurons),pcs)
tSNE_result = Rtsne::Rtsne(pcs, dims = 2, perplexity = 30, theta = 0.5, check_duplicates = F, pca = F, max_iter = 500, verbose = T)
mnnneurons$x = tSNE_result$Y[,1]
mnnneurons$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnneurons = cluster_sc(mnnneurons, dimension = "2d", method = "density", s = 2)
plot_tsne_metadata(mnnneurons, color_by = "Cluster")
plot_tsne_metadata(mnnneurons, color_by = "beads")
mnnneurons = calc_agg_bulk(mnnneurons, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

cellcut = min(table(mnnneurons$Cluster)/ncol(mnnneurons))
rawctxp5 = readRDS("results/mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            pd = pData(mnnneurons),
            DEgroup="Cluster",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            minCells = cellcut,
            outdir="results/DEgenes_markers/Neurons/")

colnames(pData(mnnneurons)) = gsub("^","neurons_",colnames(pData(mnnneurons)))
merged = merge(pData(rawctxp5), pData(mnnneurons), by=0, all=T)
rownames(merged) = merged[,1]
merged = merged[,-1]
pData(rawctxp5) = merged[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")

rawctxp5 = calc_agg_bulk(rawctxp5, aggregate_by = "neurons_Cluster")
markers = scan("results/DEgenes_markers/Neurons/markers", what = character())
cluster_order = c("Cluster6",
                  "Cluster1",
                  "Cluster11",
                  "Cluster5",
                  "Cluster7",
                  "Cluster13",
                  "Cluster2",
                  "Cluster12",
                  "Cluster14",
                  "Cluster3",
                  "Cluster9",
                  "Cluster10",
                  "Cluster8",
                  "Cluster15",
                  "Cluster4")
cluster_order = gsub("$","_bulk",cluster_order)
colnames(fData(rawctxp5)) = gsub("_num_genes_.*","_bulk",colnames(fData(rawctxp5)))
fData(rawctxp5) = fData(rawctxp5)[,cluster_order]
heat = plot_heatmap(rawctxp5, genes = markers, 
             type = "bulk", cluster_by = "row", 
             ceiling = 2.5, cluster_type = "kmeans",
             k=15, show_k = T)
reordered = c(names(heat[[2]]$cluster)[heat[[2]]$cluster==10],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==11],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==1],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==14],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==13],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==8],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==6],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==7],
              rev(names(heat[[2]]$cluster)[heat[[2]]$cluster==15]),
              names(heat[[2]]$cluster)[heat[[2]]$cluster==5],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==3],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==2],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==12],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==9],
              names(heat[[2]]$cluster)[heat[[2]]$cluster==4])
plot_heatmap(rawctxp5, genes = reordered, 
             type = "bulk", cluster_by = "none",
             gene_labels = name, 
             gene_labels_nudge = 1,
             ceiling = 2)
ggsave("results/DEgenes_markers/Neurons/markers.pdf",h=12,w=7)
plot_tsne_metadata(rawctxp5, 
                   color_by = "neurons_Cluster", 
                   xcol = "neurons_x", ycol = "neurons_y",
                   facet_by = "neurons_Cluster")
### Label subtypes
rawctxp5$subtype[rawctxp5$neurons_Cluster%in%c("Cluster4","Cluster15","Cluster8")] = "Ganglionic"
rawctxp5$subtype[rawctxp5$neurons_Cluster%in%c("Cluster3","Cluster9","Cluster10","Cluster14")] = "Interneuron"
rawctxp5$subtype[rawctxp5$neurons_Cluster=="Cluster12"] = "SVZ"
rawctxp5$subtype[rawctxp5$neurons_Cluster=="Cluster6"] = "CR"
rawctxp5$subtype[rawctxp5$neurons_Cluster%in%c("Cluster1","Cluster11","Cluster5","Cluster7","Cluster13","Cluster2")] = "Excitatory"
plot_tsne_metadata(rawctxp5, 
                   color_by = "subtype", 
                   xcol = "neurons_x", ycol = "neurons_y")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")
###
normctxp5 = readRDS("results/mnn_normcounts_exsc.rds")
pData(normctxp5) = pData(rawctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")

#' ## Markers for main subtypes

rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")
mainneurons = normctxp5[,normctxp5$subtype%in%c("Excitatory","Interneuron","Ganglionic")]
table(mainneurons$neurons_Cluster,mainneurons$subtype)
cellcut = min(table(mainneurons$subtype)/ncol(mainneurons))
findDEmarkers(rawctxp5, 
              pd = pData(mainneurons),
              DEgroup="subtype",
              sizefactor = "sizefactor",
              lib_size = "lib_size",
              batchID="beads", 
              minCells = cellcut,
              outdir="results/DEgenes_markers/Neurons/")

