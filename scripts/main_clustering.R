#' ---
#' title: Single cell analysis of WT and FX mice P5 cortex
#' author: Elisa Donnard
#' date: December 8th, 2018
#' output:
#'    html_document:
#'      toc: true
#'      highlight: zenburn
#' ---


#' ## Our pipeline

# most of the packages we use have been wrapped in easier to use functions in our lab analysis package
# devtools::install_github("garber-lab/SignallingSingleCell")
library(SignallingSingleCell)


#' ## Input

# raw_counts for all samples
load("merged_data/filt_full_results_umis.RData")

#' ## Initial clustering

# generate expression set 
ctxp5 = construct_ex_sc(raw_counts)

# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(ctxp5, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(ctxp5)),  
                        cutoff = 0.7)

# Reduce dimensions with PCA and tSNE in the top PCs that explain 90% of the variance
ctxp5 = dim_reduce(ctxp5, 
                   genelist = vargenes, 
                   pre_reduce = "vPCA", 
                   iterations = 500, 
                   nVar = 0.9)

# Define initial clusters by density peak
ctxp5 = cluster_sc(ctxp5, dimension = "2d", method = "density", s = 1)
plot_tsne_metadata(ctxp5, color_by = "Cluster", facet_by = "Cluster")
ctxp5 = calc_agg_bulk(ctxp5, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

# Microglia
plot_tsne_gene(ctxp5, "Aif1", log_scale = T) # known marker
plot_tsne_gene(ctxp5, "Csf1r", log_scale = T) # known marker
pop = "Cluster7"
ctxp5$celltype[ctxp5$Cluster%in%pop] = 1

# Pericytes and Endothelial
plot_tsne_gene(ctxp5, "Cldn5", log_scale = T) # endothelial specific (known)
plot_tsne_gene(ctxp5, "Vtn", log_scale = T) # pericyte specific (known)
plot_tsne_gene(ctxp5, "Ets1", log_scale = T) # endothelial cell migration / angiogenesis
plot_tsne_gene(ctxp5, "Apcdd1", log_scale = T) 
pop = c("Cluster4","Cluster8","Cluster15")
ctxp5$celltype[ctxp5$Cluster%in%pop] = 2

# Astrocytes
plot_tsne_gene(ctxp5, "Aldoc", log_scale = T) # Astro
plot_tsne_gene(ctxp5, "Slc1a3", log_scale = T) # Astro
plot_tsne_gene(ctxp5, "Dbi", log_scale = T) # Ependymal too
pop = c("Cluster2","Cluster6","Cluster13")
ctxp5$celltype[ctxp5$Cluster%in%pop] = 3

# Oligodendrocytes
plot_tsne_gene(ctxp5, "Sox10", log_scale = T)
pop = "Cluster3"
ctxp5$celltype[ctxp5$Cluster%in%pop] = 4

# Neuronal
plot_tsne_gene(ctxp5, "Sox11", log_scale = T)
plot_tsne_gene(ctxp5, "Bcl11b", log_scale = T)
plot_tsne_gene(ctxp5, "Mapt", log_scale = T)
plot_tsne_gene(ctxp5, "Dlx1", log_scale = T)
plot_tsne_gene(ctxp5, "Gad2", log_scale = T)
plot_tsne_gene(ctxp5, "Cdca7", log_scale = T)
plot_tsne_gene(ctxp5, "Reln", log_scale = T)
pop = c("Cluster5","Cluster10","Cluster11")
ctxp5$celltype[ctxp5$Cluster%in%pop] = 5

# Remove erithrocytes and empty droplets/doublets
plot_tsne_gene(ctxp5, "Hba-a1", log_scale = T)
plot_tsne_gene(ctxp5, "Alas2", log_scale = T)
plot_tsne_gene(ctxp5, "Rplp1", log_scale = T)
pop = c("Cluster1","Cluster9","Cluster12","Cluster14","Cluster16")
resultsctxp5 = ctxp5[,!ctxp5$Cluster%in%pop]
#remove leaky erythrocyte genes
rmgenes = c("Hba-a1","Hba-a2","Hbb-b1","Hbb-bh1","Hbb-bh2","Hbb-bs","Hbb-bt","Hbb-y","Hbq1a","Hbq1b")
resultsctxp5 = resultsctxp5[!rownames(resultsctxp5)%in%rmgenes,]

#Add metadata 
resultsctxp5$library = stringr::str_split_fixed(colnames(resultsctxp5),"_",2)[,1]
resultsctxp5$genotype = substr(resultsctxp5$library, 1, 2)
resultsctxp5$sample= substr(resultsctxp5$library, 1, 3)
resultsctxp5$beads = "V2"
resultsctxp5$beads[grep("WT2|WT3|FX3", resultsctxp5$sample)] = "V3"
resultsctxp5$lib_size = colSums(exprs(resultsctxp5))
saveRDS(resultsctxp5, "results/raw_exsc.rds")

#' ## Normalize UMI counts

# genes to use for estimating size factor batch V2
genes = rowMeans(exprs(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                                  resultsctxp5$beads=="V2"]))
meancut = quantile(genes, 0.8)
clusterlist = pData(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                               resultsctxp5$beads=="V2"])$celltype
# calculate size factors batch V2
sceV2 = SingleCellExperiment::SingleCellExperiment(list(
  counts = exprs(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                            resultsctxp5$beads=="V2"])))
sceV2 = scran::computeSumFactors(sceV2, 
                                 sizes=c(20, 40, 60, 80, 100), 
                                 positive=T, 
                                 clusters=clusterlist, 
                                 min.mean=meancut)
# genes to use for estimating size factor batch V3
genes = rowMeans(exprs(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                                  resultsctxp5$beads=="V3"]))
meancut = quantile(genes, 0.8)
clusterlist = pData(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                               resultsctxp5$beads=="V3"])$celltype
# calculate size factors batch V3
sceV3 = SingleCellExperiment::SingleCellExperiment(list(
  counts = exprs(resultsctxp5[rowSums(exprs(resultsctxp5))>0,
                            resultsctxp5$beads=="V3"])))
sceV3 = scran::computeSumFactors(sceV3, 
                                 sizes=c(20, 40, 60, 80, 100), 
                                 positive=T, 
                                 clusters=clusterlist, 
                                 min.mean=meancut)
# merge size factors
sf_cells = as.data.frame(c(sizeFactors(sceV2),
                           sizeFactors(sceV3)))
rownames(sf_cells) = c(colnames(sceV2),
                       colnames(sceV3))
colnames(sf_cells) = "sizefactor"
save(sf_cells, file="results/scran_size_factors.RData")

## limit the size factors to 0.1*gMean > SF < 10*gMean (100X range):
## compute the geometric mean of size factors
sfs <- sf_cells[sf_cells$sizefactor>0,]
gMean <- exp(mean(log(sfs)))
nLt <- length(sf_cells[sf_cells$sizefactor<(0.1*gMean),])
nGt <- length(sf_cells[sf_cells$sizefactor>(10.0*gMean),])
# remove cells with size factor <0.1*gMean or >10*gMean
keepcells = rownames(sf_cells[sf_cells$sizefactor>(0.1*gMean) & 
                                sf_cells$sizefactor<(10*gMean),,
                              drop=F])
# add size factor to pData
density_clusters_raw = merge(pData(resultsctxp5), as.data.frame(sf_cells), by=0)
rownames(density_clusters_raw) = density_clusters_raw[,1]
density_clusters_raw = density_clusters_raw[,-1]
pData(resultsctxp5) = density_clusters_raw[colnames(resultsctxp5),]
saveRDS(resultsctxp5, file = "results/raw_exsc.rds")

# normalize counts and filter out bad cells
norm_counts = sweep(exprs(resultsctxp5),
                    2,
                    pData(resultsctxp5)[colnames(resultsctxp5),"sizefactor"],
                    FUN='/')
norm_counts[!is.finite(norm_counts)] = 0
norm_counts = norm_counts[,keepcells]
norm_counts = as.matrix(norm_counts[rowSums(norm_counts)>0,
                                    colSums(norm_counts)>0])
dim(norm_counts)
normctxp5 = construct_ex_sc(norm_counts)
pData(normctxp5) = pData(resultsctxp5)[colnames(normctxp5),]
saveRDS(normctxp5, "results/norm_exsc.rds")

#' ## Batch effect correction and reclustering

normctxp5 = readRDS("results/norm_exsc.rds")
# select genes by minimum expression in 1% of cells and highest cv
vargenes = subset_genes(normctxp5, 
                        method = "CV", 
                        threshold = 3, 
                        minCells = (0.01*ncol(normctxp5)),  
                        cutoff = 0.8)


# Reduce dimensions with PCA and tSNE in the top PCs that explain 90% of the variance
normctxp5 = dim_reduce(normctxp5, 
                   genelist = vargenes, 
                   pre_reduce = "vPCA", 
                   iterations = 500, 
                   nComp = 50,
                   nVar = 0.9)

# Define initial clusters by density peak
normctxp5 = cluster_sc(normctxp5, dimension = "2d", method = "density", s = 1)
plot_tsne_metadata(normctxp5, color_by = "Cluster")
plot_tsne_metadata(normctxp5, color_by = "beads")

# run mnn correction
mnn = batchelor::fastMNN(log1p(exprs(normctxp5[,normctxp5$beads=="V2"])),
                         log1p(exprs(normctxp5[,normctxp5$beads=="V3"])),
                         subset.row = vargenes,
                         d=12)

mnnctxp5 = normctxp5
pData(mnnctxp5) = pData(mnnctxp5)[,grep("iPC",colnames(pData(mnnctxp5)), invert = T)]
pcs = SingleCellExperiment::reducedDim(mnn, "corrected")
colnames(pcs) = colnames(pData(normctxp5))[grep("iPC", colnames(pData(normctxp5)))]
pData(mnnctxp5) = cbind(pData(mnnctxp5),pcs)
tSNE_result = Rtsne::Rtsne(pcs, dims = 2, perplexity = 30, theta = 0.5, check_duplicates = F, pca = F, max_iter = 500, verbose = T)
mnnctxp5$x = tSNE_result$Y[,1]
mnnctxp5$y = tSNE_result$Y[,2]

# Define clusters by density peak
mnnctxp5 = cluster_sc(mnnctxp5, dimension = "2d", method = "density", s = 1)
plot_tsne_metadata(mnnctxp5, color_by = "Cluster")
plot_tsne_metadata(mnnctxp5, color_by = "beads")
mnnctxp5 = calc_agg_bulk(mnnctxp5, aggregate_by = "Cluster")

#' ## Identify celltypes and results up

# Microglia
plot_tsne_gene(mnnctxp5, "Aif1", log_scale = T) # known marker
plot_tsne_gene(mnnctxp5, "Csf1r", log_scale = T) # known marker
pop = c("Cluster7","Cluster17")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Microglia"

# Pericytes 
plot_tsne_gene(mnnctxp5, "Vtn", log_scale = T) # endothelial specific (known)
pop = c("Cluster8")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Endothelial"
# Endothelial
plot_tsne_gene(mnnctxp5, "Cldn5", log_scale = T) # pericyte specific (known)
pop = c("Cluster6")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Pericytes"

# Astrocytes
plot_tsne_gene(mnnctxp5, "Aldoc", log_scale = T) # Astro
plot_tsne_gene(mnnctxp5, "Slc1a3", log_scale = T) # Astro
pop = c("Cluster3","Cluster1","Cluster9","Cluster12")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Astrocytes"

#Ependymal
plot_tsne_gene(mnnctxp5, "Dbi", log_scale = T)
pop = c("Cluster18")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Ependymal"

# Oligodendrocytes and OPCs
plot_tsne_gene(mnnctxp5, "Sox10", log_scale = T)
plot_tsne_gene(mnnctxp5, "Olig1", log_scale = T)
pop = c("Cluster14","Cluster16","Cluster2")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Oligodendrocytes"

# Neuronal
plot_tsne_gene(mnnctxp5, "Sox11", log_scale = T)
plot_tsne_gene(mnnctxp5, "Bcl11b", log_scale = T)
plot_tsne_gene(mnnctxp5, "Mapt", log_scale = T)
plot_tsne_gene(mnnctxp5, "Dlx1", log_scale = T)
plot_tsne_gene(mnnctxp5, "Gad2", log_scale = T)
plot_tsne_gene(mnnctxp5, "Cdca7", log_scale = T)
plot_tsne_gene(mnnctxp5, "Reln", log_scale = T)
plot_tsne_gene(mnnctxp5, "Stmn1", log_scale = T)
pop = c("Cluster10","Cluster11","Cluster13","Cluster15","Cluster4","Cluster5")
mnnctxp5$celltype[mnnctxp5$Cluster%in%pop] = "Neurons"

saveRDS(mnnctxp5, file="results/mnn_exsc.rds")

#' ## Cell cycle classification per cell
library(org.Mm.eg.db)
rawctxp5 = readRDS("results/raw_exsc.rds")
set.seed(100)
sce = SingleCellExperiment::SingleCellExperiment(list(counts = as.matrix(exprs(rawctxp5[,colnames(mnnctxp5)]))))
mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
cyclegenes = AnnotationDbi::mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="SYMBOL", column="ENSEMBL")
assignments = scran::cyclone(sce, mm.pairs, gene.names=cyclegenes)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)
table(assignments$phases)
mnnctxp5$phases = assignments$phases
plot_tsne_metadata(mnnctxp5, color_by = "phases")
plot_tsne_gene(mnnctxp5[,mnnctxp5$phases=="G1"], gene = "H2afv", facet_by = "genotype")

#' ## Add the pdata from the corrected set back on the raw and normalized expression matrix
pData(normctxp5) = pData(mnnctxp5)[colnames(normctxp5),]
rawctxp5 = rawctxp5[,colnames(normctxp5)]
pData(rawctxp5) = pData(mnnctxp5)[colnames(rawctxp5),]
saveRDS(rawctxp5, "results/mnn_rawcounts_exsc.rds")
saveRDS(normctxp5, "results/mnn_normcounts_exsc.rds")

### Continue after running each individual subclustering script
normctxp5 = readRDS("results/subclusters_mnn_normcounts_exsc.rds")
normctxp5$subtype[normctxp5$celltype=="Ependymal"] = "Ependymal"
rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
pData(rawctxp5) = pData(normctxp5)[colnames(rawctxp5),]
saveRDS(normctxp5, "results/subclusters_mnn_normcounts_exsc.rds")
saveRDS(rawctxp5, "results/subclusters_mnn_rawcounts_exsc.rds")

#' ## Find marker genes in broad celltypes

rawctxp5 = readRDS("results/mnn_rawcounts_exsc.rds")
findDEmarkers(rawctxp5, 
            DEgroup="celltype",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            outdir="results/DEgenes_markers/", 
            minCells = 0.01)

#' ## Find DE genes in broad celltypes

findDEgenes(rawctxp5, 
            DEgroup="genotype",
            contrastID = "FX",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            outdir="results/DEgenes_genotype/batch_all_cells/", 
            facet_by="celltype",
            minCells = 0.15)

#' ## Find DE genes in subtypes

rawctxp5 = readRDS("results/subclusters_mnn_rawcounts_exsc.rds")
findDEgenes(rawctxp5, 
            DEgroup="genotype",
            contrastID = "FX",
            sizefactor = "sizefactor",
            lib_size = "lib_size",
            batchID="beads", 
            outdir="results/DEgenes_genotype/batch_subtype/", 
            facet_by="subtype",
            minCells = 0.15)
