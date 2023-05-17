# Requirements ----
# Seurat 4.3.0
# scran 1.24.1
# glmGamPoi 1.8.0
# tradeSeq 1.10.0
# scSorter 0.0.2
# ggplot2 3.4.0
# drosophila2.db 3.13.0
# sctransform 0.3.5 (indirectly by Seurat)
# slingshot 2.4.0
# Data: Move processed data files into directory structure matching 10X Genomics
# filtered matrices:
# H3/filtered_feature_bc_matrix/barcodes.tsv.gz
# H3/filtered_feature_bc_matrix/features.tsv.gz
# H3/filtered_feature_bc_matrix/matrix.mtx.gz
# T3A/filtered_feature_bc_matrix/barcodes.tsv.gz
# T3A/filtered_feature_bc_matrix/features.tsv.gz
# T3A/filtered_feature_bc_matrix/matrix.mtx.gz
library(AnnotationDbi)
library(dplyr)
library(drosophila2.db)
library(glmGamPoi)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(grDevices)
library(openssl)
library(scales)
library(scran)
library(scSorter)
library(Seurat)
library(SingleCellExperiment)
library(slingshot)
library(tradeSeq)

# FlyBase gene symbol lookup ----
MatchFlybaseOldNames <- function(gene.names) {
  # If present (Dmel 6.47 out-of-date genome), then use the name 'esg' for this
  # marker gene.
  gene.names[grep('FBgn0287768', gene.names)] = 'esg'
  stopifnot(length(unique(gene.names)) == length(gene.names))
  gene.names
}

# mapIds helper, only inserting unambiguous and non-NA symbol lookups into the
# key name array (10X features). The 10X features were originally FlyBase key
# (FBgn####), and we want to update as many as possible to readable symbols.
# This also applies this to the Read10X SingleCellExperiment.
ReadFlyBase10X <- function(dirname, db = drosophila2.db, keyname = 'FLYBASE') {
  sce = Read10X(dirname)
  key.name <- rownames(sce)
  alias.name <- mapIds(db, keys=rownames(sce), column='SYMBOL', keytype=keyname)
  alias.name <- alias.name[!is.na(alias.name) & !duplicated(alias.name)]
  updated.name <- key.name
  updated.name[match(names(alias.name), updated.name)] = alias.name
  rownames(sce) <- MatchFlybaseOldNames(updated.name)
  sce
}

# 10X data intake and preprocessing ----
h3 = (ReadFlyBase10X("H3/filtered_feature_bc_matrix")
      %>% CreateSeuratObject %>% RenameCells(add.cell.id = "1")
      %>% SCTransform(vst.flavor = "v2"))
t3a = (ReadFlyBase10X("T3A/filtered_feature_bc_matrix")
       %>% CreateSeuratObject %>% RenameCells(add.cell.id = "2")
       %>% SCTransform(vst.flavor = "v2"))

# Selecting min unique features with an unbiased threshold. We can select by
# percentile. Generally, the cluster labeled as ISC/EB later (esg+) is fairly
# comparable, but may have different levels. Therefore, after initial
# clustering, we a quantile in that cluster and apply threshold to all cells.
filter.seq.depth.by.cluster = function(seur) {
  # Identify the esg+ cluster (there may be multiple esg+ cluster with very low
  # read depth, and one good esg+ cluster).
  esg.cluster = Idents(seur)[
    which.max(AverageExpression(seur, assay = 'SCT', features = 'esg')$SCT)
  ]
  seur = subset(
    seur,
    nFeature_RNA >= (
      quantile(subset(seur, idents = esg.cluster)$nFeature_RNA, 0.5)
    ))
}

# Initial clustering, num principal components selected by ElbowPlot.
# Judiciously, this is similar between H3 and T3A and we can avoid bias from
# selecting a different num pcs.
h3 = (
  h3 %>% RunPCA(verbose = F)
  %>% FindNeighbors(reduction = "pca", dims = 1:8, verbose = F)
  %>% FindClusters(resolution = 0.25) %>% filter.seq.depth.by.cluster)
t3a = (
  t3a %>% RunPCA(verbose = F)
  %>% FindNeighbors(reduction = "pca", dims = 1:8, verbose = F)
  %>% FindClusters(resolution = 0.25) %>% filter.seq.depth.by.cluster)
# Clustering after quality control threshold applied. The appearance of the PCA
# has changed on the smaller population of cells.
h3 = (h3 %>% RunPCA(verbose = F)
      %>% RunUMAP(reduction = "pca", dims = 1:6, verbose = F)
      %>% FindNeighbors(reduction = "pca", dims = 1:6, verbose = F)
      %>% FindClusters(resolution = 0.25))
t3a = (t3a %>% RunPCA(verbose = F)
       %>% RunUMAP(reduction = "pca", dims = 1:6, verbose = F)
       %>% FindNeighbors(reduction = "pca", dims = 1:6, verbose = F)
       %>% FindClusters(resolution = 0.5))
# Have both UMAPs place the esg+ (stem cell) population on the left hand side.
# This required t3a to be flipped horizontally (dim number 1).
t3a@reductions$umap@cell.embeddings[, 1] = -t3a@reductions$umap@cell.embeddings[, 1]

# Size factors (needed for regression) ----
sf.merge.full = calculateSumFactors(
  GetAssayData(merge(h3, t3a), assay = 'RNA', slot = 'counts'),
  # Separate analysis of Seurat clusters is not needed, but separate analysis of
  # batches as the "clusters" parameter is useful.
  clusters = c(rep(1, ncol(h3)), rep(2, ncol(t3a))))
names(sf.merge.full) = Cells(merge(h3, t3a))
h3$scran = sf.merge.full[Cells(h3)]
t3a$scran = sf.merge.full[Cells(t3a)]

# Call each cluster as a particular cell population ----
# Labels applied to Seurat clusters on the basis of marker genes (Supplemental
# Figure 5).
h3.called = matrix(
  c(
    '0', 'ISC/EB', # esg+
    '1', 'aEC', # betaTry+, CG12374+
    '2', 'ISC/EB', # esg+
    '3', 'ISC/EB', # esg+
    '4', 'EE', # IA-2+, AstC+
    '5', 'LFC', # mag+
    '6', 'mEC', # Vha16-1+, weakly MtnC+
    '7', 'pEC', # LManVI+
    '8', 'pEC', # Gs2+
    '9', 'aEC', # betaTry+, CG12374+
    '10', 'aEC' # betaTry+, Npc2f+
  ),
  ncol=2,
  byrow=T)

t3a.called = matrix(
  c(
    '0', 'aEC', # betaTry+, CG12374+
    '1', 'ISC/EB', # esg+
    '2', 'ISC/EB', # esg+
    '3', 'EE', # IA-2+, AstC+
    '4', 'ISC/EB',
    '5', 'ISC/EB',
    '6', 'aEC', # betaTry+, CG12374+
    '7', 'LFC', # mag+
    '8', 'pEC', # LManVI+, Gs2+
    '9', 'mEC', # Vha16-1+
    '10', 'unknown', # mixed expression pattern and unknown marker IM33
    '11', 'pEC' # Gs2+
  ),
  ncol=2,
  byrow=T)

# Levels and color scale to be used in figures.
call.levels = c(
  'ISC/EB',
  'aEC',
  'mEC',
  'pEC',
  'LFC',
  'EE',
  'unknown')
# Omission of one color from the palette is due to an initial level name
# (EC-like) which later was removed.
call.col = setNames(
  c(hue_pal()(7)[c(1, 3:7)], '#cccccc'),
  call.levels)

# After applying levels, automatically remove the 'unknown' label (not used)
h3$call = factor(
  factor(h3.called[match(levels(Idents(h3)), h3.called[,1]), 2][as.numeric(Idents(h3))],
  levels = call.levels))
t3a$call = factor(
  t3a.called[match(levels(Idents(t3a)), t3a.called[,1]), 2][as.numeric(Idents(t3a))],
  levels = call.levels)

# UMAP plot (Figure 6). ----
umap.title = element_text(size=32)
umap.legend = element_text(size=32)
umap.axis = element_text(size=28)
umap.axis.title = element_text(size=28)
h3.umap = DimPlot(h3, group.by = 'call', cols = call.col[levels(h3$call)]) + ggtitle('Cell populations (H3)') + theme(axis.title=umap.axis.title, axis.text=umap.axis, legend.text=umap.legend, title=umap.title) + guides(color=F)
t3a.umap = DimPlot(t3a, group.by = 'call', cols = call.col[levels(t3a$call)]) + ggtitle('Cell populations (H3T3A)') + theme(axis.title=umap.axis.title, axis.text=umap.axis, legend.text=umap.legend, title=umap.title)
ggsave(file='Cell-Population-UMAP.pdf', plot=ggarrange(h3.umap, t3a.umap, ncol=2, common.legend = T, legend = 'right'),
       width=20, height=8)

# Cell type enrichment; printable cell type % table. ----
ec_count = function(seur) sum(seur$call %in% c('aEC', 'mEC', 'pEC'))
ee_count = function(seur) sum(seur$call == 'EE')
esg_count = function(seur) sum(seur$call == 'ISC/EB')
total_count = function(seur) sum(seur$call != 'unknown')
cell.types = sapply(
  list(H3=h3, H3T3A=t3a),
  function(seur) list(
    EC=ec_count(seur) / total_count(seur),
    EE=ee_count(seur) / total_count(seur),
    ISCEB = esg_count(seur) / total_count(seur)))

# Seurat-SCTransform integration of batches ----
integ = list(H3 = subset(h3, call == 'ISC/EB'), H3T3A = subset(t3a, call == 'ISC/EB'))
integ[['H3']] = integ[['H3']] %>% FindVariableFeatures(nfeatures = 5000) %>% RunPCA(npcs = 10, verbose = F)
integ[['H3T3A']] = integ[['H3T3A']] %>% FindVariableFeatures(nfeatures = 5000) %>% RunPCA(npcs = 10, verbose = F)
Idents(integ[['H3']]) = 'ISC/EB'
Idents(integ[['H3T3A']]) = 'ISC/EB'
integ.features = intersect(VariableFeatures(integ[['H3']]), VariableFeatures(integ[['H3T3A']]))
integ = integ %>% PrepSCTIntegration(anchor.features = integ.features)
integ.anchors = integ %>% FindIntegrationAnchors(normalization.method = "SCT", anchor.features = integ.features)
seur = integ.anchors %>% IntegrateData(normalization.method = "SCT")
seur = seur %>% RunPCA(npcs = 10, verbose = F) %>% RunUMAP(reduction = 'pca', dims = 1:6)
seur@meta.data$condition = factor(NA, levels=c('H3', 'H3T3A'))
seur$condition[match(Cells(integ[['H3']]), Cells(seur))] = 'H3'
seur$condition[match(Cells(integ[['H3T3A']]), Cells(seur))] = 'H3T3A'

# Slingshot pseudotime ----
# Try 6 dim (used above for UMAP): reduced performance to yield DE on H3
# wild-type, along with H3T3A mutant. Therefore, we retained the default 50 dim
dim.pca.pt = 50
sce.h3 = slingshot(
  as.SingleCellExperiment(subset(h3, call == 'ISC/EB'), assay = 'SCT'),
  reducedDim = h3@reductions$pca@cell.embeddings[h3$call == 'ISC/EB', 1:dim.pca.pt])
sce.t3a = slingshot(as.SingleCellExperiment(subset(t3a, call == 'ISC/EB'), assay = 'SCT'))
sce.t3a = slingshot(
  as.SingleCellExperiment(subset(t3a, call == 'ISC/EB'), assay = 'SCT'),
  reducedDim = t3a@reductions$pca@cell.embeddings[t3a$call == 'ISC/EB', 1:dim.pca.pt])
# T3A was found to be reversed compared to H3.
sce.t3a$slingPseudotime_1 = max(sce.t3a$slingPseudotime_1) - sce.t3a$slingPseudotime_1

# Heatmap DE gene selection ----
# The nonlinear tradeSeq model is hugely expensive. First, fit a linear model to
# pseudotime separately for H3 and H3T3A, and retain genes where either sample
# has p.adj < 0.01.
gene.keep.1 = intersect(
  FindVariableFeatures(subset(t3a, call == 'ISC/EB'), assay = 'SCT', nfeatures = 5000) %>% VariableFeatures,
  FindVariableFeatures(subset(h3, call == 'ISC/EB'), assay = 'SCT', nfeatures = 5000) %>% VariableFeatures)
gene.keep.1 = intersect(VariableFeatures(h3), VariableFeatures(t3a))
glm.heatmap.signif = glm_gp(
  cbind(
    GetAssayData(subset(h3, call == 'ISC/EB'), assay = 'RNA', slot = 'counts'),
    GetAssayData(subset(t3a, call == 'ISC/EB'), assay = 'RNA', slot = 'counts'))[
      gene.keep.1,
    ],
  ~ condition * pt,
  data.frame(
    condition = rep(c('H3', 'H3T3A'), c(sum(h3$call == 'ISC/EB'), sum(t3a$call == 'ISC/EB'))),
    pt = c(sce.h3$slingPseudotime_1, sce.t3a$slingPseudotime_1),
    row.names = c(subset(colnames(h3), h3$call == 'ISC/EB'), subset(colnames(t3a), t3a$call == 'ISC/EB'))),
  size_factors=subset(
    c(h3$scran, t3a$scran),
    c(h3$call, t3a$call) == 'ISC/EB'),
  on_disk = F)
colnames(glm.heatmap.signif$Beta)
glm.heatmap.signif.de = as.data.frame(
  cbind(
    subset(test_de(glm.heatmap.signif, c(0, 0, 1, 0)), select = c('pval', 'lfc')),
    subset(test_de(glm.heatmap.signif, c(0, 0, 1, 1)), select = c('pval', 'lfc'))
  ),
  row.names = rownames(glm.heatmap.signif$Beta))
gene.keep.2 = subset(
  gene.keep.1,
  glm.heatmap.signif.de[,1] < 0.05 & glm.heatmap.signif.de[,3] < 0.05
)
sce.h3 = sce.h3 %>% fitGAM(genes = gene.keep.2)
sce.t3a = sce.t3a %>% fitGAM(genes = gene.keep.2)
slingshot.assoc = list(H3 = sce.h3, H3T3A = sce.t3a) %>% lapply(
  function(sce) associationTest(sce, contrastType = 'consecutive', l2fc=0.5)[
    match(gene.keep.2, rownames(sce)),
  ]
)

# Use p = 0.01 associationTest (not p adjusted). gene.keep.3 - actual list of
# genes for Supplemental Figure 6.
gene.keep.3 = subset(
  gene.keep.2,
  apply(as.data.frame(lapply(slingshot.assoc, function(df) df$pvalue)), 1, mean) < 0.01)
# DTW using normal filter ----
seur = RunPCA(seur, assay = 'integrated', reduction.name = 'pca.integrated', quiet=T)
pseudotime.smooth = function(pt, reducedDim, bandwidth = NULL) {
  if (is.null(bandwidth)) bandwidth = sd(pt)
  for (pc in colnames(reducedDim)) {
    pca.fit = ksmooth(pt, reducedDim[, pc], k = 'normal', b = bandwidth)
    reducedDim[, pc] = approx(pca.fit, xout = pt)$y
  }
  reducedDim
}
library(dtw, quiet=T)
seur.dtw = dtw(
  subset(seur, condition == 'H3')@reductions$pca.integrated@cell.embeddings[order(sce.h3$slingPseudotime_1), 1:dim.pca.pt] %>% pseudotime.smooth(sce.h3$slingPseudotime_1),
  subset(seur, condition == 'H3T3A')@reductions$pca.integrated@cell.embeddings[order(sce.t3a$slingPseudotime_1), 1:dim.pca.pt] %>% pseudotime.smooth(sce.t3a$slingPseudotime_1))

# Heatmap gene log-levels data, per-cell (Supplemental Figure 6). ----
# Remove extreme cells where the density in pseudotime is very low.
sce.h3 = sce.h3[, order(sce.h3$slingPseudotime_1)]
sce.t3a = sce.t3a[, order(sce.t3a$slingPseudotime_1)]
# Boolean to select which x locations to graph. TRUE - plot every x location in
# pseudotime, even the extrema of pseudotime which were too noisy.
seur.dtw.heatmap = rep(TRUE, length(seur.dtw$index1))
heatmap.h3.cells = colnames(sce.h3)[seur.dtw$index1[seur.dtw.heatmap]]
heatmap.t3a.cells = colnames(sce.t3a)[seur.dtw$index2[seur.dtw.heatmap]]
heatmap.expr = GetAssayData(seur, assay = 'SCT', slot = 'data')[gene.keep.3, unique(c(heatmap.h3.cells, heatmap.t3a.cells))]
for (cells in list(Cells(h3), Cells(t3a))) {
  cells = intersect(cells, colnames(heatmap.expr))
  heatmap.expr[, cells] = t(scale(t(heatmap.expr[, cells]), center=F))
}

# Heatmap data - Pearson correlation distance
seur.dist = cor(t(GetAssayData(seur, assay = 'SCT', slot = 'scale.data')[gene.keep.3,]))
seur.dist = as.dist(1 - seur.dist)
seur.hclust = hclust(seur.dist)
seur.dend = as.dendrogram(seur.hclust)

heatmap.data = array(
  c(as.numeric(heatmap.expr[gene.keep.3, heatmap.h3.cells]),
    as.numeric(heatmap.expr[gene.keep.3, heatmap.t3a.cells])),
  dim = c(length(gene.keep.3), length(heatmap.h3.cells), 2),
  dimnames = list(gene.keep.3, 1:length(heatmap.h3.cells), c('H3', 'H3T3A')))

# Kernel smoother using pseudotime variable as the x axis.
heatmap.data.smooth = heatmap.data
heatmap.bandwidth = 1
heatmap.data.smooth[,, 'H3'] = t(apply(
  heatmap.data.smooth[,, 'H3'],
  1,
  function(v)
    ksmooth(
      sce.h3$slingPseudotime_1[seur.dtw$index1[seur.dtw.heatmap]],
      v,
      kernel = 'normal',
      bandwidth = heatmap.bandwidth,
      x.points = sce.h3$slingPseudotime_1[seur.dtw$index1[seur.dtw.heatmap]])$y))
heatmap.data.smooth[,, 'H3T3A'] = t(apply(
  heatmap.data.smooth[,, 'H3T3A'],
  1,
  function(v)
    ksmooth(
      sce.t3a$slingPseudotime_1[seur.dtw$index2[seur.dtw.heatmap]],
      v,
      kernel = 'normal',
      bandwidth = heatmap.bandwidth,
      x.points = sce.t3a$slingPseudotime_1[seur.dtw$index2[seur.dtw.heatmap]])$y))

# Try swapping dendrogram nodes based on the weighted mean of gene expression.
# Moment 0 is the magnitude of the gene levels, summed across the bar of graph.
# Moment 1 is the weighted mean of locations across the bar on the graph (1:n).
# This allows us to sum within a dendrogram, and compare the weighted mean (peak
# expression - early late) for a collection of genes (hierarchical clustering).
heatmap.data.moment.0 = apply(heatmap.data, c(1, 3), function(v) sum(exp(v)))
heatmap.data.moment.1 = apply(heatmap.data, c(1, 3), function(v) exp(v) %*% 1:length(v))
dend.identity = function(node) {
  if (typeof(node) == 'list') rapply(node, as.integer)
  else as.integer(node)
}
seur.dend = seur.dend %>% dendrapply(
  function(node) {
    if (typeof(node) == 'list') {
      nodes.left = dend.identity(node[[1]])
      mean.left = sum(heatmap.data.moment.1[nodes.left, 'H3']) / sum(heatmap.data.moment.0[nodes.left, 'H3'])
      nodes.right = dend.identity(node[[2]])
      mean.right = sum(heatmap.data.moment.1[nodes.right, 'H3']) / sum(heatmap.data.moment.0[nodes.right, 'H3'])
      if (mean.left < mean.right) {
        subtree.left = node[[1]]
        node[[1]] = node[[2]]
        node[[2]] = subtree.left
      }
      node
    }
    else node
  })
seur.gene.order = rank(apply(heatmap.data.moment.1, 1, sum) / apply(heatmap.data.moment.0, 1, sum))
seur.dend = reorder(seur.dend, seur.gene.order,
                    function(lookups) {
                      genes = match(lookups, seur.gene.order)
                      sum(heatmap.data.moment.1[genes, ]) / sum(heatmap.data.moment.0[genes, ])
                    })
heatmap.order = order.dendrogram(seur.dend)

# Apply upper threshold to genes with high expression, we only need to graph
# log-levels with an upper cap of 3.
level.threshold = 3
heatmap.table.h3 = within(
  expand.grid(x = 1:length(heatmap.h3.cells), y = heatmap.order),
  level <- pmin(level.threshold, heatmap.data.smooth[cbind(y, x, rep(1, length(x)))]))
# heatmap.table.h3$y = rep(seur.hclust$labels[heatmap.order], rep(length(heatmap.h3.cells), length(heatmap.order)))
heatmap.table.h3$y = rep(1:length(heatmap.order), rep(length(heatmap.h3.cells), length(heatmap.order)))
heatmap.table.h3$x = (heatmap.table.h3$x - 1) / (length(heatmap.h3.cells) - 1)

heatmap.table.h3t3a = within(
  expand.grid(x = 1:length(heatmap.h3.cells), y = heatmap.order),
  level <- pmin(level.threshold, heatmap.data.smooth[cbind(y, x, rep(2, length(x)))]))
heatmap.table.h3t3a$y = rep(1:length(heatmap.order), rep(length(heatmap.h3.cells), length(heatmap.order)))
heatmap.table.h3t3a$x = heatmap.table.h3$x

# Display Heatmap
heatmap.theme = scale_fill_viridis_c(option = "plasma", limits = c(0, level.threshold))
heatmap.x = label_percent()(seq(0, 1, length.out=length(heatmap.h3.cells)))
heatmap.x.breaks = seq(0,1, length.out = length(heatmap.h3.cells))[seq(1, length(heatmap.h3.cells), length.out = 6)[1:5]]
print.genes = seur.hclust$labels[heatmap.order]
print.genes[grep('lncRNA:alphagamma-element:CR32865', print.genes)] = 'CR32865'
heatmap.plot = (ggplot(heatmap.table.h3, aes(x, y, fill=level))
                + ggtitle("H3 wild-type")
                + geom_raster()
                + heatmap.theme
                + scale_x_continuous(labels=percent, breaks=heatmap.x.breaks, expand=c(0,0))
                + scale_y_continuous(labels=NULL, breaks=NULL, expand=c(0,0))
                + xlab(NULL)
                + ylab(NULL)
                + theme(legend.position = "none", axis.text = element_text(size = 6))
              | ggplot(heatmap.table.h3t3a, aes(x, y, fill=level))
                + geom_raster()
                + ggtitle("H3T3A mutant")
                + scale_x_continuous(labels=percent, breaks=heatmap.x.breaks, expand=c(0,0))
                + scale_y_continuous(labels=print.genes, breaks = 1:length(heatmap.order), position = 'right', expand=c(0,0))
                + xlab(NULL)
                + ylab(NULL)
                + heatmap.theme
                + theme(axis.text = element_text(size = 6)))

# Color plot of Dl and klu levels (detail) (Figure 6). ----
background.cell = col2rgb('lightgray')
background.cell.plot = rgb(background.cell[1]/255., background.cell[2]/255., background.cell[3]/255.)
klu.cell = col2rgb(hsv(0.4))
klu.cell.plot = rgb(klu.cell[1]/255., klu.cell[2]/255., klu.cell[3]/255.)
dl.cell = col2rgb(hsv(0.75))
dl.cell.plot = rgb(dl.cell[1]/255., dl.cell[2]/255., dl.cell[3]/255.)
mixed.cell = (klu.cell+dl.cell)/2.
mixed.cell.plot = rgb(mixed.cell[1]/255., mixed.cell[2]/255., mixed.cell[3]/255.)

color.cell = function(klu, dl) {
  # Color comes from weighted mean of klu.cell, dl.cell, and remaining share.
  background.contribution = pmax(0, 1 - klu - dl)
  color.mult = Matrix::Diagonal(x=1. / (background.contribution + klu + dl))
  color.table = ((background.cell %*% background.contribution
                  + klu.cell %*% klu + dl.cell %*% dl)
                 %*% color.mult) / 255.
  rgb(color.table['red',], color.table['green',], color.table['blue',])
}

GeneNormalizedLevel <- function(assay, gene) {
  gene.levels = as.numeric(GetAssayData(assay, 'data')[gene, ])
  gene.levels / max(gene.levels)
}
h3$DlNorm = GeneNormalizedLevel(GetAssay(h3, 'SCT'), 'Dl')
h3$KluNorm = GeneNormalizedLevel(GetAssay(h3, 'SCT'), 'klu')
t3a$DlNorm = GeneNormalizedLevel(GetAssay(t3a, 'SCT'), 'Dl')
t3a$KluNorm = GeneNormalizedLevel(GetAssay(t3a, 'SCT'), 'klu')

h3.marker.data = data.frame(
  x = h3@reductions$umap@cell.embeddings[,1],
  y = h3@reductions$umap@cell.embeddings[,2],
  color = color.cell(h3$KluNorm / max(h3$KluNorm), h3$DlNorm / max(h3$DlNorm)))
h3.marker.data = subset(h3.marker.data, h3$call == 'ISC/EB')

t3a.marker.data = data.frame(
  x = t3a@reductions$umap@cell.embeddings[,1],
  y = t3a@reductions$umap@cell.embeddings[,2],
  Dl = seq(0, 1, length.out = ncol(t3a)),
  klu = seq(0, 1, length.out = ncol(t3a)),
  color = color.cell(t3a$KluNorm / max(t3a$KluNorm), t3a$DlNorm / max(t3a$DlNorm)))
t3a.marker.data = subset(t3a.marker.data, t3a$call == 'ISC/EB')

# UMAP for H3 will be transposed due to the presentation having a clear
# orientation to it.
marker.adjacent = (
  ggplot(h3.marker.data, aes(x = y, y = x, colour = color)) + xlim(-6,1.5)+ylim(-10,0) + geom_point(size = 1) + scale_colour_identity() + ggtitle('H3 (ISC/EB population)')
  + labs(y = 'UMAP_1 (rotated)', x = 'UMAP_2 (rotated)')
  + theme(axis.text = element_blank())
  |
    ggplot(t3a.marker.data, aes(x, y = (7-2.5)-y, colour = color)) + xlim(-11, -4)+ylim(-2.5, 7) + geom_point(size = 1) + scale_colour_identity() + ggtitle('H3T3A (ISC/EB pop.)')
  + labs(x = 'UMAP_1 (rotated)', y = 'UMAP_2 (rotated)')
  + theme(axis.text = element_blank())
  # Now insert custom scales into right-hand graph
  +
    new_scale_colour() +
    geom_point(aes(colour = Dl), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = dl.cell.plot, labels=NULL, guide = guide_colorbar(barheight=3))
  +
    new_scale_colour() +
    geom_point(aes(colour = klu), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = klu.cell.plot, labels=NULL, guide = guide_colorbar(barheight=3))
  +
    new_scale_colour() +
    geom_point(aes(colour = Dl), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = color.cell(1,1), labels=NULL, guide = guide_colorbar(title='Dl+klu+', barheight=3))
)
marker.adjacent

marker.title = element_text(size=16)
marker.legend = element_text(size=16)
marker.axis.title = element_text(size=14)
marker.adjacent = ggarrange(
  ggplot(h3.marker.data, aes(x = y, y = x, colour = color)) + xlim(-6,1.5)+ylim(-10,0) + geom_point(size = 1) + scale_colour_identity() + ggtitle('H3 (esg+ cluster)')
  + labs(y = 'UMAP_1 (rotated)', x = 'UMAP_2 (rotated)')
  + theme(axis.text = element_blank(), title = marker.title, legend.text = marker.legend, axis.title = marker.axis.title),
  ggplot(t3a.marker.data, aes(x, y = (7-2.5)-y, colour = color)) + xlim(-11, -4)+ylim(-2.5, 7) + geom_point(size = 1) + scale_colour_identity() + ggtitle('H3T3A (esg+ cluster)')
  + labs(x = 'UMAP_1 (rotated)', y = 'UMAP_2 (rotated)')
  + theme(axis.text = element_blank(), title = marker.title, legend.text = marker.legend, axis.title = marker.axis.title)
  # Now insert custom scales into right-hand graph
  +
    new_scale_colour() +
    geom_point(aes(colour = Dl), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = dl.cell.plot, labels=NULL, guide = guide_colorbar(barheight=5))
  +
    new_scale_colour() +
    geom_point(aes(colour = klu), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = klu.cell.plot, labels=NULL, guide = guide_colorbar(barheight=5))
  +
    new_scale_colour() +
    geom_point(aes(colour = Dl), shape = NA) +
    scale_colour_gradient(low = background.cell.plot, high = color.cell(1,1), labels=NULL, guide = guide_colorbar(title='Dl+klu+', barheight=5)),
  ncol=1,
  nrow=2,
  common.legend = T,
  legend = 'right'
)

# scSorter - used for outlier removal (Unknown class). ----
sctype = scSorter(
  GetAssayData(seur), data.frame(Type=c('ISC','EB'), Marker=c('Dl','klu')))
names(sctype$Pred_Type) = colnames(seur)

# GLM regression analysis of DE in pseudotime-aligned samples. ----
# Create a new SCE with repeated cells (dtw index). This allows us to adjust the
# observations (two observations at the same location, from two samples).
align.sce = function(sce, dtw.lookup) {
  sce = sce[, order(sce$slingPseudotime_1)]
  sce.cells = paste0(colnames(sce)[dtw.lookup], "_DTW", 1:length(dtw.lookup))
  sce = sce[, dtw.lookup]
  colnames(sce) = sce.cells
  sce
}
pt.glm.features = intersect(rownames(GetAssay(h3, 'SCT')), rownames(GetAssay(t3a, 'SCT')))
# Create a new SingleCellExperiment, repeating some cells to align in
# pseudotime, then the column at index i in each SCE may be considered as a
# paired observation (and the independent variable may be updated to be their
# mean - pt.combined pseudotime).
sce.h3.pt = as.SingleCellExperiment(h3[pt.glm.features, colnames(sce.h3)], assay = 'RNA')
sce.h3.pt$slingPseudotime_1 = sce.h3$slingPseudotime_1
sce.h3.pt$scsorter = factor(
  sctype$Pred_Type[colnames(sce.h3.pt)],
  levels = c('ISC', 'EB', 'Unknown'))
sce.h3.pt$is.labeled = (sce.h3.pt$scsorter != 'Unknown') %>% replace(is.na(.), F)
sce.h3.pt = align.sce(sce.h3.pt, seur.dtw$index1)
sce.t3a.pt = as.SingleCellExperiment(t3a[pt.glm.features, colnames(sce.t3a)], assay = 'RNA')
sce.t3a.pt$slingPseudotime_1 = sce.t3a$slingPseudotime_1
sce.t3a.pt$scsorter = factor(
  sctype$Pred_Type[colnames(sce.t3a.pt)],
  levels = c('ISC', 'EB', 'Unknown'))
sce.t3a.pt$is.labeled = (sce.t3a.pt$scsorter != 'Unknown') %>% replace(is.na(.), F)
sce.t3a.pt = align.sce(sce.t3a.pt, seur.dtw$index2)
pt.combined = (
  sce.h3.pt$slingPseudotime_1 + sce.t3a.pt$slingPseudotime_1
) / 2.
pt.combined = scale(pt.combined)
sce.h3.pt$pt.combined = pt.combined
sce.t3a.pt$pt.combined = pt.combined
sce.h3.pt$condition = 'H3'
sce.t3a.pt$condition = 'H3T3A'

# GLM of linear change of gene levels in pseudotime (for Dl and klu)
pt.linear.glm = glm_gp(
  cbind(assay(sce.h3.pt, 'counts'), assay(sce.t3a.pt, 'counts'))[
    ,
    c(sce.h3.pt$is.labeled, sce.t3a.pt$is.labeled)],
  ~ condition * pt,
  subset(
    data.frame(
      condition = factor(c(rep('H3', ncol(sce.h3.pt)), rep('H3T3A', ncol(sce.t3a.pt)))),
      pt = rep(pt.combined, 2),
      row.names = c(colnames(sce.h3.pt),colnames(sce.t3a.pt))),
    c(sce.h3.pt$scsorter != 'Unknown', sce.t3a.pt$scsorter != 'Unknown')),
  size_factors=subset(
    c(sce.h3.pt$scran, sce.t3a.pt$scran),
    c(sce.h3.pt$scsorter != 'Unknown', sce.t3a.pt$scsorter != 'Unknown')),
  on_disk = F)
# Outliers GLM - for regression in Supplemental Figure 5 retaining all cells.
pt.linear.glm.outliers = glm_gp(
  cbind(assay(sce.h3.pt, 'counts'), assay(sce.t3a.pt, 'counts'))[
    intersect(
      h3@assays$SCT@var.features,
      t3a@assays$SCT@var.features),
  ],
  ~ condition * pt,
  data.frame(
    condition = factor(c(rep('H3', ncol(sce.h3.pt)), rep('H3T3A', ncol(sce.t3a.pt)))),
    pt = rep(pt.combined, 2),
    row.names = c(colnames(sce.h3.pt),colnames(sce.t3a.pt))),
  size_factors=c(sce.h3.pt$scran, sce.t3a.pt$scran),
  # overdispersion_shrinkage = F,
  on_disk = F)

# glmGamPoi p-value calculation
# coef: Intercept conditionH3T3A         pt conditionH3T3A:pt
# Our ln-fold change is conditionH3T3A:pt (difference in slope of pt lines).
pt.linear.glm.de = test_de(pt.linear.glm, c(0,0,0,1))
pt.linear.glm.h3.de = test_de(pt.linear.glm, c(0,0,1,0))
pt.linear.glm.h3t3a.de = test_de(pt.linear.glm, c(0,0,1,1))

pt.linear.glm.de.outliers = test_de(pt.linear.glm.outliers, c(0,0,0,1))
pt.linear.glm.h3.de.outliers = test_de(pt.linear.glm.outliers, c(0,0,1,0))
pt.linear.glm.h3t3a.de.outliers = test_de(pt.linear.glm.outliers, c(0,0,1,1))

# Data frame for ggplot dot plot.
plot_expr_data = function(gene) {
  condition = paste0(
    c(sce.h3.pt$condition, sce.t3a.pt$condition),
    ifelse(c(sce.h3.pt$scsorter, sce.h3.pt$scsorter) == 'Unknown', '(Outlier)', ''))
  df = data.frame(
    levels = log1p(LOG1P_SCALE * c(assay(sce.h3.pt, 'counts')[gene, ] / sce.h3.pt$scran, assay(sce.t3a.pt, 'counts')[gene, ] / sce.t3a.pt$scran)),
    x = rep(pt.combined, 2),
    condition = condition
  )
  barcode_length = as.numeric(gregexpr('_', rownames(df))[[1]])
  df = df[
    sapply(
      split(seq_along(rownames(df)), substr(rownames(df), 1, barcode_length)),
      function(unique_cells) unique_cells[1]),
  ]
  # shuffle df deterministically, so that the x-axis (y = 0) does not obscure
  # one of the conditions by covering with the other sample.
  df[order(openssl::md5(rownames(df))), ]
}
# Data frame for log-linear line from regression model (plotted as log1p).
plot_fit_data = function(intercept, slope) {
  xmin = -2
  xmax = 2
  xs = seq(xmin, xmax, length.out = 100)
  data.frame(
    x = xs,
    y = log1p(exp(intercept + slope * xs))
  )
}

# Add a multiplier before log1p, if needed (scoot the points up along the y-axis
# so that they are not distorted by being close to 0). This was not needed.
LOG1P_SCALE = 1
gene_breaks = c(0, 1, 2, 5, 10, 20, 50)
# Significance line (annotation) and label.
plot_significance = function(plot, x, ymin, ymax, yshift = 0, label) {
  label.shift = ifelse(x > 0, -0.2, 0.2)
  (plot + geom_line(
    data = data.frame(
      x = x,
      y = c(log1p(LOG1P_SCALE * ymin), log1p(LOG1P_SCALE * ymax))
    ),
    aes(x, y),
    col = 'black'
  ) + annotate(
    'text',
    x = x + label.shift,
    y = mean(c(log1p(LOG1P_SCALE * ymin), log1p(LOG1P_SCALE * ymax))) + yshift,
    label = label,
    size = glm.sig.text.size
  ))
}
# Supplemental Outlier Plot (Supplemental Figure 5). ----
outlier_plot = function(g, sample) {
  data = subset(plot_expr_data(g), condition %in% c(sample, paste0(sample, '(Outlier)')))
  (ggplot(data, aes(x, levels, col = condition))
    + geom_point(shape = 1, size = 2)
    + scale_color_manual(
      values = c(
        hsv(0.65, 0.6, 0.75, 0.5),
        hsv(0.5, 1, 1, 0.5)
      )
    )
    # Override legend: Filled square
    + guides(color = guide_legend(override.aes = list(shape = 15)))
    + xlim(-2,2)
    + scale_y_continuous(breaks = log1p(LOG1P_SCALE * gene_breaks), labels = gene_breaks)
    + ylab('Normalized Count')
    + xlab('Pseudotime (sd from mean pseudotime)')
    + ggtitle(paste0(g, ' levels in ', sample, ' condition')))
}
glm.pt.regression.outlier.identity = ggarrange(
  outlier_plot('Dl', 'H3'),
  outlier_plot('Dl', 'H3T3A'),
  outlier_plot('klu', 'H3'),
  outlier_plot('klu', 'H3T3A'),
  nrow=2,
  ncol=2
)
# Regression Plot (Figure 6), actual log-linear pseudotime plot. ----
glm.title = element_text(size=16)
glm.legend = element_text(size=16)
glm.axis = element_text(size=14)
glm.axis.title = element_text(size=14)
glm.sig.text.size = 16/.pt
# Plot outliers, or remove them?
# glm.regression.plot.condition = setNames(c('H3','H3T3A','H3','H3T3A'), c('H3','H3T3A','H3(Outlier)','H3T3A(Outlier)'))
glm.regression.plot.condition = c(H3='H3', H3T3A='H3T3A')
glm.regression.plot.beta = pt.linear.glm$Beta
glm.regression.plot.de = pt.linear.glm.de
glm_regression_plot = function(g, sig_x, sig_ymin, sig_ymax, sig_yshift, sig_text, sig_format) {
  expr_data = plot_expr_data(g)
  expr_data = subset(expr_data, condition %in% names(glm.regression.plot.condition))
  expr_data$condition = glm.regression.plot.condition[expr_data$condition]
  ((ggplot(
    expr_data,
    aes(x, levels, col = condition))
    + geom_point(shape = 1, size = 2)
    + geom_line(data = plot_fit_data(log(LOG1P_SCALE)+glm.regression.plot.beta[g,'Intercept'], glm.regression.plot.beta[g,'pt']),
                aes(x, y, col = 'H3'), size=1.5)
    + geom_line(data = plot_fit_data(log(LOG1P_SCALE)+glm.regression.plot.beta[g,'Intercept']+glm.regression.plot.beta[g,'conditionH3T3A'], glm.regression.plot.beta[g,'pt']+glm.regression.plot.beta[g,'conditionH3T3A:pt']), aes(x, y, col = 'H3T3A'),
                size=1.5)
    + scale_color_manual(
      values = c(
        hsv(0.55, 0.65, 0.75, 0.4),
        hsv(0.05, 0.5, 0.9, 0.6)
      )
    )
    + xlim(-2,2)
    + scale_y_continuous(breaks = log1p(LOG1P_SCALE * gene_breaks), labels = gene_breaks)
    + ylab('Normalized Count')
    + xlab('Pseudotime (sd from mean pseudotime)')
    + ggtitle(paste0(g, ' levels in DTW-aligned pseudotime'),
              subtitle = sprintf(sig_format, glm.regression.plot.de$pval[match(g, glm.regression.plot.de$name)]))
    + theme(title = glm.title, legend.text = glm.legend, axis.text = glm.axis, axis.title = glm.axis.title)
  )
    %>% plot_significance(sig_x, sig_ymin, sig_ymax, sig_yshift, sig_text))
}

glm.pt.regression = ggarrange(
  glm_regression_plot('Dl', -2, 5.8, 8, sig_yshift=0.25, sig_text='n.s.', sig_format = 'p = %.04f'),
  glm_regression_plot('klu', 2, 1.4, 4, sig_yshift=0, '****', sig_format = 'p = %.04g'),
  ncol=1,
  common.legend = T,
  legend = 'bottom')
