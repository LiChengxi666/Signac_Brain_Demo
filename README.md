# Signac for Human Brain Single Cell Multiomics
### 数据集与工具简介
#### Signac 概述
Signac 是 2021 年由 Tim 等发表于 Nature Methods 上的一种用于单细胞染色质数据分析的框架（原文链接：https://pubmed.ncbi.nlm.nih.gov/34725479/ ），基于R语言搭建。Signac 作为一种独立的解决方案，专门用于单细胞染色质数据的分析，其一大特点为可以与多模态数据处理的 Seurat 软件包的无缝对接，以实现多模态单细胞数据集的分析。

Signac 主要实现的功能包含识别非细胞含量条码中的背景细胞、峰值调用、基因组区域的计数定量、细胞的质量控制筛选、降维、聚类、与单细胞基因表达数据的整合、交互式基因组浏览器风格的数据可视化、差异可达峰的识别、富集DNA序列基序的检测、转录因子结合位点的分析以及将峰与潜在调控目标基因关联等。此外，Signac还提供了一个框架，用于从单细胞DNA可及性实验中鉴定线粒体基因组变异，从而实现单细胞中克隆关系与DNA可及性联合分析。
#### Signac 基本工作流
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/1.png)
Seurat包采用Seurat对象作为其核心数据结构。Seurat对象由多个Assay对象组成，每个Assay对象包含单细胞的相关数据。Assay对象最初用于单细胞基因表达数据的分析，支持原始及处理后单细胞测量值和与每个特征相关的元数据的存储与检索。为了在Seurat框架内便于单细胞染色质数据的分析，Signac 加入了一种专用的 `ChromatinAssay` 对象类。`ChromatinAssay` 允许存储和检索分析单细胞染色质数据所需的信息，包括与每个实验特征相关的基因组范围、基因注释、基因组版本信息、DNA基序信息，以及作为tabix索引片段文件存储的单细胞数据。Signac 框架将在用 Seurat 读取单细胞测序数据的基础上，以fragments 形式读取 ATAC数据，之后进行质控、聚类、连接计算等一系列计算任务。
#### 数据集概述
本次实验分析了 10x Genomics 上公开的一份 Human Brain Single Cell Multiomics 数据集，具体连接点击[此处](https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-1-0-0)查看。该数据集为人类健康脑组织的 ATAC 和基因表达单细胞测序，均采用由 10x Genomics 处理好的输出文件（区别于原始的 `.bam` 或 `.fq`，等格式，具体如下：
* scATAC-seq: [fragments.tsv.gz](https://cf.10xgenomics.com/samples/cell-arc/1.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz)
* scRNA-seq: [filtered_feature_bc_matrix.tar.gz](https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.tar.gz)

其中 RNA 数据需通过 `tar -xzf filtered_feature_bc_matrix.tar.gz` 解压后使用，为 10x omics 标准格式，包含 `matrix.mtx`、`genes.tsv`、`barcodes.tsv` 三个文件。

此外，本实验还使用了标注好的细胞数据，来自 `SCENIC+` 文档，保存为 `cell_data.tsv`。最终数据文件夹结构如下：
```
.
├── cell_data.tsv
├── filtered_feature_bc_matrix
│   ├── barcodes.tsv.gz
│   ├── features.tsv.gz
│   └── matrix.mtx.gz
└── fragments.tsv.gz
```
### 分析流程
#### 软件安装
软件可以通过 R 直接安装，同时需按照 Seurat 等相关工具。
```
install.packages("Signac")
install.packages("Seurat")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("EnsDb.Hsapiens.v86")
```
安装后引入即可
```
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(readr)
```
#### 数据导入
RNAseq 数据的三个文件可以通过 Seurat `Read10X()` 函数直接导入，并创建 Seurat 对象。
```
# 导入数据（输入文件夹地址）
counts <- Read10X("data/filtered_feature_bc_matrix")

# 创建包含 RNAseq 的对象
brain <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)
```
ATACseq 数据可以在添加注释后通过 Signac 提供的函数导入，并以 ChromatinAssay 形式加入 Seurat 对象中。
```
# 从 hg38 获取注释
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotation) <- paste0('chr', seqlevels(annotation))

# 读取并将 ATAC assay 添加至对象
brain[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
```

#### 细胞注释
Signac 官方文档中的一篇教程给出了通过参考数据进行标注的方法，具体可参考[此处](https://stuartlab.org/signac/articles/pbmc_multiomic)。这里我们使用标注好的数据直接导入即可。
```
# 读取细胞注释文件
cell_metadata <- read_tsv("E:/brain/SCENIC+/test/data/cell_data.tsv", col_types = cols())

# 准备元数据
cell_meta_df <- data.frame(cell_metadata)
rownames(cell_meta_df) <- cell_meta_df[,1]  # 第一列作为行名（细胞ID）
cell_meta_df <- cell_meta_df[,-1]  # 移除第一列

# 确保 id 匹配
cell_ids_annotation <- rownames(cell_meta_df)
barcode_from_annotation <- gsub("-10x_multiome_brain$", "", cell_ids_annotation)

# 创建映射关系
mapping_df <- data.frame(
  original_id = cell_ids_annotation,
  barcode = barcode_from_annotation,
  stringsAsFactors = FALSE
)

# 找到共同的 barcode
brain_barcodes <- colnames(brain)
common_barcodes <- intersect(brain_barcodes, barcode_from_annotation)
print(paste("共同细胞数量:", length(common_barcodes)))

# 创建匹配的元数据
matched_metadata <- cell_meta_df[mapping_df$original_id[mapping_df$barcode %in% common_barcodes], ]
rownames(matched_metadata) <- mapping_df$barcode[mapping_df$barcode %in% common_barcodes]

# 只保留有注释的细胞
brain <- subset(brain, cells = common_barcodes)
matched_metadata <- matched_metadata[common_barcodes, ]

# 添加所有注释信息到Seurat对象
brain <- AddMetaData(brain, metadata = matched_metadata)

# 设置主要的细胞类型列
brain$cell_type <- brain$VSN_cell_type
Idents(brain) <- "VSN_cell_type"
```
通过 `print(unique(brain$VSN_cell_type))` 查看所有细胞类型，包含：
* 神经元
    * GC - Granule Cells (颗粒细胞)
    * PURK - Purkinje cells (浦肯野细胞)
    * GP - Globose/Golgi cells (球状/高尔基细胞)
* 抑制性中间神经元 (INH)
    * INH_VIP - VIP+ Interneurons (血管活性肠肽阳性中间神经元)
    * INH_SST - SST+ Interneurons (生长抑素阳性中间神经元)
    * INH_PVALB - Parvalbumin+ Interneurons (小白蛋白阳性中间神经元)
* 胶质细胞
    * MOL_A/MOL_B - Oligodendrocytes A/B (少突胶质细胞A型/B型)
    * OPC - Oligodendrocyte Precursor Cells (少突胶质细胞前体细胞)
    * AST - Astrocytes (星形胶质细胞)
    * AST_CER - Cerebellar Astrocytes (小脑星形胶质细胞)
    * ASTP - Astrocyte Progenitors (星形胶质细胞前体)
    * MGL - Microglia (小胶质细胞)
* 其他
    * ENDO - Endothelial cells (内皮细胞)
#### 质量控制
使用 DNA 可及性数据计算每个细胞的质量控制指标，并去除这些指标的异常细胞，以及 RNA 或 ATAC 检测中计数较低或异常高的细胞。
```
DefaultAssay(brain) <- "ATAC"

brain <- NucleosomeSignal(brain)
brain <- TSSEnrichment(brain)
```
之后通过密度散点图和提琴图可视化结果：
```
DensityScatter(brain, 
               x = 'nCount_ATAC', 
               y = 'TSS.enrichment', 
               log_x = TRUE, 
               quantiles = TRUE)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot01.png)
```
VlnPlot(
  object = brain,
  features = c("nCount_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
  ncol = 2,
  pt.size = 0
)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot.png)
之后筛选出符合条件的数据：
```
brain <- subset(
  x = brain,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1800 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)
brain
```
此外还可以查看 ATAC 片段长度分布以确保质量
```
FragmentHistogram(object = brain, group.by = "VSN_cell_type")
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot04.png)
#### 数据处理
使用 SCTransform 对基因表达数据进行归一化，并使用 PCA 进行降维。
```
DefaultAssay(brain) <- "RNA"
brain <- SCTransform(brain)
brain <- RunPCA(brain)
```
使用潜在语义索引（LSI）处理 ATAC 数据，包含词频逆文档频率 (TF-IDF) 归一化和奇异值分解 （SVD）。
```
DefaultAssay(brain) <- "ATAC"
brain <- FindTopFeatures(brain, min.cutoff = 5)
brain <- RunTFIDF(brain)
brain <- RunSVD(brain)
```
#### 联合 UMAP 可视化
使用Seurat v4中的加权最近邻方法，计算出一个代表基因表达和 DNA 可及性测量的联合邻域图。
```
# 构建联合邻域图
brain <- FindMultiModalNeighbors(
  object = brain,
  reduction.list = list("pca", "lsi"), 
  dims.list = list(1:50, 2:40),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

# 运行 umap 可视化
brain <- RunUMAP(
  object = brain,
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

DimPlot(brain, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot02.png)
#### 峰与基因链接计算
对于每个基因，通过计算基因表达与邻近峰的可及性之间的相关性，并校正由 GC 含量、总体可及性和峰大小引起的偏差，找到可能调控该基因的一组峰。

在整个基因组上运行此步骤可能非常耗时，因此这里选取了一些细胞系的代表基因进行运算。省略 `genes.use` 参数，即可使用相同的函数查找所有基因的链接：
```
DefaultAssay(brain) <- "ATAC"

# 计算GC含量
brain <- RegionStats(brain, genome = BSgenome.Hsapiens.UCSC.hg38)

# 选择代表基因
brain_genes <- c(
  "RBFOX3",    # 神经元标记 (NeuN)
  "SYP",       # 突触蛋白
  "MAP2",      # 神经元标记
  "GFAP",      # 星形胶质细胞
  "MBP",       # 少突胶质细胞
  "AIF1",      # 小胶质细胞 (Iba1)
  "NEFL",      # 神经丝蛋白
  "ENO2",      # 神经元特异性烯醇化酶
  "SLC1A2",    # 胶质细胞谷氨酸转运蛋白
  "PLP1"       # 髓鞘蛋白
)
available_brain_genes <- brain_genes[brain_genes %in% rownames(brain[["RNA"]])]
print("可用的脑组织标记基因:")
print(available_brain_genes)


# 计算基因和峰值的关联
brain <- LinkPeaks(
  object = brain,
  peak.assay = "ATAC",
  expression.assay = "SCT",
  genes.use = brain_genes
)
```
使用 `CoveragePlot()` 函数将这些链接可视化
```
# 使用全部细胞系绘图
idents.plot <- brain$VSN_cell_type

# brain <- SortIdents(brain)
DefaultAssay(brain) <- "ATAC"
p1 <- CoveragePlot(
  object = brain ,
  region = "MBP",
  features = "MBP",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 5000,
  extend.downstream = 5000,
  links = TRUE,  # 明确设置显示links
)

p2 <- CoveragePlot(
  object = brain,
  region = "RBFOX3",
  features = "RBFOX3",
  expression.assay = "SCT",
  idents = idents.plot,
  extend.upstream = 5000,
  extend.downstream = 5000
)

patchwork::wrap_plots(p1, p2, ncol = 1)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot03.png)
这里选取了`MBP`和`RBFOX3`两个基因，可以看到二者分别在与之相关的胶质细胞和神经元细胞中高表达，相关区域染色质可及性信号也较强。
#### 特定基因特征可视化
可以通过 `DotPlot()` 函数和 `FeaturePlot()` 函数以点状热图和降维散点图的形式可视化差异基因表达。
```
# 点状热图
DefaultAssay(brain) <- "SCT"
DotPlot(brain, 
        features = c("MBP", "RBFOX3", "GAD1", "SLC1A2"),
        group.by = "VSN_cell_type") + RotatedAxis()

# umap散点图
FeaturePlot(
  object = brain,
  features = brain_genes,
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 5
)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot05.png)
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot14.png)
#### 寻找差异化标记
除去人为选择，`Signac` 也提供了 `FindMarkers()` 函数，用于寻找所有基因的差异化标记。下面以比较神经元细胞和胶质细胞的 ATAC peaks 为例进行展示：
```
celltype1 <- c("GC", "GP", "PURK")
celltype2 <- c("MOL_A", "MOL_B", "OPC", "AST_CER", "ASTP", "AST", "MGL")   

da_peaks <- FindMarkers(
  object = brain,
  ident.1 = celltype1,
  ident.2 = celltype2, # 设置为 NULL 则为除 ident.1 外的所有细胞
  min.pct = 0.05,
  test.use = 'LR',  # 逻辑回归
  latent.vars = 'nCount_ATAC')  # 控制测序深度
  
  print("Differential accessibility analysis results:")
  print(paste("Total peaks tested:", nrow(da_peaks)))
  print(paste("Significant peaks (p_val_adj < 0.05):", sum(da_peaks$p_val_adj < 0.05)))
  head(da_peaks)

  # 选出 top10 的peaks
  top_da_peaks <- rownames(da_peaks)[1:min(10, nrow(da_peaks))]
```
之后可以用与前面类似的方法进行可视化：
```

DefaultAssay(brain) <- "ATAC"
peak_features <- rownames(brain)

# 小提琴图
plot1 <- VlnPlot(
  object = brain,
  features = top_da_peaks[1],  # 这里选取第一个peak为例子
  pt.size = 0.1,
  idents = target_idents,
  assay = "ATAC"
) 

# 特征图
plot2 <- FeaturePlot(
  object = brain,
  features = top_da_peaks[1],
  pt.size = 0.1,
  max.cutoff = 'q95',
  reduction = "umap"
) 

# 点状热图
plot3 <- DotPlot(brain, 
  features = top_da_peaks,
  group.by = "VSN_cell_type"
) + RotatedAxis()                 

# 组合图
combined_plot <- plot1 | plot2 | plot3
print(combined_plot)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot07.png)
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot09.png)
同时，使用 `CoveragePlot()` 也可以可视化这些peaks相关的区域信息。
```
CoveragePlot(
  object = brain,
  region = top_da_peaks[1],
  features = "HEY2", # 根据图中的注释选定
  expression.assay = "SCT",
  annotation = TRUE,
  peaks = TRUE,
  tile = TRUE,
  links = TRUE,
  extend.upstream = 10000,
  extend.downstream = 10000
)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot11.png)
相似的代码可以用于寻找基因表达的标记：
```
DefaultAssay(brain) <- "SCT" 

# RNA差异表达分析
de_genes <- FindMarkers(
    object = brain,
    ident.1 = celltype1,
    ident.2 = celltype2,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = 'wilcox'
)

print("Differential gene expression results:")
print(paste("Total genes tested:", nrow(de_genes)))
print(paste("Significant genes (p_val_adj < 0.05):", sum(de_genes$p_val_adj < 0.05)))
head(de_genes)
```
此外，`Signac` 还提供了 `FindAllMarkers()` 函数，用于寻找所有细胞类型共同的差异化标记：
```
all_markers <- FindAllMarkers(
    object = brain,
    only.pos = TRUE,
    min.pct = 0.05,
    logfc.threshold = 0.2,
    test.use = 'LR',
    latent.vars = 'nCount_ATAC'
)
```
这部分代码耗时较长，因此后两部分未实际运行。
#### Motif （基序）分析
`Signac` 提供了两种互补的基序分析方法：一种是在一组差异可及峰中查找高表达基序，另一种是在细胞组间进行差异基序活性分析。后者计算量较大，因此此处选择前者进行测试。
为了方便在 `Signac` 中进行基序分析，需要通过 `AddMotifs()` 创建一个 `Motif` 类来存储所有必需信息，包括位置权重矩阵 (PWM) 或位置频率矩阵 (PFM) 列表以及基序出现矩阵,并将其加入 Seurat 对象中：
```
library(JASPAR2020)
library(TFBSTools)

# 获取motif信息
pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(species = 9606, all_versions = FALSE))

# 添加motif信息到Seurat对象
brain <- AddMotifs(brain, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pwm_set)
```
为了识别潜在重要的细胞类型特异性调控序列，可以搜索在一组在不同细胞类型之间具有差异性的峰中过度表达的 DNA 基序，并通过超几何检验分析差异。这里采用上一部分中找到的差异峰。
```
# 找到全部符合要求的峰
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005 & da_peaks$pct.1 > 0.2, ])

# 富集分析
enriched.motifs <- FindMotifs(
  object = brain,
  features = top.da.peak
)

# 可视化
MotifPlot(
  object = brain,
  motifs = head(rownames(enriched.motifs))
)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot08.png)
#### 转录因子足迹分析
依据之前找到的基序，使用 `Footprint()` 可以计算绘制其对应的转录因子足迹信息，`PlotFootprint()` 函数可以进行可视化:
```
brain <- Footprint(
  object = brain,
  motif.name = c("TCFL5", "NRF1", "EGR2", "ZBTB14", "ZBTB33", "HINFP"),
  genome = BSgenome.Hsapiens.UCSC.hg38
)

# 对每种细胞都绘制转录因子足迹图
p2 <- PlotFootprint(brain, features = c("NRF1"))
p2 + patchwork::plot_layout(ncol = 1)
```
![](https://github.com/LiChengxi666/Signac_Brain_Demo/blob/main/plots/Rplot13.png)
### 总结
作为一个单细胞数据分析框架，Signac 的一大特点是其与通用的 Seurat 兼容性高，并且提供了较为丰富的功能。同时其环境支持性好，部署相对简单，可以作为分析 ATAC-seq 数据的参考工具。但是其相比更为常用的 ArchR，也存在运行速度慢、对原始数据支持较差的缺点
