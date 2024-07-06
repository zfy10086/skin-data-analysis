library(Seurat)
library(GetoptLong)
library(EnhancedVolcano)
set_python_path("/public/home/zhangfy/anaconda3/envs/r/bin/python")

draw_path <- "draw/"

groups <- c(
    "DSQ-23W" = "Fetus",
    "DSQ-C1" = "Adult",
    "GSM5494441" = "Keloid"
)

fetus <- Load10X_Spatial(
    data.dir = "data/ST/DSQ-23W",
    filename = "filtered_feature_bc_matrix.h5",
    filter.matrix = FALSE
)
fetus$orig.ident <- "DSQ-23W"
fetus$group <- "Fetus"
f_features <- Features(fetus)

adult <- Load10X_Spatial(
    data.dir = "data/ST/DSQ-C1",
    filename = "filtered_feature_bc_matrix.h5",
    filter.matrix = FALSE
)
adult$orig.ident <- "DSQ-23W"
adult$group <- "Adults"
adult <- adults[f_features, ]

keloid <- Load10X_Spatial(
    data.dir = "data/ST/GSM5494441",
    filename = "filtered_feature_bc_matrix.h5",
    filter.matrix = FALSE
)
keloid$orig.ident <- "GSM5494441"
keloid$group <- "Keloid"
keloid <- keloid[f_features, ]


# int.skin <- list(fetus, adult, keloid)
# int.skin <- merge(x = int.skin[[1]], y = int.skin[-1], add.cell.ids = NULL)
# int.skin <- NormalizeData(int.skin)
# int.skin <- FindVariableFeatures(int.skin, selection.method = "vst", nfeatures = 2000)
# int.skin <- ScaleData(int.skin, features = NULL)
# int.skin <- RunPCA(int.skin, npcs = 50, verbose = FALSE)
# int.skin <- IntegrateLayers(
#     object = int.skin,
#     method = HarmonyIntegration,
#     orig.reduction = "pca",
#     new.reduction = "harmony",
#     npcs = 50, # 默认
#     verbose = FALSE
# )
# int.skin[["Spatial"]] <- JoinLayers(int.skin[["Spatial"]])

integrated_objs <- function(
    objslist,
    n_features = 2000,
    npcs = 50
){
   
    int.objs <- merge(x = objslist[[1]], y = objslist[-1], add.cell.ids = NULL)
    
    int.objs <- NormalizeData(int.objs)
    int.objs <- FindVariableFeatures(int.objs, selection.method = "vst", nfeatures = n_features)
    int.objs <- ScaleData(int.objs, features = NULL)

    int.objs <- RunPCA(int.objs[[1]], npcs = npcs, verbose = FALSE)

    int.objs <- IntegrateLayers(
        object = int.objs,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        npcs = npcs,
        verbose = FALSE
    )

    int.objs[["Spatial"]] <- JoinLayers(int.objs[["Spatial"]])

    return(int.objs)
}

run_analysis <- function(
    obj,
    assay = "Spatial",
    use.reduction = "harmony",
    npcs = 50,
    cluster_method = "leiden",
    resolution = 0.5
) {

    DefaultAssay(obj) <- assay

    obj <- RunUMAP(obj, reduction = use.reduction, dims = 1:npcs)
    obj <- RunTSNE(obj, reduction = use.reduction, dims = 1:npcs)

    obj <- FindNeighbors(obj, reduction = use.reduction, dims = 1:npcs)
    if (cluster_method == "louvain") cluster_method <- 1
    if (cluster_method == "leiden") cluster_method <- 4

    obj <- FindClusters(
        obj,
        resolution = resolution,
        algorithm = cluster_method,
        method = "igraph"
    )

    return(obj)
}

draw_analysis <- function(
    obj,
    split.by = "group",
    group.by = "group",
    markers = NULL,
    marker_split = FALSE,
    img_savedir = NULL,
    spotsize = NULL,
    show_labels = TRUE
) {

    if(!is.null(img_savedir)){
        if(!dir.exists(img_savedir)) dir.create(img_savedir)
    }

    p1 <- DimPlot(obj,
        reduction = "umap", split.by = split.by,
        label = show_labels, pt.size = spotsize
    )
    p2 <- DimPlot(obj,
        reduction = "tsne", split.by = split.by,
        label = show_labels, pt.size = spotsize
    )
    # 画一个左边为标记实验组与对照组, 右边为聚类结果的umap和tsne
    p3.1 <- DimPlot(obj,
        reduction = "umap", group.by = group.by,
        pt.size = spotsize
    )
    p3.2 <- DimPlot(obj,
        reduction = "umap", label = show_labels,
        repel = TRUE, pt.size = spotsize
    )
    p3 <- p3.1 + p3.2
    p4.1 <- DimPlot(obj,
        reduction = "tsne", group.by = group.by,
        pt.size = spotsize
    )
    p4.2 <- DimPlot(obj,
        reduction = "tsne", label = show_labels,
        repel = TRUE, pt.size = spotsize
    )
    p4 <- p4.1 + p4.2

    savename1 <- paste0(
        img_savedir, "/split_umap.jpg"
    )
    savename2 <- paste0(
        img_savedir, "/split_tsne.jpg"
    )
    savename3 <- paste0(
        img_savedir, "/group_umap&cluster.jpg"
    )
    savename4 <- paste0(
        img_savedir, "/group_tsne&cluster.jpg"
    )
    ggplot2::ggsave( # split umap
        savename1,
        plot = p1,
        width = 5.5 * 2,
        height = 5,
        dpi = 480
    )
    ggplot2::ggsave( # split tsne
        savename2,
        plot = p2,
        width = 5.5 * 2,
        height = 5,
        dpi = 480
    )
    ggplot2::ggsave( # group umap
        savename3,
        plot = p3,
        width = 13,
        height = 5,
        dpi = 480
    )
    ggplot2::ggsave( # group tsne
        savename4,
        plot = p4,
        width = 13,
        height = 5,
        dpi = 480
    )
}

int.skin <- list(fetus, adult, keloid)
int.skin <- integrated_objs(int.skin)
int.skin <- run_analysis(int.skin)

fetus <- subset(int.skin, subset = group == "Fetus")
fetus@project.name <- "DSQ-23W"
# fetus@meta.data <- fetus@meta.data[, -c(6, 7)]
fetus@images <- fetus@images["slice1"]
fetus <- run_analysis(fetus, resolution = 0.8, assay = "Spatial", use.reduction = "harmony", npcs = 50, cluster_method = "leiden")
draw_analysis(
    fetus,
    split.by = "orig.ident",
    group.by = NULL,
    img_savedir = paste0(draw_path, fetus@project.name)
)



draw_markers <- c("KRT1", "TP63", "DCN", "MRC1", "TRAC", "S100B", "ADGRL4", "MMRN1")

cols <- scale_fill_gradientn(colors = c("white", "#ffbd00", "#b38400"))
for(obj in c(fetus, adult, keloid)){
    for(i in draw_markers){
        p <- SpatialPlot(obj, features = i) + cols
        ggsave(plot = p, filename = qq("@{draw_path}/@{fetus@project.name}/@{i}.pdf"))
    }
}

# manual annotation
fetus_temp <- fetus
Idents(fetus_temp) <- "seurat_clusters"
fetus_temp <- RenameIdents(
    object = fetus_temp,  
    "1" = "Dermis", "2" = "Subcutaneous tissue",
    "3" = "Subcutaneous tissue",
    "4" = "Subcutaneous tissue", "5" = "Fibroblasts",
    "6" = "Fibroblasts", "7" = "Epidermis",
    "8" = "Hair follicle"
)
Idents(fetus_temp) <- factor(
    Idents(fetus_temp),
    levels = c("Epidermis", "Dermis", "Hair follicle", "Fibroblasts", "Subcutaneous tissue")
)

# p <- SpatialPlot(
#     fetus_temp,
#     label = TRUE,
#     # label.size = 12, # legend尺寸
#     stroke = 0
# )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{fetus@project.name}/spatial_FBs56_clst.pdf"
#     )
# )

# p <- DimPlot(fetus_temp,
#         reduction = "tsne",
#         label = TRUE, pt.size = NULL
#     )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{fetus@project.name}/spatial_FBs56_tsne.pdf"
#     ), width = 8, height = 7
# )

# ======
adult <- subset(int.skin, subset = group == "Adult")
adult@project.name <- "DSQ-C1"
# adult@meta.data <- adult@meta.data[, -c(6, 7)]
adult@images <- adult@images["slice1.DSQ.C1"]
adult <- run_analysis(adult, resolution = 0.8, assay = "Spatial", use.reduction = "harmony", npcs = 50, cluster_method = "leiden")
draw_analysis(
    adult,
    split.by = "orig.ident",
    group.by = "orig.ident",
    img_savedir = paste0(draw_path, adult@project.name)
)

# manual annotation
adult_temp <- adult
Idents(adult_temp) <- "seurat_clusters"
adult_temp <- RenameIdents(
    object = adult_temp,  
    "1" = "Fibroblasts", "2" = "Dermis 1", "3" = "Epidermis",
    "4" = "Fibroblasts", "5" = "Dermis 2", "6" = "Hair follicle",
    "7" = "Fibroblasts", "8" = "Dermis 3", "9" = "Dermis 4",
    "10" = "Dermis 5", "11" = "Dermis 6", "12" = "Dermis 6",
    "13" = "Dermis 6"
)
Idents(adult_temp) <- factor(
    Idents(adult_temp),
    levels = c(
        "Epidermis", "Dermis 1", "Dermis 2", "Dermis 3",
        "Dermis 4", "Dermis 5", "Dermis 6", "Hair follicle",
        "Fibroblasts"
    )
)

# p <- SpatialPlot(
#     adult_temp,
#     label = TRUE,
#     # label.size = 12, # legend尺寸
#     stroke = 0
#     )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{adult@project.name}/spatial_FBs147_clst.pdf"
#     )
# )
# p <- DimPlot(adult_temp,
#         reduction = "tsne",
#         label = TRUE, pt.size = NULL
#     )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{adult@project.name}/spatial_FBs147_tsne.pdf"
#     ), width = 7.5, height = 7
# )

# ======
keloid <- subset(int.skin, subset = group == "Keloid")
keloid@project.name <- "GSM5494441"
# keloid@meta.data <- keloid@meta.data[, -c(6, 7)]
keloid@images <- keloid@images["slice1.GSM5494441"]
keloid <- run_analysis(keloid, resolution = 0.5, assay = "Spatial", use.reduction = "harmony", npcs = 50, cluster_method = "leiden")
draw_analysis(
    keloid,
    split.by = "orig.ident",
    group.by = "orig.ident",
    img_savedir = paste0(draw_path, keloid@project.name)
)

# manual annotation

# keloid_temp <- keloid
# Idents(keloid_temp) <- "seurat_clusters"
# keloid_temp <- RenameIdents(
#     object = keloid_temp,  
#     "1" = "Dermis 1", "2" = "Fibroblasts", "3" = "Dermis 2",
#     "4" = "Dermis 3", "5" = "Dermis 4", "6" = "Epidermis"
# )
# Idents(keloid_temp) <- factor(
#     Idents(keloid_temp),
#     levels = c(
#         "Epidermis", "Dermis 1", "Dermis 2",
#         "Dermis 3", "Dermis 4", "Fibroblasts"
#     )
# )

# p <- SpatialPlot(
#     keloid_temp,
#     label = TRUE,
#     # label.size = 12, # legend尺寸
#     stroke = 0
#     )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{keloid@project.name}/spatial_FBs2_clst.pdf"
#     ) 
# )
# p <- DimPlot(keloid_temp,
#         reduction = "tsne",
#         label = TRUE, pt.size = NULL
#     )
# ggsave(plot = p, filename = qq(
#     "@{draw_path}/@{keloid@project.name}/spatial_FBs2_tsne.pdf"
#     ), width = 7.5, height = 7
# )

