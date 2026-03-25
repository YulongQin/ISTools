

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Convert between RDS and H5AD files
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#' Convert AnnData h5ad file to Seurat rds file
#'
#' Converts a Python AnnData h5ad file to a Seurat rds file, preserving spatial
#' coordinates, dimensional reductions, and optionally SCT-transformed data.
#'
#' @param file_path Character, path to the input h5ad file
#' @param python_path Character, path to Python executable with scanpy installed
#'        (default: "D:/APP/anaconda3/envs/scanpy/python.exe")
#' @param data_type Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
#'        (default: "stRNA")
#' @param convert_mode Character, conversion method - "scanpy" or "seurat"
#'        (default: "scanpy")
#' @param assay_id Character, assay name for the Seurat object (default: NULL,
#'        auto-detected based on data_type)
#' @param X_index Character, which matrix to use as counts - "X" or "rawX"
#'        (default: "rawX")
#' @param binsize Numeric, bin size for spatial coordinate scaling (default: 1)
#' @param SCT_index Logical, whether to convert SCT results from stereopy
#'        (default: FALSE)
#' @param reduction_index Logical, whether to preserve dimensional reductions
#'        (default: FALSE)
#' @param image_index Logical, whether to create spatial image object
#'        (default: NULL, auto-detected based on data_type)
#' @param return_object Logical, whether to return the Seurat object
#'        (default: TRUE)
#' @param grp_nm Character, group name for output organization (default: "sample1")
#' @param dir_nm Character, directory name for output (default: "M0_annh5ad2rds")
#'
#' @return If return_object = TRUE, returns a Seurat object; otherwise returns NULL
#'
#' @import reticulate
#' @import Matrix
#' @import rjson
#' @import Seurat
#' @importFrom anndata read_h5ad
#' @importFrom SeuratDisk Convert LoadH5Seurat
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert spatial transcriptomics h5ad to Seurat
#' seurat_obj <- annh5ad2rds(
#'   file_path = "data/spatial_data.h5ad",
#'   python_path = "~/anaconda3/envs/scanpy/bin/python",
#'   data_type = "stRNA",
#'   binsize = 100,
#'   reduction_index = TRUE,
#'   image_index = TRUE
#' )
#' }
annh5ad2rds <- function(
    file_path = NULL, python_path = "D:/APP/anaconda3/envs/scanpy/python.exe",
    data_type = "stRNA", convert_mode = "scanpy",
    assay_id = NULL, X_index = "rawX", binsize = 1,
    SCT_index = FALSE, reduction_index = FALSE, image_index = NULL, return_object = TRUE,
    grp_nm = "sample1",dir_nm = "M0_annh5ad2rds"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  clog_start()

  # >>> Check input patameter
  clog_args("Check input parameter")
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm) # 默认不带最后一个斜杠
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir,recursive = T, showWarnings =F)
  # dir.create(photo_dir,recursive = T, showWarnings = F)
  clog_normal(paste0("Your output_dir is ",output_dir))
  # clog_normal(paste0("Your photo_dir is ",photo_dir))

  if (packageVersion("Seurat") < "5.0") {
    clog_error("Please install Seurat version >= 5.0")
  }
  if(data_type == "stRNA"){
    assay_id <- ifelse(is.null(assay_id),"Spatial",assay_id)
    image_index <- ifelse(is.null(image_index),TRUE,image_index)
  }else if(data_type == "scRNA"){
    assay_id <- ifelse(is.null(assay_id),"RNA",assay_id)
    image_index <- FALSE
  }else{
    clog_error("data_type is not correct")
  }
  clog_normal(paste0("Your binsize is ",binsize))
  if(!is.null(python_path)){
    use_python(python_path)
  }else{
    clog_warn("No python_path, use the default python")
    clog_normal(paste0("Your miniconda_path is ",miniconda_path()))
    install_miniconda(path = miniconda_path(),
                      update = TRUE,
                      force = FALSE)
  }
  clog_normal(paste0("Your python is ",py_config()$python))
  # >>> End check

  # >>> Start main pipeline
  clog_step("Step1/5: Read the h5ad file")
  clog_normal("Processing the h5ad")
  file_nm <- gsub(".h5ad", ".rds", basename(file_path))
  clog_normal(paste0("Converting ", basename(file_path)," -> ", file_nm))

  #
  if(convert_mode == "seurat"){
    clog_normal("convert h5ad to rds by seurat mode")
    Convert(file_path, dest = "h5seurat", assay = assay_id, overwrite = TRUE)

    h5file <- paste(paste(unlist(strsplit(file_path, "h5ad", fixed = TRUE)), collapse='h5ad'), "h5seurat", sep="")
    clog_normal(paste(c("Finished! Converting h5ad to h5seurat file at:", h5file), sep=" ", collapse=NULL))

    object <- LoadH5Seurat(h5file,meta.data = F, misc = F) # 不读入misc中的内容，默认只有count数据
    clog_normal(paste(c("Successfully load h5seurat:", h5file), sep=" ", collapse=NULL))

    if(image_index){
      if(!("x" %in% colnames(object@meta.data)) & !("y" %in% colnames(object@meta.data)) &
         !is.null(object@reductions[[assay_id]])){
        clog_normal("Add x and y to meta.data")
        object@meta.data[c("x","y")] <- object@reductions[[assay_id]]@cell.embeddings # seurat格式meta中的x和y一般来源于这里
      }
    }

  }else if(convert_mode == "scanpy"){
    clog_normal("Converting h5ad to rds by scanpy mode")
    ad <- read_h5ad(file_path)
    clog_normal('if an invalid class dgRMatrix object: x slot is not of type double error occurs,
          please use adata.X = adata.raw.X.astype("float32") in python')

    if(is.null(ad$raw) | X_index =="X") {
      clog_normal(paste0("your X_index is ",X_index))
      clog_warn("No raw data, use X as counts")
      X_index <- "X"
      counts <- ad$X
      X_type <- .determine_X_type(counts)
      clog_normal(paste("X matrix is ", X_type))
      if(X_type == "scale.data"){
        clog_normal("Convert scale.data to dgRMatrix")
        counts <- as(counts, "dgRMatrix")
      }
    }else if(X_index == "rawX"){
      counts <- ad$raw$X
      X_type <- .determine_X_type(counts)
      clog_normal(paste("rawX matrix is ", X_type))
      if(X_type == "scale.data"){
        clog_normal("Convert scale.data to dgRMatrix")
        counts <- as(counts, "dgRMatrix")
      }
    }else{
      clog_error("X_index is not correct")
    }

    ad_dim <- dim(counts)
    metadata <- ad$obs
    if(X_index == "X"){
      features <- ad$var
    }else if(X_index == "rawX"){
      features <- ad$raw$var
    }else{
      clog_error("X_index is not correct")
    }

    if(image_index){
      if(!("x" %in% metadata) & !("y" %in% metadata)){
        clog_normal("Add x and y to meta.data")
        coord <- ad$obsm[[str_to_lower(assay_id)]] %>%
          as.data.frame() %>%
          mutate(across(everything(),as.numeric))
        colnames(coord) <- c("x","y")
        metadata <- cbind(metadata,coord)
      }
    }
    # rm(ad);gc()

    counts_dgT <- Matrix::summary(counts)
    rm(counts);gc()
    counts_dgC <- sparseMatrix(i = counts_dgT$i, j = counts_dgT$j,  # 重构dgRMatrix -> dgCMatrix
                               x = counts_dgT$x, dims = ad_dim, repr = "C") %>%

      t()
    rm(counts_dgT);gc()
    dimnames(counts_dgC) <- list(rownames(features), rownames(metadata))
    object <- CreateSeuratObject(counts = counts_dgC,
                                 meta.data = metadata,
                                 assay = assay_id,
                                 names.field = 1,
                                 names.delim = "_")
    rownames(object) <- rownames(features)
    rm(counts_dgC);gc()
  }
  print(paste0("The gene name is : ",paste0(head(rownames(object)),collapse = " ")," ..."))

  #
  clog_step("Step2/5: Add reduction")
  if( reduction_index){
    if(convert_mode == "scanpy"){ # Seurat格式会自动转
      if(!is.null(ad$obsm[["X_umap"]])){
        X_umap <- ad$obsm[["X_umap"]]
        colnames(X_umap) <- paste0("umap_",1:ncol(X_umap))
        rownames(X_umap) <- rownames(ad)
        object@reductions$umap <- CreateDimReducObject(embeddings = X_umap, key = "umap_",assay = assay_id)
      }
      if(!is.null(ad$obsm[["X_tsne"]])){
        X_tsne <- ad$obsm[["X_tsne"]]
        colnames(X_tsne) <- paste0("tsne_",1:ncol(X_tsne))
        rownames(X_tsne) <- rownames(ad)
        object@reductions$tsne <- CreateDimReducObject(embeddings = X_tsne, key = "tsne_",assay = assay_id)
      }
      if(!is.null(ad$obsm[["X_pca"]])){
        X_pca <- ad$obsm[["X_pca"]]
        colnames(X_pca) <- paste0("pca_",1:ncol(X_pca))
        rownames(X_pca) <- rownames(ad)
        object@reductions$pca <- CreateDimReducObject(embeddings = X_pca, key = "pca_",assay = assay_id)
      }
    }
  }else{
    # 不保留任何的reduction
    object@reductions[[assay_id]] <- NULL
  }

  #
  clog_step("Step3/5: Add misc")
  if(SCT_index){
    if (
      !is.null(object@misc$sct_counts) &&
      !is.null(object@misc$sct_data) &&
      !is.null(object@misc$sct_scale) &&
      !is.null(object@misc$sct_cellname) &&
      !is.null(object@misc$sct_genename) &&
      !is.null(object@misc$sct_top_features)
    ) {
      clog_normal("convert stereopy SCT result to seurat SCT result")
      sct.assay.out <- CreateAssayObject(counts=object[['Spatial']]@counts, check.matrix=FALSE)
      sct.assay.out <- SetAssayData(
        object = sct.assay.out,
        slot = "data",
        new.data = log1p(x=GetAssayData(object=sct.assay.out, slot="counts"))
      )
      sct.assay.out@scale.data <- as.matrix(object@misc$sct_scale)
      colnames(sct.assay.out@scale.data) <- object@misc$sct_cellname
      rownames(sct.assay.out@scale.data) <- object@misc$sct_scale_genename
      sct.assay.out <- Seurat:::SCTAssay(sct.assay.out, assay.orig='Spatial')
      Seurat::VariableFeatures(object = sct.assay.out) <- object@misc$sct_top_features
      object[['SCT']] <- sct.assay.out
      DefaultAssay(object=object) <- 'SCT'

      # TODO: tag the reductions as SCT, this will influence the find_cluster choice of data
      object@reductions$pca@assay.used <- 'SCT'
      object@reductions$umap@assay.used <- 'SCT'
      assay.used <- 'SCT'
      clog_normal("Got SCTransform result in object, create a new SCTAssay and set it as default assay.")
    }
  } else {
    assay.used <- assay_id
    clog_normal("Get raw counts only, auto create log-normalize data.")
    object <- NormalizeData(object, assay = assay.used)
    clog_normal("Create log-normalize data.")
  }

  #
  clog_step("Step4/5: Add image")
  if( image_index){
    clog_normal("Start add image to seurat object, This may take some minutes...")
    cell_coords <- unique(object@meta.data[, c('x', 'y')]) # 一般不会重复吧
    cell_coords['cell'] <- row.names(cell_coords)

    # meta中的坐标在seurat中是修改的，但是metadata中的是没有修改的
    clog_normal("Here we do the same three steps for the x and y coordinates:
                  1. Subtract the minimum value 2. Divide by binsize 3. Plus one")
    cell_coords$x <- (cell_coords$x - min(cell_coords$x))/binsize + 1 # 从1开始，rds的坐标都是从1开始的？
    cell_coords$y <- (cell_coords$y - min(cell_coords$y))/binsize + 1

    # object of images$slice1@image, all illustrated as 1 since no concrete pic
    tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))

    # object of images$slice1@coordinates, concrete coordinate of X and Y
    tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                        tissue = 1,
                                        row = cell_coords$y, col = cell_coords$x,
                                        imagerow = cell_coords$y, imagecol = cell_coords$x)
    # @images$slice1@scale.factors
    scalefactors_json <- toJSON(list(fiducial_diameter_fullres = 1,
                                     tissue_hires_scalef = 1,
                                     tissue_lowres_scalef = 1))

    # generate object @images$slice1
    .generate_stereo_spatial <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE) {
      if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
      }
      unnormalized.radius <- scale.factors$fiducial_diameter_fullres *
        scale.factors$tissue_lowres_scalef
      spot.radius <- unnormalized.radius / max(dim(x = image))
      return(new(Class = 'VisiumV1',
                 image = image,
                 scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                              fiducial = scale.factors$fiducial_diameter_fullres,
                                              hires = scale.factors$tissue_hires_scalef,
                                              lowres = scale.factors$tissue_lowres_scalef),
                 coordinates = tissue.positions,
                 spot.radius = spot.radius))
    }

    stereo_spatial <- .generate_stereo_spatial(image = tissue_lowres_image,
                                              scale.factors = fromJSON(scalefactors_json), # 啥作用？
                                              tissue.positions = tissue_positions_list)

    # can be thought of as a background of spatial
    # import image into seurat object
    object@images[['slice1']] <- stereo_spatial
    object@images$slice1@key <- "slice1_"
    object@images$slice1@assay <- assay.used

  }
  clog_normal("The seurat object info:")
  print(object)

  #
  clog_step("Step5/5: Save rds")
  outfile <- paste0(output_dir, "/", file_nm)
  clog_normal(paste0("Save rds to ", outfile))
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  gc()
  saveRDS(object, outfile)
  # >>> End main pipeline

  # >>> End
  clog_end()

  # >>> return
  if(return_object){
    return(object)
  }else{
    return()
  }
}

#' Determine the type of expression matrix from AnnData
#'
#' Internal function to classify the expression matrix type based on matrix class
#' and properties.
#'
#' @param counts Matrix object from AnnData (typically dgRMatrix)
#'
#' @return Character string indicating matrix type: "count", "data", "scale.data",
#'         or "unknown"
#'
#' @keywords internal
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' X_type <- determine_X_type(adata$X)
#' }
.determine_X_type <- function(counts = NULL){
  if(!inherits(counts, "dgRMatrix")){
    return("scale.data")
  }else{
    if(all(counts[1,] %% 1 == 0) & ((colSums(counts) %>% unique() %>% length())>1)){
      return("count")
    }else if((colSums(counts) %>% unique() %>% length()) == 1){
      return("data")
    }else{
      return("unknown")
    }
  }
}


#' Convert Seurat rds file to AnnData h5ad file
#'
#' Converts a Seurat rds file to a Python AnnData h5ad file for use with scanpy
#' and other Python-based tools.
#'
#' @param seurat_obj A Seurat object to convert
#' @param python_path Character, path to Python executable with scanpy installed
#'        (default: "D:/APP/anaconda3/envs/scanpy/python.exe")
#' @param data_type Character, data type - "stRNA" (spatial) or "scRNA" (single-cell)
#'        (default: "stRNA")
#' @param convert_mode Character, conversion method - "seurat" (using SeuratDisk)
#'        (default: "seurat")
#' @param assay_id Character, assay name to convert (default: NULL, auto-detected)
#' @param X_index Character, which matrix to use as counts - "X" or "rawX"
#'        (default: "rawX")
#' @param binsize Numeric, bin size for spatial coordinate scaling (default: 1)
#' @param SCT_index Logical, whether to include SCT-transformed data (default: FALSE)
#' @param reduction_index Logical, whether to include dimensional reductions
#'        (default: FALSE)
#' @param image_index Logical, whether to include spatial image data
#'        (default: NULL, auto-detected)
#' @param return_object Logical, whether to return the object (default: FALSE)
#' @param grp_nm Character, group name for output organization (default: NULL)
#' @param dir_nm Character, directory name for output (default: "M0_rds2annh5ad")
#'
#' @return NULL (invisible), saves h5ad file to disk
#'
#' @import reticulate
#' @import Matrix
#' @import dplyr
#' @import rjson
#' @import Seurat
#' @import ggplot2
#' @importFrom SeuratDisk SaveH5Seurat Convert LoadH5Seurat
#' @import SeuratObject
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Convert Seurat object to h5ad
#' rds2annh5ad(
#'   seurat_obj = seurat_object,
#'   python_path = "~/anaconda3/envs/scanpy/bin/python",
#'   data_type = "stRNA",
#'   grp_nm = "sample1"
#' )
#' }
rds2annh5ad <- function(
    seurat_obj = NULL, python_path = "D:/APP/anaconda3/envs/scanpy/python.exe",
    data_type = "stRNA", convert_mode = "seurat",
    assay_id = NULL, X_index = "rawX", binsize = 1,
    SCT_index = FALSE, reduction_index = FALSE, image_index = NULL, return_object = FALSE,
    grp_nm = NULL,dir_nm = "M0_rds2annh5ad"
){
  on.exit(while(sink.number() > 0){sink()}, add = TRUE)

  # >>> Start pipeline
  tmp_file <- tempfile()
  sink(tmp_file,split = TRUE)
  clog_start()

  # >>> Check input patameter
  clog_args("Check input arguments")
  if (packageVersion("Seurat") < "5.0") {
    clog_error("Please install Seurat version >= 5.0")
  }
  # .check_null_args(geneset_list)
  now_time <- format(Sys.time(), "%Y%m%d_%H%M%S")
  if(is.null(grp_nm)){
    grp_nm <- now_time
  }

  if(data_type == "stRNA"){
    assay_id <- "Spatial"
    image_index <- TRUE
  }else if(data_type == "scRNA"){
    assay_id <- "RNA"
    image_index <- FALSE
  }else{
    clog_error("data_type is not correct")
  }
  clog_normal(paste0("Your binsize is ",binsize))
  if(!is.null(python_path)){
    use_python(python_path)
  }else{
    clog_warn("No python_path, use the default python")
    clog_normal(paste0("Your miniconda_path is ",miniconda_path()))
    install_miniconda(path = miniconda_path(),
                      update = TRUE,
                      force = FALSE)
  }
  clog_normal(paste0("Your python is ",py_config()$python))
  # >>> End check

  # >>> dir
  clog_normal("Create outputdata and photo directories")
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm,"/")
  # photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm,"/")
  logic <- .create_directory_interactive(output_dir)
  if(!logic) {return(invisible(NULL))}
  # logic <- .create_directory_interactive(photo_dir)
  # if(!logic) {return(invisible(NULL))}
  clog_normal(paste0("Your output_dir is ",output_dir))
  # clog_normal(paste0("Your photo_dir is ",photo_dir))

  # >>> Start main pipeline
  clog_step("Step1/5: Read the seurat object")
  clog_normal("Processing the seurat object")
  file_nm <- as.character(substitute(seurat_obj))
  clog_normal(paste0("Converting ", file_nm, ".rds -> ", file_nm, ".h5ad"))

  #
  if(convert_mode == "seurat"){
    clog_normal("convert rds to h5ad by seurat mode")
    seurat_obj@meta.data <- seurat_obj@meta.data %>%
      mutate(across(where(is.factor), as.character)) # 转为character，否则h5ad中只保留数字格式
    seurat_obj[[assay_id]] <- as(seurat_obj[[assay_id]], Class = "Assay")
    seurat_obj@assays[[assay_id]]["scale.data"] <- NULL
    clog_warn("Remove the scale.data")

    clog_step("Step2/5: SaveH5Seurat")
    SaveH5Seurat(seurat_obj,filename="test.h5seurat", overwrite = TRUE)
    clog_step("Step3/5: Convert to h5ad")
    outfile <- paste0(output_dir, "/", file_nm,".h5ad")
    clog_normal(paste0("Save h5ad to ", outfile))
    Convert("test.h5seurat",
            dest = outfile,
            overwrite = TRUE)
    file.remove("test.h5seurat")

  }else if(convert_mode == "scanpy"){
    clog_normal("Converting rds to h5ad by scanpy mode")
    # 剩余代码未开发
  }
  # >>> End main pipeline

  # >>> Final
  .save_function_params("rds2annh5ad", envir = environment(), file = paste0(output_dir,"Log_function_params_(rds2annh5ad).log") )
  clog_end()
  sink()
  file.rename(tmp_file, paste0(output_dir,"Log_termial_output_(rds2annh5ad).log")) %>% invisible()
  invisible(NULL)
}





#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# co_localization
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

co_localization <- function(
    object = NULL, features = NULL, celltype_col = NULL,
    conda_path = NULL,
    Microbe = "Microbe",org = "mouse",
    grp_nm = NULL,dir_nm = "M2_Co_Localization"
){
  # suppressMessages({
  #   library(reticulate)
  #   # library(anndata) # sc读取的和R中的anndata格式不同
  # })
  if(is.null(grp_nm)){
    grp_nm <- Microbe
  }
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm,"/")
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm,"/")
  dir.create(output_dir,recursive = T, showWarnings =F)
  dir.create(photo_dir,recursive = T, showWarnings = F)
  dir.create("./db",recursive = T, showWarnings = F)

  # >>> Check input
  if(is.null(conda_path)){
    install_miniconda(path = miniconda_path(),
                      update = TRUE,
                      force = FALSE)
    use_condaenv(miniconda_path())
    conda_install(
      # envname = "my_r_python", # 默认读取的就是当前的环境，use_condaenv决定，也可通过py_config()确定
      packages = c("scanpy",
                   "squidpy",
                   "matplotlib",
                   "seaborn",
                   "pandas",
                   "numpy"
      )
    )
  }else{
    # conda_path <- "D:/APP/Anaconda3/envs/scanpy"
    use_condaenv(conda_path)
  }
  if(py_available()){
    cat("The following is your py_config: \n")
    print(py_config())
  }
  # >>> End check

  # >>> Start pipeline
  cat(paste0("\n>>> Co_localization start at ", Sys.time()),"\n\n")
  # Reference:
  # https://kayla-morrell.github.io/SquidpyR/articles/C_reproduceSquidpyFigures.html

  #
  cat(paste0(">>> Step1/4: Import python package\n"))
  cat(paste0("Import python package start at ", Sys.time()),"\n")
  sc <- reticulate::import("scanpy")
  sq <- reticulate::import("squidpy")
  plt <- reticulate::import("matplotlib.pyplot")
  sns <- reticulate::import("seaborn") # squidpy使用seaborn绘图
  pd <- reticulate::import("pandas")
  np <- reticulate::import("numpy")
  sns$set_theme(style = "white")  # 设置白色背景主题
  sc$logging$print_header()
  sc$set_figure_params(facecolor="white", figsize = c(8,8))
  sc$settings$verbosity <- 'hint'
  sc$settings$figdir <- photo_dir
  PATH_FIGURES <- photo_dir

  #
  cat(paste0(">>> Step2/4: Read h5ad data\n"))
  cat(paste0("Read h5ad data start at ", Sys.time()),"\n")
  adata <- sq$datasets$seqfish() # demo
  adata
  class(adata)
  celltype_col <- "celltype_mapped_refined"
  #

  pdf(file = paste0(photo_dir,"/sc_celltype_umap.pdf"),
      width = 6,height = 6)
  sc$pl$spatial(
    adata,
    # color = celltype_col,
    # save = "_seqfish.png",
    spot_size = 0.03,
    show = TRUE,save="pdf"
  )
  dev.off()

  #
  cat(paste0(">>> Step2/4: nhood_enrichment\n"))
  cat(paste0("nhood_enrichment start at ", Sys.time()),"\n")
  sq$gr$spatial_neighbors(adata, delaunay = TRUE) # added to adata.obsp, adata.uns
  sq$gr$nhood_enrichment(adata, cluster_key = "celltype_mapped_refined",
                         n_perms = 100) # added to adata.uns
  arr <- np$array(adata$uns["celltype_mapped_refined_nhood_enrichment"]["zscore"])
  sq$pl$nhood_enrichment(
    adata,
    cluster_key = "celltype_mapped_refined",
    cmap = "coolwarm",
    title = "",
    method = "ward",
    dpi = 300,
    figsize = c(5, 4),
    save = "nhood_seqfish.png",
    cbar_kwargs = list(label = "Z-score"),
    vmin = -50,
    vmax = 50
  )

  #
  sq$pl$nhood_enrichment(
    adata,
    cluster_key = "celltype_mapped_refined",
    cmap = "inferno",
    title = "",
    method = "ward",
    dpi = 300,
    figsize = c(5, 4),
    save = "nhod_seqfish.png",
    cbar_kwargs = list(label = "Z-score"),
    vmin = -50,
    vmax = 150
  )

  #
  ax <- plt$subplots(
    figsize = c(3, 6)
  )
  sc$pl$spatial(
    adata,
    color = "celltype_mapped_refined",
    groups = c(
      "Endothelium",
      "Haematoendothelial progenitors",
      "Allantois",
      "Lateral plate mesoderm",
      "Intermediate mesoderm",
      "Presomitic mesoderm",
      "Dermomyotome"
    ),
    save = "_endo_seqfish.png",
    frameon = FALSE,
    #    ax = ax,
    show = FALSE,
    spot_size = 0.03
  )

  # 共线概率图
  sq$gr$co_occurrence(
    adata,
    cluster_key = "celltype_mapped_refined"
  ) # added to adata.uns
  sq$pl$co_occurrence(
    adata,
    cluster_key = "celltype_mapped_refined",
    clusters = "Lateral plate mesoderm",
    figsize = c(5, 3),
    legend = FALSE,
    dpi = 300,
    save = "co_occurrence_seqfish"
  )


  pdf(file = paste0(photo_dir,"PPI_network_all_show_gene.pdf"),
      width = 6,height = 6)
  print(p2)
  dev.off()

  # >>> End pipeline
  cat(paste0("\n>>> Co_localization end at ", Sys.time()),"\n\n")

  # return()
}






