
## web upload

### 10x
# mat1 = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/V1_Mouse_Brain_Sagittal_Posterior_Section_2_filtered_feature_bc_matrix.h5'
# image1 = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/spatial/tissue_lowres_image.png'
# json1 = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/spatial/scalefactors_json.json'
# coordinates1 = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/spatial/tissue_positions.csv'
# mat_r1 = NULL
# mat_c1 = NULL
# #
# mat2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/filtered_feature_bc_matrix/matrix.mtx.gz'
# mat_r2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/filtered_feature_bc_matrix/features.tsv.gz'
# mat_c2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'
# image2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/spatial/tissue_lowres_image.png'
# json2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/spatial/scalefactors_json.json'
# coordinates2 = '/mnt/raid62/pzz/shiny_data/seurat/tosce/outs2/anterior/outs/spatial/tissue_positions_list.csv'
# obj2 <-
#   load_spatial_10x_web(
#     mat1 = mat1,
#     mat_r1 = mat_r1,
#     mat_c1 = mat_c1,
#     image1 = image1,
#     json1 = json1,
#     coordinates1 = coordinates1,
#     mat2 = mat2,
#     mat_r2 = mat_r2,
#     mat_c2 = mat_c2,
#     image2 = image2,
#     json2 = json2,
#     coordinates2 = coordinates2,
#     species = 'mm',
#     resolution = 0.4
#   )

load_spatial_10x <- function(mat1 = NULL,
                             mat_r1 = NULL,
                             mat_c1 = NULL,
                             image1 = NULL,
                             json1 = NULL,
                             coordinates1 = NULL,
                             min = 3) {
  if (!require('shinySRT')) {
    install.packages('shinySRT')
  }
  if (length(grep(pattern = '.png', image1)) == 0 &
      length(grep(pattern = '.json', json1)) == 0 &
      length(grep(pattern = 'position', coordinates1)) == 0) {
    validate('Missing key documents!!')
  }
  
  if (length(grep(pattern = '.mtx', mat1)) > 0 &
      length(mat_r1) > 0 & length(mat_c1) > 0) {
    data <- Matrix::readMM(mat1)
    rownames(data) <-
      read.csv(mat_r1, header = F, sep = '\t')$V2 %>% make.unique()
    colnames(data) <- read.csv(mat_c1, header = F, sep = '\t')$V1
  } else if (length(grep(pattern = '.h5', mat1)) > 0) {
    data <- Read10X_h5(mat1)
  } else{
    validate('Missing matrix documents!!')
  }
  
  obj <- CreateSeuratObject(data, min.cells = min)
  img <- png::readPNG(image1)
  scale_factors <- jsonlite::fromJSON(txt = json1)
  tissue.positions <- read.csv(
    file = coordinates1,
    col.names = c("barcodes", "tissue", "row", "col", "imagerow",
                  "imagecol"),
    header = ifelse(
      test = basename(coordinates1) ==
        basename(coordinates1),
      yes = TRUE,
      no = FALSE
    ),
    as.is = TRUE,
    row.names = 1
  )
  
  tissue.positions <-
    tissue.positions[which(x = tissue.positions$tissue == 1), ]
  unnormalized.radius <-
    scale_factors$fiducial_diameter_fullres *
    scale_factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(x = img))
  
  
  img_obj <- new(
    Class = "VisiumV1",
    image = img,
    scale.factors = Seurat::scalefactors(
      spot = scale_factors$spot_diameter_fullres,
      fiducial = scale_factors$fiducial_diameter_fullres,
      hires = scale_factors$tissue_hires_scalef,
      scale_factors$tissue_lowres_scalef
    ),
    coordinates = tissue.positions,
    spot.radius = spot.radius
  )
  
  img_obj <- img_obj[Seurat::Cells(x = obj)]
  Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
  obj[['slice']] <- img_obj
  return(obj)
}

load_spatial_10x_web <-
  function(mat1 = NULL,
           mat_r1 = NULL,
           mat_c1 = NULL,
           image1 = NULL,
           json1 = NULL,
           coordinates1 = NULL,
           mat2 = NULL,
           mat_r2 = NULL,
           mat_c2 = NULL,
           image2 = NULL,
           json2 = NULL,
           coordinates2 = NULL,
           species = 'hg',
           npcs = 20,
           resolution = 1) {
    if (!require('shinySRT')) {
      install.packages('shinySRT')
    }
    if (length(image1) == 0 &
        length(json1) == 0 & length(coordinates1) == 0) {
      obj1_class <- F
    } else{
      if (length(grep(pattern = 'matrix.mtx', mat1)) > 0 &
          length(mat_r1) > 0 & length(mat_c1) > 0) {
        obj1_class <- T
      } else if (length(grep(pattern = 'h5$', mat1)) > 0) {
        mat_r1 = NULL
        mat_c1 = NULL
        obj1_class <- T
      } else{
        obj1_class <- F
      }
    }
    
    if (length(image2) == 0 &
        length(json2) == 0 & length(coordinates2) == 0) {
      obj2_class <- F
    } else{
      if (length(grep(pattern = 'matrix.mtx', mat2)) > 0 &
          length(mat_r2) > 0 & length(mat_c2) > 0) {
        obj2_class <- T
      } else if (length(grep(pattern = 'h5$', mat2)) > 0) {
        mat_r2 = NULL
        mat_c2 = NULL
        obj2_class <- T
      } else{
        obj2_class <- F
      }
    }
    
    if (obj1_class) {
      object1 <-
        load_spatial_10x(mat1, mat_r1, mat_c1, image1, json1, coordinates1)
    } else{
      object1 <- NULL
    }
    
    if (obj2_class) {
      object2 <-
        load_spatial_10x(mat2, mat_r2, mat_c2, image2, json2, coordinates2)
    } else{
      object2 <- NULL
    }
    
    
    if (is.null(object1) & is.null(object2)) {
      shiny::validate('Input error')
    }
    
    if(!is.null(object1) & is.null(object2)){
      obj <- object1
    }else if(!is.null(object2) & is.null(object2)){
      obj <- object2
    }else{
      obj <- Reduce(merge, list(object1, object2))
    }
    
    # obj <- Reduce(merge, list(object1, object2))
    # return(list(object1, object2))
    
    SRT_object <- obj_list_process(
      obj_list = list(obj),
      species = species,
      resolution = resolution,
      npcs = npcs
    )
    return(SRT_object)
  }




load_spatial_obj_web <-
  function(obj1 = NULL,
           obj2 = NULL,
           species = 'hg',
           resolution = 1,
           npcs = 20) {
    if (!require('shinySRT')) {
      install.packages('shinySRT')
    }
    obj1_class <- NULL
    if (class(obj1) == 'Seurat') {
      obj1_class <- 'Seurat'
    } else if (class(obj1) == 'list') {
      if (names(obj1)[1] == 'seurat_obj') {
        obj1_class <- 'list'
      }
    }
    
    obj2_class <- NULL
    if (class(obj2) == 'Seurat') {
      obj2_class <- 'Seurat'
    } else if (class(obj2) == 'list') {
      if (names(obj2)[1] == 'seurat_obj') {
        obj2_class <- 'list'
      }
    }
    
    if (length(obj1_class) == 0 & length(obj1_class) == 0) {
      shiny::validate("Input error!")
    }
    
    if (length(obj1_class) == 0) {
      obj1 <- NULL
    }
    
    if (length(obj2_class) == 0) {
      obj2 <- NULL
    }
    
    if(!is.null(obj1) & is.null(obj2)){
      obj <- obj1
    }else if(!is.null(obj2) & is.null(obj2)){
      obj <- obj2
    }else{
      obj <- Reduce(merge, list(obj1, obj2))
    }
    
    
    SRT_object <- obj_list_process(
      obj_list = list(obj),
      species = species,
      resolution = resolution,
      npcs = npcs
    )
    return(SRT_object)
    # return(list(obj1, obj2))
  }

# matx1 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/matrix.txt'
# meta1 = NULL
# coordi1 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/position.csv'
# image1 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/image.png'
# x_reverse1 = T
# 
# 
# matx2 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/matrix.txt'
# meta2 = NULL
# coordi2 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/position.csv'
# image2 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/image.png'
# x_reverse2 = T



# matx = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/matrix.txt'
# meta = NULL
# coordi = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/position.csv'
# image = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/image.png'
# x_reverse = T


# matx2 = NULL
# meta2 = NULL
# coordi2 = NULL
# image2 = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/image.png'
# x_reverse2 = T

load_spatial_mat_web <-
  function(matx1 = NULL,
           meta1 = NULL,
           coordi1 = NULL,
           image1 = NULL,
           x_reverse1 = F,
           matx2 = NULL,
           meta2 = NULL,
           coordi2 = NULL,
           image2 = NULL,resolution = 1,
           x_reverse2 = F,
           min_cells = 3,npcs = 20,
           species = 'hg') {
    if (!require('shinySRT')) {
      install.packages('shinySRT')
    }
    if (length(coordi1) == 0 & length(matx1) == 0) {
      obj1_class <- F
    } else{
      if (length(grep(pattern = 'txt|csv|tsv|xlsx', matx1)) > 0) {
        obj1_class <- T
      } else{
        obj1_class <- F
      }
    }
    
    if (length(coordi2) == 0 & length(matx2) == 0) {
      obj2_class <- F
    } else{
      if (length(grep(pattern = 'txt|csv|tsv|xlsx', matx2)) > 0) {
        obj2_class <- T
      } else{
        obj2_class <- F
      }
    }
    
    if (obj1_class) {
      object1 <-
        load_spatial(
          matx = matx1,
          meta = meta1,
          coordi = coordi1,
          image = image1,
          x_reverse = x_reverse1,
          species = species,
          min_cells = min_cells
        )
    } else{
      object1 <- NULL
    }
    
    if (obj2_class) {
      object2 <-
        load_spatial(
          matx = matx2,
          meta = meta2,
          coordi = coordi2,
          image = image2,
          x_reverse = x_reverse2,
          species = species,
          min_cells = min_cells
        )
    } else{
      object2 <- NULL
    }
    
    if (is.null(object1) & is.null(object2)) {
      shiny::validate('Input error')
    }
    
    
    if(!is.null(object1) & is.null(object2)){
      obj <- list(object1)
    }else if(!is.null(object2) & is.null(object2)){
      obj <- list(object2)
    }else if(class(object1) != class(object2)){
      obj <- list(object1, object2)
    } else{
      obj <- Reduce(merge, list(object1, object2))
      obj <- list(obj)
    }
    
    
    SRT_object <- obj_list_process(
      obj_list = obj,
      species = species,
      resolution = resolution,
      npcs = npcs
    )
    return(SRT_object)
  }


#############################################################
#############################################################
#############################################################
# ## load SRT data(dir)
# spatial_load_dir <- function(x,
#                              species = 'hg',
#                              x_reverse = F,
#                              min_cells = 3,
#                              ...) {
#   if (length(grep(pattern = 'filtered_feature_bc_matrix.h5$', dir(x))) == 1) {
#     ## Raw data processing requires seurat
#     h5 <-
#       grep(pattern = 'filtered_feature_bc_matrix.h5$', dir(x), value = T)
#     data <- Seurat::Read10X_h5(filename = file.path(x, h5))
#     object <-
#       Seurat::CreateSeuratObject(counts = data,
#                                  assay = 'Spatial',
#                                  min.cells = min_cells)
#     img_obj <-
#       Seurat::Read10X_Image(image.dir = file.path(x, "spatial"),
#                             filter.matrix = T)
#     img_obj <- img_obj[Seurat::Cells(x = object)]
#     Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
#     object[[x]] <- img_obj
#     
#   } else if (length(which(dir(x) == "filtered_feature_bc_matrix")) == 1) {
#     data <- Seurat::Read10X(file.path(x, "filtered_feature_bc_matrix"))
#     object <-
#       Seurat::CreateSeuratObject(data, assay = 'Spatial', min.cells = min_cells)
#     
#     img_obj <-
#       Seurat::Read10X_Image(image.dir = file.path(x, "spatial"),
#                             filter.matrix = T)
#     img_obj <- img_obj[Seurat::Cells(x = object)]
#     
#     Seurat::DefaultAssay(object = img_obj) <- 'Spatial'
#     object[[x]] <- img_obj
#     
#   } else if (length(grep(pattern = 'matrix', dir(x))) > 0) {
#     mtx_file <- Sys.glob(paths = file.path(x, "*matrix*"))
#     meta <- Sys.glob(paths = file.path(x, "*meta*"))
#     scale_factors_file <-
#       Sys.glob(paths = file.path(x, "scalefactors_json.json"))
#     tissue.positions.path <-
#       Sys.glob(paths = file.path(x, "*position*"))
#     img_file <- Sys.glob(paths = file.path(x, "*image*"))
#     object <- load_spatial(
#       matx = mtx_file,
#       meta = meta,
#       coordi = tissue.positions.path,
#       image = img_file,
#       x_reverse = x_reverse,
#       species = species,
#       scale_factors_file = scale_factors_file
#     )
#   }
#   return(object)
# }

# obj3 <- multi_dir_spatial(multi_dir = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/',species = 'mm')

#####################################





# obj <-
#   single_op_dir(dir = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/', species = 'mm')
# obj_ <-
#   single_op_dir(dir = '/mnt/raid62/pzz/shiny_data/seurat/multi_raw/1/', species = 'mm',resolution = 2)




# obj2 <-
#   single_op_file(
#     matx = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/matrix.txt',
#     coordi = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/position.csv',
#     image = '/mnt/raid62/pzz/shiny_data/SPATA2data/spata/image.png',
#     species = 'mm',
#     x_reverse = T
#   )


generate_random_string <- function(length){
  charset <- c(letters,LETTERS, 0, 9)
  
  random_indices <- sample(1:length(charset),length, replace = T)
  random_str <- paste(charset[random_indices], collapse = '')
  return(random_str)
}
