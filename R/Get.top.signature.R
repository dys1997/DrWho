#' Sort matrix based on corrected p-value and log2 fold change
#'
#'
#' @param mtx Matrix of differentially expressed genes
#'
#' @return a sorted matrix
#' @export
#'
Sort_mtx <- function(mtx) {
  sorted_indices <- order(mtx[, "p_val_adj"],-mtx[, "avg_log2FC"] )
  sorted_matrix <- mtx[sorted_indices, ]
  return(sorted_matrix)
}
#' Obtain the specified number of genes as cancer signatures used for enrichment analysis
#'
#'
#' @return a matrix with cancer signature information
#' @export
Get.top.signature <- function(top.num=100) {
  ###读取数据###


  DEG.Glioblastoma  <- system.file("data", "DEG.Glioblastoma.rds", package = "ISEA")
  DEG.Glioblastoma = readRDS(DEG.Glioblastoma)
  DEG.H3K27M_gliomas  <- system.file("data", "DEG.H3K27M_gliomas.rds", package = "ISEA")
  DEG.H3K27M_gliomas = readRDS(DEG.H3K27M_gliomas)
  DEG.BREAST  <- system.file("data", "DEG.BREAST.rds", package = "ISEA")
  DEG.BREAST = readRDS(DEG.BREAST)
  DEG.PDAC  <- system.file("data", "DEG.PDAC.rds", package = "ISEA")
  DEG.PDAC=readRDS(DEG.PDAC)
  DEG.head_neck <- system.file("data", "DEG.head_neck.rds", package = "ISEA")
  DEG.head_neck=readRDS(DEG.head_neck)
  DEG.Prostate   <- system.file("data", "DEG.Prostate.rds", package = "ISEA")
  DEG.Prostate=readRDS(DEG.Prostate)
  DEG.lung   <- system.file("data", "DEG.lung.rds", package = "ISEA")
  DEG.lung=readRDS(DEG.lung)
  DEG.liver   <- system.file("data", "DEG.liver.rds", package = "ISEA")
  DEG.liver=readRDS(DEG.liver)
  DEG.Kindey   <- system.file("data", "DEG.Kindey.rds", package = "ISEA")
  DEG.Kindey= readRDS(DEG.Kindey)
  DEG.HNSCC   <- system.file("data", "DEG.HNSCC.rds", package = "ISEA")
  DEG.HNSCC=readRDS(DEG.HNSCC)
  DEG.Sarcoma_Osteosarcoma   <- system.file("data", "DEG.Sarcoma_Osteosarcoma.rds", package = "ISEA")
  DEG.Sarcoma_Osteosarcoma=readRDS(DEG.Sarcoma_Osteosarcoma)
  DEG.ALL <- system.file("data", "DEG.ALL.rds", package = "ISEA")
  DEG.ALL=readRDS(DEG.ALL)
  DEG.Multiple_myeloma   <- system.file("data", "DEG.Multiple_myeloma.rds", package = "ISEA")
  DEG.Multiple_myeloma=readRDS(DEG.Multiple_myeloma)
  DEG.CML  <- system.file("data", "DEG.CML.rds", package = "ISEA")
  DEG.CML=readRDS(DEG.CML)
  DEG.Neuroblastoma  <- system.file("data", "DEG.Neuroblastoma.rds", package = "ISEA")
  DEG.Neuroblastoma=readRDS(DEG.Neuroblastoma)
  DEg.CRC  <- system.file("data", "DEg.CRC.rds", package = "ISEA")
  DEg.CRC =readRDS(DEg.CRC)
  DEG.Myeloproliferative_neoplasms  <- system.file("data", "DEG.Myeloproliferative_neoplasms.rds", package = "ISEA")
  DEG.Myeloproliferative_neoplasms = readRDS(DEG.Myeloproliferative_neoplasms)
  DEG.Ovarian_cancer  <- system.file("data", "DEG.Ovarian_cancer.rds", package = "ISEA")
  DEG.Ovarian_cancer = readRDS(DEG.Ovarian_cancer)


  ###筛选差异表达基因###

  DEG.ALL.used = DEG.ALL[DEG.ALL$avg_log2FC>0,]
  DEG.ALL.used = Sort_mtx(DEG.ALL.used)
  if(top.num > length(rownames(DEG.ALL.used))){
    DEG.ALL.used = DEG.ALL.used
  }else{
    DEG.ALL.used = DEG.ALL.used[1:top.num,]
  }

  DEG.Multiple_myeloma.used =DEG.Multiple_myeloma[DEG.Multiple_myeloma$avg_log2FC>0,]
  DEG.Multiple_myeloma.used = Sort_mtx(DEG.Multiple_myeloma.used)
  if(top.num > length(rownames(DEG.Multiple_myeloma.used))){
    DEG.Multiple_myeloma.used = DEG.Multiple_myeloma.used
  }else{
    DEG.Multiple_myeloma.used = DEG.Multiple_myeloma.used[1:top.num,]
  }

  DEG.CML.used = DEG.CML[DEG.CML$avg_log2FC>0,]
  DEG.CML.used = Sort_mtx(DEG.CML.used)
  if(top.num > length(rownames(DEG.CML.used))){
    DEG.CML.used = DEG.CML.used
  }else{
    DEG.CML.used = DEG.CML.used[1:top.num,]
  }

  DEG.Neuroblastoma.used = DEG.Neuroblastoma[DEG.Neuroblastoma$avg_log2FC>0,]
  DEG.Neuroblastoma.used = Sort_mtx(DEG.Neuroblastoma.used)
  if(top.num > length(rownames(DEG.Neuroblastoma.used))){
    DEG.Neuroblastoma.used = DEG.Neuroblastoma.used
  }else{
    DEG.Neuroblastoma.used = DEG.Neuroblastoma.used[1:top.num,]
  }

  DEg.CRC.used = DEg.CRC[DEg.CRC$avg_log2FC>0,]
  DEg.CRC.used =Sort_mtx(DEg.CRC.used)
  if(top.num > length(rownames(DEg.CRC.used))){
    DEg.CRC.used = DEg.CRC.used
  }else{
    DEg.CRC.used = DEg.CRC.used[1:top.num,]
  }

  DEG.Myeloproliferative_neoplasms.used = DEG.Myeloproliferative_neoplasms[DEG.Myeloproliferative_neoplasms$avg_log2FC>0,]
  DEG.Myeloproliferative_neoplasms.used = Sort_mtx(DEG.Myeloproliferative_neoplasms.used)
  if(top.num > length(rownames(DEG.Myeloproliferative_neoplasms.used))){
    DEG.Myeloproliferative_neoplasms.used = DEG.Myeloproliferative_neoplasms.used
  }else{
    DEG.Myeloproliferative_neoplasms.used = DEG.Myeloproliferative_neoplasms.used[1:top.num,]
  }

  DEG.Ovarian_cancer.used  = DEG.Ovarian_cancer[DEG.Ovarian_cancer$avg_log2FC>0,]
  DEG.Ovarian_cancer.used = Sort_mtx(DEG.Ovarian_cancer.used)
  if(top.num > length(rownames(DEG.Ovarian_cancer.used))){
    DEG.Ovarian_cancer.used = DEG.Ovarian_cancer.used
  }else{
    DEG.Ovarian_cancer.used = DEG.Ovarian_cancer.used[1:top.num,]
  }
  DEG.Glioblastoma.used = DEG.Glioblastoma[DEG.Glioblastoma$avg_log2FC >0,]
  DEG.Glioblastoma.used = Sort_mtx(DEG.Glioblastoma.used)
  if(top.num > length(rownames(DEG.Glioblastoma.used))){
    DEG.Glioblastoma.used = DEG.Glioblastoma.used
  }else{
    DEG.Glioblastoma.used = DEG.Glioblastoma.used[1:top.num,]
  }
  DEG.H3K27M_gliomas.used = DEG.H3K27M_gliomas[DEG.H3K27M_gliomas$avg_log2FC >0,]
  DEG.H3K27M_gliomas.used = Sort_mtx(DEG.H3K27M_gliomas.used)
  if(top.num > length(rownames(DEG.H3K27M_gliomas.used))){
    DEG.H3K27M_gliomas.used = DEG.H3K27M_gliomas.used
  }else{
    DEG.H3K27M_gliomas.used = DEG.H3K27M_gliomas.used[1:top.num,]
  }
  DEG.head_neck.used = DEG.head_neck[DEG.head_neck$avg_log2FC >0,]
  DEG.head_neck.used =Sort_mtx(DEG.head_neck.used)
  if(top.num > length(rownames(DEG.head_neck.used))){
    DEG.head_neck.used = DEG.head_neck.used
  }else{
    DEG.head_neck.used = DEG.head_neck.used[1:top.num,]
  }
  DEG.Prostate.used = DEG.Prostate[DEG.Prostate$avg_log2FC>0,]
  DEG.Prostate.used = Sort_mtx(DEG.Prostate.used)
  if(top.num > length(rownames(DEG.Prostate.used))){
    DEG.Prostate.used = DEG.Prostate.used
  }else{
    DEG.Prostate.used = DEG.Prostate.used[1:top.num,]
  }
  DEG.lung.used = DEG.lung[DEG.lung$avg_log2FC>0,]
  DEG.lung.used = Sort_mtx(DEG.lung.used)
  if(top.num > length(rownames(DEG.lung.used))){
    DEG.lung.used = DEG.lung.used
  }else{
    DEG.lung.used = DEG.lung.used[1:top.num,]
  }
  DEG.liver.used = DEG.liver[DEG.liver$avg_log2FC>0,]
  DEG.liver.used = Sort_mtx(DEG.liver.used)
  if(top.num > length(rownames(DEG.liver.used))){
    DEG.liver.used = DEG.liver.used
  }else{
    DEG.liver.used = DEG.liver.used[1:top.num,]
  }
  DEG.Kindey.used = DEG.Kindey[DEG.Kindey$avg_log2FC >0,]
  DEG.Kindey.used = Sort_mtx(DEG.Kindey.used)
  if(top.num > length(rownames(DEG.Kindey.used))){
    DEG.Kindey.used = DEG.Kindey.used
  }else{
    DEG.Kindey.used = DEG.Kindey.used[1:top.num,]
  }
  DEG.HNSCC.used = DEG.HNSCC[DEG.HNSCC$avg_log2FC>0,]
  DEG.HNSCC.used = Sort_mtx(DEG.HNSCC.used)
  if(top.num > length(rownames(DEG.HNSCC.used))){
    DEG.HNSCC.used = DEG.HNSCC.used
  }else{
    DEG.HNSCC.used = DEG.HNSCC.used[1:top.num,]
  }
  DEG.PDAC.used = DEG.PDAC[DEG.PDAC$avg_log2FC>0,]
  DEG.PDAC.used = Sort_mtx(DEG.PDAC.used)
  if(top.num > length(rownames(DEG.PDAC.used))){
    DEG.PDAC.used = DEG.PDAC.used
  }else{
    DEG.PDAC.used = DEG.PDAC.used[1:top.num,]
  }
  DEG.Sarcoma_Osteosarcoma.used = DEG.Sarcoma_Osteosarcoma[DEG.Sarcoma_Osteosarcoma$avg_log2FC>0,]
  DEG.Sarcoma_Osteosarcoma.used = Sort_mtx(DEG.Sarcoma_Osteosarcoma.used)
  if(top.num > length(rownames(DEG.PDAC.used))){
    DEG.Sarcoma_Osteosarcoma.used = DEG.Sarcoma_Osteosarcoma.used
  }else{
    DEG.Sarcoma_Osteosarcoma.used = DEG.Sarcoma_Osteosarcoma.used[1:top.num,]
  }
  DEG.BREAST.used = DEG.BREAST[DEG.BREAST$avg_log2FC > 0,]
  DEG.BREAST.used = Sort_mtx(DEG.BREAST.used)
  if(top.num > length(rownames(DEG.BREAST.used))){
    DEG.BREAST.used = DEG.BREAST.used
  }else{
    DEG.BREAST.used = DEG.BREAST.used[1:top.num,]
  }


  ###创建gene matrix####

  cancer.list = c("Glioblastoma(GBM)",
                  "Gliomas","Breast Cancer(BRCA)",
                  "Pancreatic ductal adenocarcinoma(PDAC)","Nasopharyngeal carcinoma(NPC)",
                  "Prostate cancer(PRAD)",
                  "Lung adenocarcinoma(LUAD)",
                  "liver cancaer(HCC)",
                  "Kidney cancaer(RCC)",
                  "Head and neck squamous cell carcinoma(HNSCC)",
                  "Osteosarcoma",
                  "Acute lymphoblastic leukemia(ALL)",
                  "Multiple myeloma(MM)",
                  "Chronic myelogenous leukemia(CML)",
                  "Neuroblastoma(NB)",
                  "Colorectal cancer(CRC)",
                  "Myeloproliferative neoplasms(MPN)",
                  "Ovarian cancer(OV)")
  empty_matrix <- matrix(nrow = length(cancer.list), ncol = 2)
  colnames(empty_matrix) <- c("cancer", "gene_list")
  empty_matrix = as.data.frame(empty_matrix)
  rownames(empty_matrix) = cancer.list
  empty_matrix$cancer = cancer.list
  empty_matrix[1,2] = paste(DEG.Glioblastoma.used$X,collapse = ",")
  empty_matrix[2,2]= paste(DEG.H3K27M_gliomas.used$X,collapse = ",")
  empty_matrix[3,2]= paste(DEG.BREAST.used$X,collapse = ",")
  empty_matrix[4,2]=paste(DEG.PDAC.used$X,collapse = ",")
  empty_matrix[5,2]=paste(DEG.head_neck.used$X,collapse = ",")
  empty_matrix[6,2]=paste(DEG.Prostate.used$X,collapse = ",")
  empty_matrix[7,2]=paste(DEG.lung.used$X,collapse = ",")
  empty_matrix[8,2]=paste(DEG.liver.used$X,collapse = ",")
  empty_matrix[9,2]=paste(DEG.Kindey.used$X,collapse = ",")
  empty_matrix[10,2]=paste(DEG.HNSCC.used$X,collapse = ",")
  empty_matrix[11,2]=paste(DEG.Sarcoma_Osteosarcoma.used$X,collapse = ",")
  empty_matrix[12,2]=paste(DEG.ALL.used$X,collapse = ",")
  empty_matrix[13,2]=paste(DEG.Multiple_myeloma.used$X,collapse = ",")
  empty_matrix[14,2]=paste(DEG.CML.used$X,collapse = ",")
  empty_matrix[15,2]=paste(DEG.Neuroblastoma.used$X,collapse = ",")
  empty_matrix[16,2]=paste(DEg.CRC.used$X,collapse = ",")
  empty_matrix[17,2]=paste(DEG.Myeloproliferative_neoplasms.used$X,collapse = ",")
  empty_matrix[18,2]=paste(DEG.Ovarian_cancer.used$X,collapse = ",")
  return(empty_matrix)
}


