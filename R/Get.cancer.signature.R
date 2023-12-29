

#' get cancer signature from DEG matirx
#'
#' @param adj.pvalue Set a adjusted p-value to screen for cancer signatures
#' @param log2fc Set up a log2fc to filter cancer signatures
#' @param max.num Set the upper line value for the number of genes in the cancer signature
#'
#' @return a dataframe with cancer signature information
#' @export
#'
#'
Get.cancer.signature <- function(adj.pvalue,log2fc,max.num=FALSE) {
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
  DEG.ALL.used = DEG.ALL[DEG.ALL$avg_log2FC>log2fc,]
  DEG.ALL.used = DEG.ALL.used[DEG.ALL.used$p_val_adj<adj.pvalue,]
  DEG.Multiple_myeloma.used =DEG.Multiple_myeloma[DEG.Multiple_myeloma$avg_log2FC>log2fc,]
  DEG.Multiple_myeloma.used = DEG.Multiple_myeloma.used[DEG.Multiple_myeloma.used$p_val_adj<adj.pvalue,]
  DEG.CML.used = DEG.CML[DEG.CML$avg_log2FC>log2fc,]
  DEG.CML.used = DEG.CML.used[DEG.CML.used$p_val_adj<adj.pvalue,]
  DEG.Neuroblastoma.used = DEG.Neuroblastoma[DEG.Neuroblastoma$avg_log2FC>log2fc,]
  DEG.Neuroblastoma.used = DEG.Neuroblastoma.used[DEG.Neuroblastoma.used$p_val_adj<adj.pvalue,]
  DEg.CRC.used = DEg.CRC[DEg.CRC$avg_log2FC>log2fc,]
  DEg.CRC.used =DEg.CRC.used[DEg.CRC.used$p_val_adj<adj.pvalue,]
  DEG.Myeloproliferative_neoplasms.used = DEG.Myeloproliferative_neoplasms[DEG.Myeloproliferative_neoplasms$avg_log2FC>log2fc,]
  DEG.Myeloproliferative_neoplasms.used = DEG.Myeloproliferative_neoplasms.used[DEG.Myeloproliferative_neoplasms.used$p_val_adj<adj.pvalue,]
  DEG.Ovarian_cancer.used  = DEG.Ovarian_cancer[DEG.Ovarian_cancer$avg_log2FC>log2fc,]
  DEG.Ovarian_cancer.used = DEG.Ovarian_cancer.used[DEG.Ovarian_cancer.used$p_val_adj<adj.pvalue,]
  DEG.Glioblastoma.used = DEG.Glioblastoma[DEG.Glioblastoma$avg_log2FC >log2fc,]
  DEG.Glioblastoma.used = DEG.Glioblastoma.used[DEG.Glioblastoma.used$p_val_adj< adj.pvalue,]
  DEG.H3K27M_gliomas.used = DEG.H3K27M_gliomas[DEG.H3K27M_gliomas$avg_log2FC >log2fc,]
  DEG.H3K27M_gliomas.used = DEG.H3K27M_gliomas.used[DEG.H3K27M_gliomas.used$p_val < adj.pvalue,]
  DEG.head_neck.used = DEG.head_neck[DEG.head_neck$avg_log2FC >log2fc,]
  DEG.head_neck.used =DEG.head_neck.used[DEG.head_neck.used$p_val < adj.pvalue,]
  DEG.Prostate.used = DEG.Prostate[DEG.Prostate$avg_log2FC>log2fc,]
  DEG.Prostate.used = DEG.Prostate.used[DEG.Prostate.used$p_val_adj<adj.pvalue,]
  DEG.lung.used = DEG.lung[DEG.lung$avg_log2FC>log2fc,]
  DEG.lung.used = DEG.lung.used[DEG.lung.used$p_val<adj.pvalue,]
  DEG.liver.used = DEG.liver[DEG.liver$avg_log2FC>log2fc,]
  DEG.liver.used = DEG.liver.used[DEG.liver.used$p_val_adj<adj.pvalue,]
  DEG.Kindey.used = DEG.Kindey[DEG.Kindey$avg_log2FC >log2fc,]
  DEG.Kindey.used = DEG.Kindey.used[DEG.Kindey.used$p_val_adj<adj.pvalue,]
  DEG.HNSCC.used = DEG.HNSCC[DEG.HNSCC$avg_log2FC>log2fc,]
  DEG.HNSCC.used = DEG.HNSCC.used[DEG.HNSCC.used$p_val_adj<adj.pvalue,]
  DEG.PDAC.used = DEG.PDAC[DEG.PDAC$avg_log2FC>log2fc,]
  DEG.PDAC.used = DEG.PDAC.used[DEG.PDAC.used$p_val_adj<adj.pvalue,]
  DEG.Sarcoma_Osteosarcoma.used = DEG.Sarcoma_Osteosarcoma[DEG.Sarcoma_Osteosarcoma$avg_log2FC>log2fc,]
  DEG.Sarcoma_Osteosarcoma.used = DEG.Sarcoma_Osteosarcoma.used[DEG.Sarcoma_Osteosarcoma.used$p_val_adj<adj.pvalue,]
  DEG.BREAST.used = DEG.BREAST[DEG.BREAST$avg_log2FC > log2fc,]
  DEG.BREAST.used = DEG.BREAST.used[DEG.BREAST.used$p_val < adj.pvalue,]

  if(max.num){
    DEG.ALL.used <- DEG.ALL.used[order(DEG.ALL.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.ALL.used)[1]>max.num){
      DEG.ALL.used=DEG.ALL.used[1:max.num,]
    }

    DEG.Multiple_myeloma.used <- DEG.Multiple_myeloma.used[order(DEG.Multiple_myeloma.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Multiple_myeloma.used)[1]>max.num){
      DEG.Multiple_myeloma.used=DEG.Multiple_myeloma.used[1:max.num,]
    }
    DEG.CML.used <- DEG.CML.used[order(DEG.CML.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.CML.used)[1]>max.num){
      DEG.CML.used=DEG.CML.used[1:max.num,]
    }
    DEG.Neuroblastoma.used <- DEG.Neuroblastoma.used[order(DEG.Neuroblastoma.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Neuroblastoma.used)[1]>max.num){
      DEG.Neuroblastoma.used=DEG.Neuroblastoma.used[1:max.num,]
    }
    DEg.CRC.used <- DEg.CRC.used[order(DEg.CRC.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEg.CRC.used)[1]>max.num){
      DEg.CRC.used=DEg.CRC.used[1:max.num,]
    }
    DEG.Myeloproliferative_neoplasms.used <- DEG.Myeloproliferative_neoplasms.used[order(DEG.Myeloproliferative_neoplasms.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Myeloproliferative_neoplasms.used)[1]>max.num){
      DEG.Myeloproliferative_neoplasms.used=DEG.Myeloproliferative_neoplasms.used[1:max.num,]
    }
    DEG.Ovarian_cancer.used <- DEG.Ovarian_cancer.used[order(DEG.Ovarian_cancer.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Ovarian_cancer.used)[1]>max.num){
      DEG.Ovarian_cancer.used=DEG.Ovarian_cancer.used[1:max.num,]
    }
    DEG.Glioblastoma.used <- DEG.Glioblastoma.used[order(DEG.Glioblastoma.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Glioblastoma.used)[1]>max.num){
      DEG.Glioblastoma.used=DEG.Glioblastoma.used[1:max.num,]
    }
    DEG.H3K27M_gliomas.used <- DEG.H3K27M_gliomas.used[order(DEG.H3K27M_gliomas.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.H3K27M_gliomas.used)[1]>max.num){
      DEG.H3K27M_gliomas.used=DEG.H3K27M_gliomas.used[1:max.num,]
    }
    DEG.head_neck.used <- DEG.head_neck.used[order(DEG.head_neck.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.head_neck.used)[1]>max.num){
      DEG.head_neck.used=DEG.head_neck.used[1:max.num,]
    }
    DEG.Prostate.used <- DEG.Prostate.used[order(DEG.Prostate.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Prostate.used)[1]>max.num){
      DEG.Prostate.used=DEG.Prostate.used[1:max.num,]
    }
    DEG.lung.used <- DEG.lung.used[order(DEG.lung.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.lung.used)[1]>max.num){
      DEG.lung.used=DEG.lung.used[1:max.num,]
    }
    DEG.liver.used <- DEG.liver.used[order(DEG.liver.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.liver.used)[1]>max.num){
      DEG.liver.used=DEG.liver.used[1:max.num,]
    }
    DEG.Kindey.used <- DEG.Kindey.used[order(DEG.Kindey.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Kindey.used)[1]>max.num){
      DEG.Kindey.used=DEG.Kindey.used[1:max.num,]
    }
    DEG.HNSCC.used <- DEG.HNSCC.used[order(DEG.HNSCC.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.HNSCC.used)[1]>max.num){
      DEG.HNSCC.used=DEG.HNSCC.used[1:max.num,]
    }
    DEG.PDAC.used <- DEG.PDAC.used[order(DEG.PDAC.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.PDAC.used)[1]>max.num){
      DEG.PDAC.used=DEG.PDAC.used[1:max.num,]
    }
    DEG.Sarcoma_Osteosarcoma.used <- DEG.Sarcoma_Osteosarcoma.used[order(DEG.Sarcoma_Osteosarcoma.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.Sarcoma_Osteosarcoma.used)[1]>max.num){
      DEG.Sarcoma_Osteosarcoma.used=DEG.Sarcoma_Osteosarcoma.used[1:max.num,]
    }
    DEG.BREAST.used <- DEG.BREAST.used[order(DEG.BREAST.used[, "p_val_adj"], decreasing = FALSE), ]
    if(dim(DEG.BREAST.used)[1]>max.num){
      DEG.BREAST.used=DEG.BREAST.used[1:max.num,]
    }
  }else{

  }
  ###创建gene matrix####

 cancer.list = c("Glioblastoma(GBM)",
                  "Gliomas","Breast Cancer(BRCA)",
                  "Pancreatic ductal adenocarcinoma(PDAC)","Nasopharyngeal carcinoma(NPC)",
                  "Prostate cancer(PRAD)",
                  "Lung adenocarcinoma(LUAD)",
                  "Liver cancer(HCC)",
                  "Kidney cancer(RCC)",
                  "Head and neck squamous cell carcinoma(HNSC)",
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

#' Visual display of the number of genes included in the cancer signature
#'
#' @param sig.mtx a matrix containing cancer signatures
#'
#' @export
#'
Signature.Plot<- function(sig.mtx) {
  for ( i in 1:length(sig.mtx$gene_list)){
    gene_list = sig.mtx[i,"gene_list"]
    gene_num = length(strsplit(gene_list, ",")[[1]])
    sig.mtx[i,"gene_num"] = gene_num
  }
  ggplot2::ggplot(data = sig.mtx, mapping = ggplot2::aes(x = cancer, y = gene_num)) +
    ggplot2::geom_bar(stat = 'identity',fill="skyblue")+
    ggplot2::labs(title="Cancer signature ", x="cancer type",y = "gene number")+
    ggplot2::coord_flip() +
    ggplot2::theme_bw()
}



