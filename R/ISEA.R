#' Convert mouse genes into corresponding human homologous genes
#'
#'
#' @param x A vector containing mouse genes
#'
#' @return Matrix containing human and mouse homologous genes
#' @export

mouse2human <- function(x){
  human = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")

  genes_cyno = biomaRt::getLDS(attributes = c("external_gene_name"),
                      filters = "external_gene_name",
                      values = x , mart = mouse,
                      attributesL = c("external_gene_name"),
                      martL = human, uniqueRows=T
  )
  colnames(genes_cyno) = c("Mouse.Gene","Human.Gene")
  return(genes_cyno)
}
#' Convert zebrafish genes into corresponding human homologous genes
#'
#'
#' @return Matrix containing human and zebrafish homologous genes
#' @export
#'
zebrafish2human <- function(x){
  human = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl" ,host = "https://dec2021.archive.ensembl.org/")
  zebrafish = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "drerio_gene_ensembl",host = "https://dec2021.archive.ensembl.org/")
  genes_cyno = biomaRt::getLDS(attributes = c("external_gene_name"),
                      filters = "external_gene_name",
                      values = x , mart = zebrafish,
                      attributesL = c("external_gene_name"),
                      martL = human, uniqueRows=T
  )
  colnames(genes_cyno) = c("Zebrafish.Gene","Human.Gene")
  return(genes_cyno)
}
#' Enrichment analysis of differentially expressed genes based on cancer signature
#'
#'
#' @param DEG.gene.select Differentially expressed genes list
#' @param sig.mtx Cancer signature matirx
#' @param species Specify which species to perform gene enrichment analysis. The default is human. You can also fill in mouse or zebrafish.
#'
#' @return Matrix containing ISEA result
#' @export
#'
ISEA <- function( DEG.gene.select,sig.mtx,species = "human") {
  if(species=="human"){
    DEG.gene.select=DEG.gene.select
  }else if(species=="mouse"){

    mtx = mouse2human(DEG.gene.select)
    DEG.gene.select = mtx$Human.Gene
  }else if(species=="zebrafish"){

    mtx = zebrafish2human(DEG.gene.select)
    DEG.gene.select = mtx$Human.Gene
  }
  human.gene  <- system.file("data", "human.tsv", package = "ISEA")
  human.gene = read.table(human.gene,sep = "\t",header = T)
  human.gene = human.gene$Gene.name
  Other.gene = setdiff(human.gene,DEG.gene.select)
  for ( i in 1:length(sig.mtx$gene_list)){
    signature =  sig.mtx[i,"gene_list"]
    signature = strsplit(signature, ",")
    signature =unlist(signature)
    DEG.inter = intersect(DEG.gene.select,signature)
    DEG.diff = setdiff(DEG.gene.select,signature)
    Other.inter = intersect(Other.gene,signature)
    Other.diff = setdiff(Other.gene,signature)
    deg.inter.count = length(DEG.inter)
    deg.diff.count = length(DEG.diff)
    other.inter.count =length(Other.inter)
    other.diff.count =length(Other.diff)
    test.mat = matrix(c(deg.inter.count,other.inter.count,deg.diff.count,other.diff.count),ncol=2,nrow=2)
    result = fisher.test(test.mat)
    sig.mtx[i,"p.value"] = result$p.value
    sig.mtx[i,"intersect_gene"] = deg.inter.count

  }

  sig.mtx <- sig.mtx[order(sig.mtx[, "p.value"], decreasing = FALSE), ]
  adjusted_p_values <- p.adjust(sig.mtx$p.value, method = "BH")
  sig.mtx$BH = adjusted_p_values
  adjusted_p_values <- p.adjust(sig.mtx$p.value, method = "bonferroni")
  sig.mtx$Bonferroni = adjusted_p_values
  sig.mtx$'log10_BH' = -log10(sig.mtx$BH)
  sig.mtx$'log10_bonferroni' = -log10(sig.mtx$Bonferroni)
  return(sig.mtx)
}

#' Visually display the results of ISEA enrichment analysis
#'
#'
#' @param sig.mtx A cancer signature matrix after ISEA enrichment analysis
#' @param method BH/bonferroni
#'
#' @return Visualization of ISEA enrichment analysis results
#' @export

ISEA.Plot <- function( sig.mtx,method="BH"){
  sig.mtx$cancer = factor(sig.mtx$cancer,levels = rev(sig.mtx$cancer))
  if(method == "BH"){
    # sig.mtx$cancer = factor(sig.mtx$cancer,levels = rev(sig.mtx$cancer))
    ggplot2::ggplot(data = sig.mtx, mapping = ggplot2::aes(x = cancer, y = log10_BH)) +
      ggplot2::geom_bar(stat = 'identity',fill="skyblue")+
      ggplot2::labs(title="Cancer enrichment", x="cancer type",y = "-log10(BH)")+
      ggplot2::coord_flip() +
      ggplot2::theme_bw()
  }else if(method == "bonferroni"){

    ggplot2::ggplot(data = sig.mtx, mapping = ggplot2::aes(x = cancer, y = log10_bonferroni)) +
      ggplot2::geom_bar(stat = 'identity',fill="skyblue")+
      ggplot2::labs(title="Cancer enrichment", x="cancer type",y = "-log10(bonferroni)")+
      ggplot2::coord_flip() +
      ggplot2::theme_bw()
  }


}




