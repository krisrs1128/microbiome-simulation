#' Wrapper function for each DA_method
#' @param raw_data: an otu table with Otus on rows and Samples/Subjects on columns.
#' @param meta_data: meta data containing group information and other covariate information on subjects
#' @param DA_method: character, the DA method.
#' @return The output of the DA method (log fold change across groups, pvalue, adj.pvalue)
#' @export
#' @importFrom edgeR DGEList calcNormFactors
differential_analysis <- function(exper, DA_method) {
  
  raw_data <- assay(exper)
  meta_data <- colData(exper)
  
  if (DA_method == "ANCOMBC") {
    tseObj <- TreeSummarizedExperiment(SimpleList(counts = raw_data),
      colData = meta_data
    )
    out <- run_ancombc2(
      tseObj = tseObj,
      assay_name = "counts",
      fixed_formula = "bmi_group",
      p_adj_method = "BH"
    )
  } else if (DA_method == "LIMMA_voom") {
    dgeObj <- DGEList(raw_data, group = meta_data$bmi_group) |>
      calcNormFactors()
    out <- run_limma(
      dge = dgeObj,
      metaD = meta_data,
      formula = ~ bmi_group,
      p_adj_method = "BH"
    )
  } else if (DA_method == "DESeq2") {
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = raw_data,
      colData = meta_data,
      design = ~bmi_group
    )
    out <- run_DESeq(dds)
  }else if (DA_method == "wilcox-clr") {
    out <- run_wilcox(raw_data, clr = TRUE,  meta_data$bmi_group , pAdjustMethod = "BH")
  }
  out
}


####### ANCOM-BC2 method #######
#'
#' This function run ANCOM-BC2 method on a TreeSummarizedExperiment object.
#'
#' @param tseObj the input data. The data parameter should be either a phyloseq or a TreeSummarizedExperiment object,
#'             which consists of a feature table (microbial count table), a sample metadata table,
#'             a taxonomy table (optional), and a phylogenetic tree (optional).
#' @param assay_name character. Default is "counts". See ?SummarizedExperiment::assay for more details.
#' @param fixed_formula the character string expresses how the microbial absolute abundances for each taxon depend on the fixed effects in metadata.
#' @param p_adj_method the multiple correction procedure to be used.
#' @return a data frame with log fold changes, p_values, q_valuesh.
#'
#'

run_ancombc2 <- function(tseObj, assay_name = "counts", fixed_formula, p_adj_method = "BH") {
  output <- ancombc2(
    data = tseObj,
    assay_name = assay_name,
    fix_formula = fixed_formula,
    p_adj_method = p_adj_method
  )

  out_df <- data.frame(
    "log_FC(overweight)" = output$res$lfc_bmi_groupoverweight,
    "log_FC(lean)" = output$res$lfc_bmi_grouplean,
    "p_value(overweight)" = output$res$p_bmi_groupoverweight,
    "p_value(lean)" = output$res$p_bmi_grouplean,
    "q_value(overweight)" = output$res$q_bmi_groupoverweight,
    "q_value(lean)" = output$res$q_bmi_grouplean,
    row.names = rownames(output$feature_table), stringsAsFactors = F
  )

  return(out_df)
}
#######

########  LIMMA-VOOM method#######
#'
#' This function run LIMMA-VOOM method on a TreeSummarizedExperiment object.
#'
#' @param dge a numeric matrix containing raw counts, or an ExpressionSet containing raw counts, or a DGEList object.
#' @param metaD meta data containing group information and other covariate information on subjects
#' @param formula the character string expresses how the microbial absolute abundances for each taxon depend on the fixed effects in metadata.
#' @param p_adj_method the multiple correction procedure to be used.
#' @return a data frame with log fold changes, p_values, q_valuesh.
#'
#'
run_limma <- function(dge, metaD, formula, p_adj_method = "BH") {
  design <- model.matrix(formula, data = metaD)
  baysfit <- voom(dge, design) |> # account for differences in uncertainty
    lmFit(design) |>
    eBayes()
  res1 <- topTable(baysfit, number = Inf, coef = "bmi_groupoverweight", adjust.method = p_adj_method) |>
    rownames_to_column("ID")
  res2 <- topTable(baysfit, number = Inf, coef = "bmi_grouplean", adjust.method = p_adj_method) |>
    rownames_to_column("ID")


  out_df <- data.frame(
    "log_FC(overweight)" = res1$logFC,
    "log_FC(lean)" = res2$logFC,
    "p_value(overweight)" = res1$P.Value,
    "p_value(lean)" = res2$P.Value,
    "q_value(overweight)" = res1$adj.P.Val,
    "q_value(lean)" = res1$adj.P.Va,
    row.names = res1$ID, stringsAsFactors = F
  )

  return(out_df)
}
#######


########  DESeq method#######
#'
#' This function run DESeq method on a DESeq object.
#'
#' @param dds  a DESeqDataSet object.
#' @param sftype 	either "ratio", "poscounts", or "iterate" for the type of size factor estimation.
#' @param tidy 	whether to output the results table with rownames as a first column
#' @param format 	character, either "DataFrame", "GRanges", or "GRangesList"
#' @param pAdjustMethod the multiple correction procedure to be used.
#' @return a data frame with log fold changes, p_values, q_valuesh.
#'
#'
run_DESeq <- function(dds, sftype = "poscounts", tidy = TRUE, format = "DataFrame", pAdjustMethod = "BH") {
  dds_res <- DESeq2::DESeq(dds, sfType = sftype)

  res1 <- DESeq2::results(dds_res,
    tidy = tidy, format = format,
    pAdjustMethod = pAdjustMethod,
    contrast = c("bmi_group", "overweight", "obese")
  )
  res2 <- DESeq2::results(dds_res,
    tidy = tidy, format = format,
    pAdjustMethod = pAdjustMethod,
    contrast = c("bmi_group", "lean", "obese")
  )


  out_df <- data.frame(
    "log_FC(overweight)" = res1$log2FoldChange,
    "log_FC(lean)" = res2$log2FoldChange,
    "p_value(overweight)" = res1$pvalue,
    "p_value(lean)" = res2$pvalue,
    "q_value(overweight)" = res1$padj,
    "q_value(lean)" = res2$padj,
    row.names = res1$row, stringsAsFactors = F
  )

  return(out_df)
}

########   Wilcox method#######
#'
#' his function run Wilcoxon Rank Sum test on otu table.
#'
#' @param rawData  a numeric matrix containing raw counts 
#' @param clr 	 whether to perform clr transformation. Default=F
#' @param group_labels character vector of the same length as the columns of 'raw_data'. the labels of each subject.
#' @param pAdjustMethod the multiple correction procedure to be used.
#' @return a data frame with log fold changes, p_values, q_valuesh.
#'
#'
run_wilcox <- function(rawData, clr = FALSE, group_labels , pAdjustMethod = "BH") {
  
  if(clr){
    rawData[rawData == 0] <- 1e-6
    clr_counts <- apply(rawData, 1, compositions::clr)
    counts <- t(clr_counts)  
  }else{
    counts <- rawData
  }
  
  p_values <- numeric(nrow(counts))
  
  # Perform Kruskal-Wallis test for each row
  for (i in 1:nrow(counts)) {
    test_result <- kruskal.test(counts[i, ] ~ group_labels)
    p_values[i] <- test_result$p.value
  }
  q_value <- p.adjust(p_values, method=pAdjustMethod)
  
  
  out_df <- data.frame(
    "p_value(bmi_group)" = p_values,
    "q_value(bmi_group)" = q_value,
    row.names = rownames(rawData)
  )
  
  return(out_df)
}


#######

#' FDR and Power for Differential Analysis
#' @export
da_metrics <- function(results, null, level = 0.1, focus_col = ncol(results)) {
  flagged <- results[results[[focus_col]] < level, ] |>
    rownames()
  
  nonnull <- setdiff(rownames(results), null)
  data.frame(
    metric = c("FDR", "power"),
    value = c(length(intersect(null, flagged)) / max(1, length(flagged)), length(intersect(nonnull, flagged)) / length(nonnull))
  )
}
