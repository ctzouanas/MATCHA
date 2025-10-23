#' Multi/tandem-omic workflow - Map a gene program of interest to co-accessible peaks
#' @param rna.atac.combo.list.in List of all needed datasets, formatted as described in Github example notebook
#' @param program.name.in Name of desired gene program
#' @param genes.in Genes comprising the program
#' @param save.folder.in Folder for saving outputs
#' @param seqlevelsStyle.in Style of genomic chromosome names
#' @import Matrix
#' @export
genestopeaks.matcha <- function(
    rna.atac.combo.list.in,
    program.name.in,
    genes.in,
    save.folder.in,
    seqlevelsStyle.in = "UCSC"){

  dataset.name.in = rna.atac.combo.list.in[["dataset.name"]]
  hep.atac.obj.in = rna.atac.combo.list.in[["clean.atac.obj"]]
  conns.in = rna.atac.combo.list.in[["conns"]]
  search.window = rna.atac.combo.list.in[["search.window"]]
  strength.filter = rna.atac.combo.list.in[["strength.filter"]]
  ensdb.in = GetGRangesFromEnsDb(rna.atac.combo.list.in[["ensdb.chosen"]])
  gene.name.conversion.df <- rna.atac.combo.list.in[["gene.name.conversion"]]
  if(!is.null(gene.name.conversion.df)){
    genes.in <- gene.name.conversion.df %>%
      dplyr::filter(!is.na(name.universal) & !is.na(name.datasetspecific)) %>%
      dplyr::filter(name.universal %in% genes.in) %>%
      dplyr::pull(name.datasetspecific) %>% unique()
  }

  annotations <- Signac::Annotation(hep.atac.obj.in)
  gene.ranges <- ensdb.in
  GenomeInfoDb::seqlevelsStyle(gene.ranges) <- seqlevelsStyle.in

  hep.closestfeature <- Signac::ClosestFeature(hep.atac.obj.in,
                                               regions = rownames(hep.atac.obj.in))

  peaks.in.genes.of.int <- hep.closestfeature %>%
    dplyr::filter(gene_name %in% genes.in) %>%
    dplyr::filter(distance == 0)

  conns.filter.gene.of.int <- conns.in %>%
    dplyr::filter(coaccess > strength.filter) %>%
    dplyr::filter(Peak1 %in% peaks.in.genes.of.int$query_region | Peak2 %in% peaks.in.genes.of.int$query_region) %>%
    unique()

  peaks.of.int.final <- unique(c(as.character(conns.filter.gene.of.int$Peak1), as.character(conns.filter.gene.of.int$Peak2)))

  annotations.subset <- annotations[annotations$gene_name %in% genes.in]

  hep.ranges <- hep.atac.obj.in@assays[[rna.atac.combo.list.in[["atac.assay.name"]]]]@ranges

  hep.ranges.within.window <- Signac::findOverlaps(query = hep.ranges,
                                                   subject = annotations.subset,
                                                   maxgap = search.window,
                                                   ignore.strand = TRUE)

  hep.overlap.peak.indices <- hep.ranges.within.window@from %>% unique()

  hep.peaks.in.proximity <- rownames(hep.atac.obj.in)[hep.overlap.peak.indices]

  peaks.of.int.final <- peaks.of.int.final[peaks.of.int.final %in% hep.peaks.in.proximity]

  peaks.export.jj <- peaks.of.int.final %>% data.frame() %>%
    tidyr::separate(col = ".", into = c("chr", "start", "end"), sep = "-")

  motif.matrix <- Signac::GetMotifData(object = hep.atac.obj.in, slot = "data")
  peak.matrix <- Seurat::GetAssayData(object = hep.atac.obj.in, slot = "counts")
  idx.keep <- Matrix::rowSums(x = peak.matrix) > 0
  peak.matrix <- peak.matrix[idx.keep, , drop = FALSE]
  motif.matrix <- motif.matrix[idx.keep, , drop = FALSE]
  peak.ranges <- Signac::granges(x = hep.atac.obj.in)
  peak.ranges <- peak.ranges[idx.keep]
  chromvar.obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak.matrix),
    rowRanges = peak.ranges
  )
  chromvar.obj <- chromVAR::addGCBias(
    object = chromvar.obj,
    genome = BSgenome.Mmusculus.UCSC.mm10
  )
  bg <- chromVAR::getBackgroundPeaks(
    object = chromvar.obj
  )

  bg.subset.indices <- bg[rownames(hep.atac.obj.in) %in% peaks.of.int.final, ] %>% as.vector() %>% unique()
  background.export.jj <- rownames(hep.atac.obj.in)[bg.subset.indices] %>%
    data.frame() %>%
    tidyr::separate(col = ".", into = c("chr", "start", "end"), sep = "-")

  write.table(x = peaks.export.jj,
              file = paste0(save.folder.in, "/",
                            dataset.name.in, "_", program.name.in, "_", search.window/1000, "kbWindow_", strength.filter, "CoaccessFilter_Interest.bed"),
              col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE)

  write.table(x = background.export.jj,
              file = paste0(save.folder.in, "/",
                            dataset.name.in, "_", program.name.in, "_", search.window/1000, "kbWindow_", strength.filter, "CoaccessFilter_GCandAccessibilityMatched.bed"),
              col.names = FALSE, row.names = FALSE,
              sep = "\t", quote = FALSE)
}

#' Multi/tandem-omic workflow - Quantify accessibility of transcription factor motifs within program-linked co-accessible peaks
#' @param rna.atac.combo.list.in List of all needed datasets, formatted as described in Github example notebook
#' @param save.folder.in Folder in which the output of genestopeaks.matcha was saved
#' @param program.name.in Name of desired gene program
#' @import Matrix
#' @return Signac object with a new chromVAR assay corresponding to transcription factor motif accessibility within program-linked co-accessible peaks
#' @export
peakstoTFs.matcha <- function(
    rna.atac.combo.list.in,
    save.folder.in,
    program.name.in){
  peakset.choice <- program.name.in
  dataset.name.in = rna.atac.combo.list.in[["dataset.name"]]
  hep.atac.obj.in = rna.atac.combo.list.in[["clean.atac.obj"]]
  search.window = rna.atac.combo.list.in[["search.window"]]
  strength.filter = rna.atac.combo.list.in[["strength.filter"]]
  genome.in = rna.atac.combo.list.in[["bsgenome.chosen"]]
  gene.name.conversion.df <- rna.atac.combo.list.in[["gene.name.conversion"]]
  atac.assay.name = rna.atac.combo.list[["atac.assay.name"]]
  if(!is.null(gene.name.conversion.df)){
    genes.in <- gene.name.conversion.df %>%
      dplyr::filter(!is.na(name.universal) & !is.na(name.datasetspecific)) %>%
      dplyr::filter(name.universal %in% genes.in) %>%
      dplyr::pull(name.datasetspecific) %>% unique()
  }

  peaks.of.int.ii <- read.table(
    file = paste0(save.folder.in, "/",
                  dataset.name.in, "_", program.name.in, "_", search.window/1000, "kbWindow_", strength.filter, "CoaccessFilter_Interest.bed"),
    sep = "\t", header = FALSE) %>%
    dplyr::mutate(full_peakname = paste0(V1, "-", V2, "-", V3)) %>% pull(full_peakname)

  Seurat::DefaultAssay(hep.atac.obj.in) <- atac.assay.name
  motif.matrix <- Signac::GetMotifData(object = hep.atac.obj.in)
  motif.matrix.mod <- Matrix::Matrix(data = 0,
                                     nrow = nrow(motif.matrix),
                                     ncol = ncol(motif.matrix),
                                     sparse = TRUE,
                                     dimnames = list(rownames(motif.matrix), colnames(motif.matrix)))
  peaks.to.keep <- rownames(motif.matrix.mod)[rownames(motif.matrix.mod) %in% peaks.of.int.ii]
  motif.matrix.mod[peaks.to.keep, ] <- motif.matrix[peaks.to.keep, ]

  hep.atac.obj.in <- Signac::RunChromVAR(
    object = hep.atac.obj.in,
    genome = genome.in,
    motif.matrix = motif.matrix.mod,
    new.assay.name = paste0("chromVAR.", peakset.choice)
  )
  hep.atac.obj.in@assays[[paste0("chromVAR.", peakset.choice)]]@meta.features <- hep.atac.obj.in@misc$chromvar.meta.features

  return(hep.atac.obj.in)
}

#' Multi/tandem-omic workflow - Identify transcription factors whose expression and motif accessibility within program-linked co-accessible peaks associates with program expression
#' @param rna.atac.combo.list.in List of all needed datasets, formatted as described in Github example notebook
#' @param save.folder.in Folder in which the output of genestopeaks.matcha was saved
#' @param program.name.in Name of desired gene program
#' @param genes.in Genes comprising the gene program
#' @import Matrix
#' @return Dataframe of prioritized transcription factors
#' @export
TFstoRankedTFs.atac.rna.matcha <- function(
    rna.atac.combo.list.in,
    save.folder.in,
    program.name.in,
    genes.in){

  dataset.name.in = rna.atac.combo.list.in[["dataset.name"]]
  hep.atac.obj.in <- rna.atac.combo.list.in[["clean.atac.obj"]]
  atac.sample.field <- rna.atac.combo.list.in[["atac.sample.field"]]
  hep.rna.obj.in <- rna.atac.combo.list.in[["clean.rna.obj"]]
  rna.sample.field <- rna.atac.combo.list.in[["rna.sample.field"]]
  gene.name.conversion.df <- rna.atac.combo.list.in[["gene.name.conversion"]]
  motif.name.conversion.df <- rna.atac.combo.list.in[["motif.name.conversion"]]
  if(!is.null(gene.name.conversion.df)){
    genes.in <- gene.name.conversion.df %>%
      dplyr::filter(!is.na(name.universal) & !is.na(name.datasetspecific)) %>%
      dplyr::filter(name.universal %in% genes.in) %>%
      pull(name.datasetspecific) %>% unique()
  }

  hep.atac.obj.in$cell.barcode <- colnames(hep.atac.obj.in)
  hep.atac.obj.in$sample.name <- hep.atac.obj.in@meta.data[, atac.sample.field]
  hep.rna.obj.in$cell.barcode <- colnames(hep.rna.obj.in)
  hep.rna.obj.in$sample.name <- hep.rna.obj.in@meta.data[, rna.sample.field]

  hep.atac.chromvar.scores.long <- hep.atac.obj.in@assays[[paste0("chromVAR.", program.name.in)]]@data %>% data.frame() %>%
    rownames_to_column("TF.num") %>%
    tidyr::pivot_longer(!TF.num, names_to = "cell.barcode", values_to = "chromVAR.score") %>%
    dplyr::mutate(cell.barcode = gsub(pattern = "\\.", replacement = "-", x = .$cell.barcode))
  hep.atac.chromvar.scores.long <- hep.atac.chromvar.scores.long %>%
    dplyr::filter(!is.na(chromVAR.score)) %>%
    dplyr::inner_join(x = .,
                      y = hep.atac.obj.in@meta.data[, c("cell.barcode", "sample.name")],
                      by = "cell.barcode")
  hep.atac.chromvar.summarise <- hep.atac.chromvar.scores.long %>%
    dplyr::group_by(sample.name, TF.num) %>%
    dplyr::summarise(mean.TF.chromVAR.score = mean(chromVAR.score)) %>% ungroup()
  hep.atac.chromvar.summarise.wide <- hep.atac.chromvar.summarise %>%
    tidyr::pivot_wider(names_from = TF.num, values_from = mean.TF.chromVAR.score, values_fill = 0) %>%
    tibble::column_to_rownames("sample.name") %>%
    as.matrix()

  hep.rna.obj.in <- Seurat::AddModuleScore(object = hep.rna.obj.in, features = list(genes.in), name = program.name.in,
                                           assay = rna.atac.combo.list.in[["rna.assay.name"]],
                                           seed = 42, nbin = 24, ctrl = 50)
  hep.rna.module.df <- hep.rna.obj.in@meta.data[, c("cell.barcode", "sample.name", paste0(program.name.in, "1"))] %>%
    dplyr::rename(module.choice = 3) %>%
    dplyr::group_by(sample.name) %>%
    dplyr::summarise(mean.module.score = mean(module.choice)) %>%
    tibble::column_to_rownames("sample.name") %>% as.matrix()

  hep.rna.tfavgexpr <- Seurat::AverageExpression(hep.rna.obj.in, assays = rna.atac.combo.list.in[["rna.assay.name"]],
                                                 features = unique(motif.name.conversion.df$name.datasetspecific),
                                                 group.by = "sample.name")[[1]] %>% t()

  if((sum(rownames(hep.atac.chromvar.summarise.wide) == rownames(hep.rna.module.df)) == nrow(hep.atac.chromvar.summarise.wide)) &
     nrow(hep.atac.chromvar.summarise.wide) == nrow(hep.rna.module.df)){
    cor.TFchromvar.modulescore.mat <- cor(x = hep.atac.chromvar.summarise.wide,
                                          y = hep.rna.module.df, method = "spearman") %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "motif.name") %>%
      dplyr::relocate(motif.name) %>%
      dplyr::rename(TFchromVAR.ExpressionModule.Cor = 2) %>%
      dplyr::inner_join(.,
                        y = hep.atac.obj.in@misc$chromvar.meta.features,
                        by = c("motif.name" = "orig.colnames")) %>%
      dplyr::arrange(-TFchromVAR.ExpressionModule.Cor) %>% dplyr::relocate(motif.name)
  } else{
    stop("Mismatch on sample names!")
  }

  if((sum(rownames(hep.rna.tfavgexpr) == rownames(hep.rna.module.df)) == nrow(hep.rna.tfavgexpr)) &
     nrow(hep.rna.tfavgexpr) == nrow(hep.rna.module.df)){
    cor.TFexpress.modulescore.mat <- cor(x = as.matrix(hep.rna.tfavgexpr[, colSums(hep.rna.tfavgexpr != 0) > 0]),
                                         y = hep.rna.module.df,
                                         method = "spearman") %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "name.datasetspecific") %>%
      dplyr::relocate(name.datasetspecific) %>%
      dplyr::rename(TFexpress.ExpressionModule.Cor = 2) %>%
      dplyr::left_join(x = .,
                       y = motif.name.conversion.df,
                       by = c("name.datasetspecific")) %>%
      dplyr::arrange(-TFexpress.ExpressionModule.Cor) %>%
      dplyr::select(name.datasetspecific, TFexpress.ExpressionModule.Cor) %>%
      dplyr::group_by(name.datasetspecific) %>%
      dplyr::summarise(TFexpress.ExpressionModule.Cor = mean(TFexpress.ExpressionModule.Cor, na.rm = TRUE)) %>% ungroup()
  } else{
    stop("Mismatch on sample names!")
  }

  # Need to futz with ranking depending on whether gene set is going up/down - that'll affect whether agreement means positive or negative ranking of Govaere dataset
  TF.prioritization.chromVARandExpression <-
    dplyr::inner_join(x = motif.name.conversion.df,
                      y = cor.TFchromvar.modulescore.mat %>% dplyr::select(motif.name, TFchromVAR.ExpressionModule.Cor),
                      by = c("motif.name")) %>%
    dplyr::inner_join(x = .,
                      y = cor.TFexpress.modulescore.mat,
                      by = c("name.datasetspecific")) %>%
    dplyr::mutate(
      TFchromVAR.ExpressionModule.Cor.Scale = (TFchromVAR.ExpressionModule.Cor - min(TFchromVAR.ExpressionModule.Cor, na.rm = TRUE)) /
        (max(TFchromVAR.ExpressionModule.Cor, na.rm = TRUE) - min(TFchromVAR.ExpressionModule.Cor, na.rm = TRUE)),
      TFexpress.ExpressionModule.Cor.Scale = (TFexpress.ExpressionModule.Cor - min(TFexpress.ExpressionModule.Cor, na.rm = TRUE)) /
        (max(TFexpress.ExpressionModule.Cor, na.rm = TRUE) - min(TFexpress.ExpressionModule.Cor, na.rm = TRUE))) %>%
    dplyr::mutate(overall.scale = 0.5*(TFchromVAR.ExpressionModule.Cor.Scale + TFexpress.ExpressionModule.Cor.Scale)) %>%
    dplyr::arrange(-overall.scale) %>% dplyr::relocate(overall.scale, .after = 3)

  write.table(x = TF.prioritization.chromVARandExpression,
              file = paste0(save.folder.in, "/", dataset.name.in, "_", program.name.in, "_TFPrioritization_ATACandRNA.tsv"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  pdf(file = paste0(save.folder.in, "/",
                    dataset.name.in, "_", program.name.in, "_TFPrioritization_ATACandRNA.pdf"), width = 16, height = 16)
  p1 <- ggplot(TF.prioritization.chromVARandExpression, aes(x = TFexpress.ExpressionModule.Cor, y = TFchromVAR.ExpressionModule.Cor, label = name.universal)) +
    xlab("Correlation between\nTF Expression and\nGene Module Expression") +
    ylab("Correlation between\nAccessibility at TF Binding Sites and\nGene Module Expression") +
    ggtitle(program.name.in) +
    geom_hline(yintercept = 0) + geom_vline(xintercept = 0) +
    geom_point() + geom_text_repel(max.overlaps = 10) + theme_classic() +
    theme(text = element_text(size = 16), plot.title = element_text(hjust = 0.5))

  plot(p1)
  dev.off()

  return(TF.prioritization.chromVARandExpression)

}

#' Transcriptomic-only workflow - Identify transcription factors whose expression associates with program expression
#' @param rna.only.list.in List of all needed datasets, formatted as described in Github example notebook
#' @param save.folder.in Folder in which the output of genestopeaks.matcha was saved
#' @param program.name.in Name of desired gene program
#' @param genes.in Genes comprising the gene program
#' @import Matrix
#' @return Dataframe of prioritized transcription factors
#' @export
TFstoRankedTFs.rna.only.matcha <- function(
    rna.only.list.in,
    save.folder.in,
    program.name.in,
    genes.in){

  dataset.name.in = rna.only.list.in[["dataset.name"]]
  hep.rna.obj.in <- rna.only.list.in[["clean.rna.obj"]]
  rna.sample.field <- rna.only.list.in[["rna.sample.field"]]
  gene.name.conversion.df <- rna.only.list.in[["gene.name.conversion"]]
  motif.name.conversion.df <- rna.only.list.in[["motif.name.conversion"]]
  if(!is.null(gene.name.conversion.df)){
    genes.in <- gene.name.conversion.df %>%
      dplyr::filter(!is.na(name.universal) & !is.na(name.datasetspecific)) %>%
      dplyr::filter(name.universal %in% genes.in) %>%
      dplyr::pull(name.datasetspecific) %>% unique()
  }

  hep.rna.obj.in$cell.barcode <- colnames(hep.rna.obj.in)
  hep.rna.obj.in$sample.name <- hep.rna.obj.in@meta.data[, rna.sample.field]

  hep.rna.obj.in <- Seurat::AddModuleScore(object = hep.rna.obj.in, features = list(genes.in), name = program.name.in,
                                           assay = rna.only.list.in[["rna.assay.name"]],
                                           seed = 42, nbin = 24, ctrl = 50)
  hep.rna.module.df <- hep.rna.obj.in@meta.data[, c("cell.barcode", "sample.name", paste0(program.name.in, "1"))] %>%
    dplyr::rename(module.choice = 3) %>%
    dplyr::group_by(sample.name) %>%
    dplyr::summarise(mean.module.score = mean(module.choice)) %>%
    tibble::column_to_rownames("sample.name") %>% as.matrix()

  hep.rna.tfavgexpr <- Seurat::AverageExpression(hep.rna.obj.in, assays = rna.only.list.in[["rna.assay.name"]],
                                                 features = unique(motif.name.conversion.df$name.datasetspecific),
                                                 group.by = "sample.name")[[1]] %>% t()
  rownames(hep.rna.tfavgexpr) <- janitor::make_clean_names(rownames(hep.rna.tfavgexpr))

  if((sum(rownames(hep.rna.tfavgexpr) == rownames(hep.rna.module.df)) == nrow(hep.rna.tfavgexpr)) &
     nrow(hep.rna.tfavgexpr) == nrow(hep.rna.module.df)){
    cor.TFexpress.modulescore.mat <- cor(x = as.matrix(hep.rna.tfavgexpr[, colSums(hep.rna.tfavgexpr != 0) > 0]),
                                         y = hep.rna.module.df,
                                         method = "spearman") %>%
      data.frame() %>%
      tibble::rownames_to_column(var = "name.datasetspecific") %>%
      dplyr::relocate(name.datasetspecific) %>%
      dplyr::rename(TFexpress.ExpressionModule.Cor = 2) %>%
      dplyr::left_join(x = .,
                       y = motif.name.conversion.df,
                       by = c("name.datasetspecific")) %>%
      dplyr::arrange(-TFexpress.ExpressionModule.Cor) %>%
      dplyr::select(name.datasetspecific, TFexpress.ExpressionModule.Cor) %>%
      dplyr::group_by(name.datasetspecific) %>%
      dplyr::summarise(TFexpress.ExpressionModule.Cor = mean(TFexpress.ExpressionModule.Cor, na.rm = TRUE)) %>% ungroup()
  } else{
    stop("Mismatch on sample names!")
  }

  # Need to futz with ranking depending on whether gene set is going up/down - that'll affect whether agreement means positive or negative ranking of Govaere dataset
  TF.prioritization.expression <-
    dplyr::inner_join(x = motif.name.conversion.df,
                      y = cor.TFexpress.modulescore.mat,
                      by = c("name.datasetspecific")) %>%
    dplyr::mutate(
      TFexpress.ExpressionModule.Cor.Scale = (TFexpress.ExpressionModule.Cor - min(TFexpress.ExpressionModule.Cor, na.rm = TRUE)) /
        (max(TFexpress.ExpressionModule.Cor, na.rm = TRUE) - min(TFexpress.ExpressionModule.Cor, na.rm = TRUE))) %>%
    dplyr::mutate(overall.scale = TFexpress.ExpressionModule.Cor.Scale) %>%
    dplyr::arrange(-overall.scale) %>% dplyr::relocate(overall.scale, .after = 3)

  write.table(x = TF.prioritization.expression,
              file = paste0(save.folder.in, "/", dataset.name.in, "_", program.name.in, "_TFPrioritization_RNAonly.tsv"),
              row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")

  return(TF.prioritization.expression)

}

#' For a single gene program, prioritize consensus transcription factors that drive the program across all datasets
#' @param folder.in Folder in which outputs of all datasets are saved
#' @param name.in Name of gene program
#' @param plot.width Desired width of output plot
#' @param plot.height Desired height of output plot
#' @param plot.font.size Desired font size of output plot
#' @import Matrix
#' @export
singleprogram.crossdataset.TFprioritize.matcha <- function(
    folder.in,
    name.in,
    plot.width = 6,
    plot.height = 6,
    plot.font.size = 12){

  dataset.folders.ii <- list.dirs(folder.in, full.names = FALSE)
  dataset.folders.ii <- dataset.folders.ii[dataset.folders.ii != ""]
  datasets.ii <- gsub(pattern = name.ii, replacement = "", x = dataset.folders.ii)
  datasets.ii <- gsub(pattern = "__[A-z0-9]*$", replacement = "", x = datasets.ii)

  collapsetfs.compile.df <- data.frame()
  for(jj in 1:length(dataset.folders.ii)){
    dataset.ii.jj <- datasets.ii[[jj]]
    dataset.path.ii.jj <- paste0(folder.in, "/", dataset.folders.ii[[jj]])
    files.ii.jj <- list.files(dataset.path.ii.jj)
    chosen.file.ii.jj <- files.ii.jj[grepl(pattern = "TFPrioritization_[A-z0-9]*\\.tsv$", x = files.ii.jj)]
    df.ii.jj <- read.table(
      file = paste0(dataset.path.ii.jj, "/", chosen.file.ii.jj),
      sep = "\t",
      header = TRUE
    ) %>%
      dplyr::select(!overall.scale) %>%
      tidyr::pivot_longer(!c(motif.name, name.universal, name.datasetspecific), names_to = "Comparison", values_to = "Value") %>%
      dplyr::mutate(Comparison.Type = if_else(
        condition = grepl(pattern = "\\.Scale$", x = Comparison),
        true = "Scale", false = "Cor"
      )) %>%
      dplyr::relocate(Comparison.Type, .after = Comparison)
    df.collapsetfs.ii.jj <- df.ii.jj %>%
      dplyr::select(!c(motif.name, name.datasetspecific)) %>%
      dplyr::group_by(name.universal, Comparison, Comparison.Type) %>%
      dplyr::summarise(Value = mean(Value)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dataset = dataset.ii.jj) %>%
      dplyr::relocate(dataset, .after = name.universal)
    collapsetfs.compile.df <- rbind(collapsetfs.compile.df, df.collapsetfs.ii.jj) %>%
      dplyr::arrange(name.universal)
  }

  prioritize.compile.df <- collapsetfs.compile.df %>%
    dplyr::filter(Comparison.Type == "Scale") %>%
    dplyr::group_by(name.universal) %>%
    dplyr::summarise(mean.Value = mean(Value), n = n()) %>%
    dplyr::filter(n == max(n)) %>%
    dplyr::arrange(-mean.Value)
  extreme.tfs <- c(
    prioritize.compile.df %>% dplyr::slice_max(order_by = mean.Value, n = 10) %>% dplyr::pull(name.universal),
    prioritize.compile.df %>% dplyr::slice_min(order_by = mean.Value, n = 10) %>% dplyr::pull(name.universal)
  )

  plot.compile.df <- collapsetfs.compile.df %>%
    dplyr::filter(name.universal %in% extreme.tfs) %>%
    dplyr::filter(Comparison.Type == "Cor")
  plot.summary.df <- plot.compile.df %>%
    dplyr::group_by(name.universal) %>%
    dplyr::summarise(Value = mean(Value)) %>%
    dplyr::mutate(dataset = "Overall.Summary", Comparison = "Overall.Summary", Comparison.Type = "Cor") %>%
    dplyr::relocate(Value, .after = Comparison.Type) %>%
    dplyr::arrange(-Value)
  plot.combo.df <- rbind(plot.compile.df, plot.summary.df) %>%
    dplyr::mutate(name.universal = factor(name.universal,
                                          levels = rev(plot.summary.df$name.universal))) %>%
    dplyr::arrange(name.universal) %>%
    dplyr::mutate(dataset.comparison = paste0(dataset, ".", Comparison)) %>%
    dplyr::mutate(dataset.comparison = if_else(
      condition = dataset.comparison == "Overall.Summary.Overall.Summary",
      true = "Overall.Summary", false = dataset.comparison
    )) %>%
    dplyr::mutate(is.summary = if_else(
      condition = dataset.comparison == "Overall.Summary",
      true = "Overall.Summary", false = "Not"
    ))

  n.comparisons = length(unique(plot.combo.df$dataset.comparison))-1
  plot.colors <- paletteer::paletteer_c("ggthemes::Sunset-Sunrise Diverging", n.comparisons) %>% as.character()
  names(plot.colors) <- plot.combo.df %>%
    dplyr::filter(dataset.comparison != "Overall.Summary") %>%
    dplyr::pull(dataset.comparison) %>%
    unique()
  plot.colors["Overall.Summary"] <- "#555555"

  p.tf.prioritization.dotplot <-
    ggplot2::ggplot(plot.combo.df, ggplot2::aes(x = Value, y = name.universal, color = dataset.comparison, shape = dataset.comparison)) +
    ggplot2::geom_vline(xintercept = 0, color = "#dddddd") +
    ggplot2::scale_color_manual(values = plot.colors) +
    ggplot2::geom_point(aes(shape = is.summary, size = is.summary)) +
    ggplot2::scale_shape(guide = "none") +
    ggplot2::scale_size_manual(values = c("Overall.Summary" = 2, "Not"=1, 1), guide = "none") +
    ggplot2::xlab("Spearman Correlation") +
    ggplot2::theme_classic() +
    ggplot2::ggtitle(name.ii) +
    ggplot2::theme(legend.position = "bottom", axis.title.y = element_blank(),
                   plot.title = ggplot2::element_text(hjust = 0.5),
                   legend.box = "vertical", legend.justification = "center", legend.direction = "vertical",
                   text = ggplot2::element_text(size = plot.font.size))

  pdf(paste0(folder.in, "/", name.ii, "_CrossDatasetTFPrioritization.pdf"), width = plot.width, height = plot.height)
  plot(p.tf.prioritization.dotplot)
  dev.off()

  df.export <- list(
    "CrossDataset.Prioritization" = prioritize.compile.df %>% dplyr::select(!n) %>% dplyr::rename(mean.Scaled.Value = mean.Value),
    "Cor.Value" = collapsetfs.compile.df %>% dplyr::filter(Comparison.Type == "Cor"),
    "Scaled.Value" = collapsetfs.compile.df %>% dplyr::filter(Comparison.Type == "Scale")
  )

  writexl::write_xlsx(
    x = df.export,
    path = paste0(folder.in, "/", name.ii, "_CrossDatasetTFPrioritization.xlsx")
  )
}

#' Across gene programs, prioritize consensus transcription factors that drive multiple programs across all datasets
#' @param folder.in Folder in which outputs of all datasets are saved
#' @param gene.list.in List object containing all gene programs, formatted as in example notebook
#' @param n.tfs.per.program Top and bottom n transcription factors to retain for each program
#' @param min.tf.degree Minimum number of programs to which a transcription factor has to be linked in order to be retained
#' @param graph.layout Layout algorithm for plotting transcription factor - gene program network. Can be any layout compatible with ggraph
#' @param plot.width Desired width of output plot
#' @param plot.height Desired height of output plot
#' @param plot.font.size Desired font size of output plot
#' @import Matrix
#' @export
multiprogram.crossdataset.TFprioritize.matcha <- function(
    folder.in,
    gene.list.in,
    n.tfs.per.program = 10,
    min.tf.degree = 2,
    graph.layout = "fr",
    plot.width = 8,
    plot.height = 8,
    plot.font.size = 12
){
  extreme.tfs <- data.frame()
  n.tfs <- n.tfs.per.program
  degree.display.label <- min.tf.degree

  for(ii in 1:length(names(gene.list.in))){
    name.ii <- names(gene.list.in)[[ii]]
    priority.df.ii <-
      readxl::read_excel(
        path = paste0(folder.in, "/", name.ii, "/", name.ii, "_CrossDatasetTFPrioritization.xlsx"),
        sheet = "CrossDataset.Prioritization"
      )
    tfs.chosen.up <- priority.df.ii %>%
      dplyr::slice_max(order_by = mean.Scaled.Value, n = n.tfs) %>% dplyr::pull(name.universal)
    tfs.chosen.down <- priority.df.ii %>%
      dplyr::slice_min(order_by = mean.Scaled.Value, n = n.tfs) %>% dplyr::pull(name.universal)
    extreme.tfs.ii <- data.frame(
      name.universal = c(tfs.chosen.up, tfs.chosen.down),
      value = c(rep(1, times = length(tfs.chosen.up)), rep(-1, times = length(tfs.chosen.down))),
      program = name.ii
    )
    extreme.tfs <- rbind(extreme.tfs, extreme.tfs.ii)
  }

  links <- extreme.tfs %>%
    dplyr::rename(source = name.universal, target = program) %>%
    dplyr::relocate(value, .after = last_col()) %>%
    dplyr::mutate(value = factor(as.character(value),
                                 levels = c("-1", "1"),
                                 labels = c("Predict\nRepress", "Predict\nActivate")))

  nodes <- data.frame(
    node.name = unique(c(links$source, links$target))
  ) %>%
    dplyr::mutate(node.type = if_else(node.name %in% unique(links$target),
                                      "Module",
                                      "TF"))
  node.degree <- links$source %>% table() %>%
    as.data.frame() %>%
    dplyr::arrange(-Freq) %>%
    dplyr::rename(source = 1, degree = 2) %>%
    dplyr::filter(degree >= degree.display.label)

  links <- dplyr::right_join(x = links,
                             y = node.degree,
                             by = c("source"))

  nodes <- dplyr::left_join(x = nodes,
                            y = node.degree,
                            by = c("node.name" = "source")) %>%
    dplyr::mutate(node.label = node.name) %>%
    dplyr::mutate(node.label = gsub(pattern = "_", replacement = " ", x = node.label)) %>%
    dplyr::mutate(node.label = stringr::str_wrap(node.label, width = 20)) %>%
    dplyr::filter(!(is.na(degree) & (node.type == "TF"))) %>%
    dplyr::mutate(degree = if_else(is.na(degree),
                                   as.integer(9),
                                   degree))

  network <- igraph::graph_from_data_frame(d=links, vertices = nodes, directed=TRUE)

  edge.colors <- c("Predict\nRepress" = "#125175aa",
                   "Predict\nActivate" = "#ad2524aa")
  node.colors = c("Module" = "#eb7424",
                  "TF" = "#ec9c22aa")

  writexl::write_xlsx(x = links,
                      path = paste0(folder.in, "/", "CrossDataset_MultiProgram_TFPrioritization.xlsx"))

  p.tf.graph <- ggraph::ggraph(network, layout = graph.layout) +
    ggraph::geom_edge_fan(
      ggplot2::aes(colour = value,
          start_cap = ggraph::circle(0, "pt"),
          end_cap = ggraph::circle(node2.degree + 10, "pt")),
      arrow = grid::arrow(
        angle = 10,
        length = unit(0.125, "inches"),
        ends = "last",
        type = "closed"
      )) +
    ggraph::geom_node_point(aes(shape = node.type, color = node.type, size = 3)) +
    ggraph::geom_node_text(aes(label = node.label), repel=FALSE, lineheight = 0.75) +
    ggraph::scale_edge_width_continuous(range = c(0.5, 1.5)) +
    ggraph::scale_edge_color_manual(
      values = edge.colors) +
    ggplot2::scale_size(guide = "none") +
    ggplot2::scale_color_manual(values = node.colors) +
    ggraph::theme_graph(base_family="sans") +
    ggplot2::theme(legend.position = "bottom", legend.justification = "center",
                   text = element_text(size = plot.font.size))

  pdf(file = paste0(folder.in, "/", "CrossDataset_MultiProgram_TFPrioritization.pdf"),
      width = plot.width, height = plot.height)
  plot(p.tf.graph)
  dev.off()
}
