#
#' plot_abundance
#'
#' plots the abundance of selected taxa in each sample vs. each sample's DNA concentration after PCR amplification
#'
#' @param seqtab (Required). \code{Integer matrix} or \code{phyloseq} object.
#' A feature table recording the observed abundances of each sequence variant (or OTU) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#' If a phyloseq object is provided, the otu-table component will be extracted.
#'
#' @param taxa_to_plot (Required). \code{character}.
#' The names of the taxa to include in this plot.
#'
#' @param conc (Required). \code{numeric}. Required if performing frequency-based testing.
#' A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
#' All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
#' If \code{seqtab} was prodivded as a phyloseq object, the name of the appropriate sample-variable in that
#' phyloseq object can be provided.
#'
#' @param normalize (Optional). Default TRUE.
#' If TRUE, the input \code{seqtab} is normalized so that each row sums to 1 (converted to frequency).
#' If FALSE, no normalization is performed (the data should already be frequencies or counts from equal-depth samples).
#'
#' @param showModels (Optional). Default TRUE.
#' If TRUE, the contaminant (red) and non-contaminant (black) models are shown in the plot.
#'
#' @param log (Optional). Default TRUE.
#' If TRUE, the axes are log10-scaled.
#'
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @import ggplot2
#'
#' @export
#'
#' @examples
#'
plot_abundance <- function(seqtab, taxa_to_plot, conc, normalize=TRUE, showModels=TRUE, log=TRUE){
  # Validate input
  if(is(seqtab, "phyloseq")) {
    ps <- seqtab
    seqtab <- as(ps@otu_table, "matrix")
    if(ps@otu_table@taxa_are_rows) { seqtab <- t(seqtab) }
    if(is.character(conc) && length(conc)==1) { conc <- getFromPS(ps, conc) }
  } else {
    ps <- NULL # No phyloseq object
  }
  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  if(!(is.numeric(conc) && all(conc>0))) stop("conc must be positive numeric.")
  if(is.character(taxa_to_plot)) {
    seqtab <- seqtab[,colnames(seqtab) %in% taxa_to_plot,drop=FALSE]
  } else {
    stop("taxa must be a vector of taxa names.")
  }
  if(ncol(seqtab) == 0) stop("None of the provided taxa were present.")
  # Prepare plotting data.frame
  if(is.null(ps)) {
    plotdf <- cbind(data.frame(seqtab), DNA_conc=conc)
  } else {
    plotdf <- cbind(data.frame(seqtab), data.frame(sample_data(ps)), DNA_conc=conc)
  }
  plot_melt <- melt(plotdf, measure.vars=1:ncol(seqtab), variable.name="taxa_to_plot", value.name="taxon_abundance")

  taxon_levels <- taxa_to_plot
  plot_melt$taxa_to_plot <- factor(plot_melt$taxa_to_plot, levels = taxon_levels)

  p1 <- ggplot(data=plot_melt, aes(DNA_conc, taxon_abundance)) + xlab("DNA Concentration") + ylab("Relative Abundance")
  if(log) p1 <- p1 + scale_x_log10() + scale_y_log10()
  p1 + geom_point()
}

##Examples
#plot_abundance(fungal,"denovo7","qPCR_copies",taxa_are_rows=TRUE)

#plot_abundance(MUC,"Seq1","quant_reading", taxa_are_rows=FALSE)
#plot_abundance(MUC,c("Seq152","Seq1"),"quant_reading",taxa_are_rows=FALSE)
#plot_abundance(MUC,c("Seq152","Seq1"),"quant_reading",taxa_are_rows=FALSE) + facet_wrap(~taxa_to_plot)

#plot_abundance(cdiff, c("42372","9710"), conc="Fluorescence")
#p <- plot_abundance(cdiff, c("42372","9710"), conc="Fluorescence")
#p + facet_wrap(~OTU)
#p + geom_point(aes(color=Treatment)) + facet_wrap(~OTU)
