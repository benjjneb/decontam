#' Plot frequencies as a function of input DNA concentration
#'
#' Plots the frequencies of selected sequence features vs. each sample's DNA concentration.
#'
#' @param seqtab (Required). \code{Integer matrix} or \code{phyloseq} object.
#' A feature table recording the observed abundances of each sequence feature (e.g. OTUs or ASVs or
#' or genus or ortholog or...) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#' If a phyloseq object is provided, the otu-table component will be extracted.
#'
#' @param taxa (Required). \code{character}.
#' The names of the sequence features to include in this plot. Should match \code{colnames(setab)}
#' if a \code{matrix} was provided, or \code{taxa_names(seqtab)} if a \code{phyloseq} object was provided.
#'
#' @param conc (Required). \code{numeric}.
#' A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
#' All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
#' If \code{seqtab} was provided as a \code{phyloseq} object, the name of the sample variable in the
#' \code{phyloseq} object can be provided.
#'
#' @param neg (Optional). \code{logical}. Default NULL.
#' TRUE if sample is a negative control, and FALSE if not.
#' If \code{seqtab} was provided as a phyloseq object, the name of the appropriate sample-variable in that
#' phyloseq object can be provided. NULL indicates no samples should be condired negative controls.
#'
#' @param normalize (Optional). Default TRUE.
#' If TRUE, the input \code{seqtab} is normalized so that each row sums to 1 (converted to frequency).
#' If FALSE, no normalization is performed (the data should already be frequencies or counts from
#' equal-depth samples).
#'
#' @param showModels (Optional). Default TRUE.
#' If TRUE, the contaminant (red, dashed line) and non-contaminant (black, solid line) models are
#' shown in the plot.
#'
#' @param log (Optional). Default TRUE.
#' If TRUE, the axes are log10-scaled.
#'
#' @param facet (Optional). Default TRUE.
#' If TRUE, multiple sequence features will be plotted in separate facets.
#'
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom reshape2 melt
#' @importFrom stats lm
#' @import ggplot2
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # MUC is a phyloseq object
#'   plot_frequency(MUC,"Seq1",conc="quant_reading")
#'   plot_frequency(MUC,c("Seq1", "Seq10", "Seq33"),conc=sample_data(MUC)$quant_reading)
#'   plot_frequency(MUC,"Seq1",conc="quant_reading", normalize=FALSE, log=FALSE)
#' }
plot_frequency <- function(seqtab, taxa, conc, neg=NULL, normalize=TRUE, showModels=TRUE, log=TRUE, facet=TRUE){
  # Validate input
  if(is(seqtab, "phyloseq")) {
    ps <- seqtab
    seqtab <- as(ps@otu_table, "matrix")
    if(ps@otu_table@taxa_are_rows) { seqtab <- t(seqtab) }
    if(is.character(conc) && length(conc)==1) { conc <- getFromPS(ps, conc) }
    if(is.character(neg) && length(neg)==1) { neg <- getFromPS(ps, neg) }
  } else {
    ps <- NULL # No phyloseq object
  }
  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  if(!(is.numeric(conc) && all(conc>0))) stop("conc must be positive numeric.")
  if(is.null(neg)) neg <- rep(FALSE, length(conc)) # Don't ignore any samples
  if(is.character(taxa)) {
    seqtab <- seqtab[,colnames(seqtab) %in% taxa,drop=FALSE]
  } else {
    stop("taxa must be a vector of taxa names.")
  }
  ntax.plot <- ncol(seqtab)
  if(ntax.plot == 0) stop("None of the provided taxa were present in seqtab.")
  # Prepare plotting data.frame
  if(is.null(ps)) {
    plotdf <- cbind(data.frame(seqtab, check.names=FALSE), DNA_conc=conc, Type=ifelse(neg, "Negative", "Sample"))
  } else {
    plotdf <- cbind(data.frame(seqtab, check.names=FALSE), data.frame(sample_data(ps), check.names=FALSE),
                    DNA_conc=conc, Type=ifelse(neg, "Negative", "Sample"))
  }
  plot_melt <- melt(plotdf, measure.vars=1:ntax.plot, variable.name="taxa", value.name="taxon_abundance")

  taxon_levels <- taxa
  plot_melt$taxa <- factor(plot_melt$taxa, levels = taxon_levels)

  if(showModels) {
    #    moddf <- data.frame(DNA_conc=seq(min(plotdf$DNA_conc), max(plotdf$DNA_conc), length.out=200))
    mod_melts <- split(plot_melt, plot_melt$taxa)
    logc <- log(seq(min(plotdf$DNA_conc), max(plotdf$DNA_conc), length.out=1000))
    for(tax in names(mod_melts)) {
      # Code copied in from isContaminantFrequency
      # Should really consider functionalizing this to avoid issues with the duplicated code
      newdata <- data.frame(logc=logc, taxa=tax, DNA_conc=exp(logc))
      freq <- mod_melts[[tax]]$taxon_abundance
      conc <- mod_melts[[tax]]$DNA_conc
      df <- data.frame(logc=log(conc), logf=log(freq))
      df <- df[!neg | is.na(neg),]
      df <- df[freq>0,]
      if(sum(freq>0)>1) {
        lm1 <- lm(logf~offset(-1*logc), data=df)
        lm0 <- lm(logf~1, data=df)
        newdata$contam <- exp(predict(lm1, newdata=newdata))
        newdata$non.contam <- exp(predict(lm0, newdata=newdata))
      } else {
        newdata$contam <- NA
        newdata$non.contam <- NA
      }
      mod_melts[[tax]] <- newdata
    }
    mod_melt <- do.call(rbind, mod_melts)
  }

  p1 <- ggplot(data=plot_melt, aes(DNA_conc, taxon_abundance)) + xlab("DNA Concentration")
  p1 <- p1 + ylab(ifelse(normalize, "Frequency", "Relative Abundance"))
  if(log) p1 <- p1 + scale_x_log10()
  if(log) p1 <- p1 + scale_y_log10(limits=c(NA, ifelse(normalize || all(plot_melt$DNA_conc <= 1), 1, NA)))
  if(nlevels(factor(neg))>1) p1 <- p1 + aes(color=Type)
  if(facet && ntax.plot > 1) p1 <- p1 + facet_wrap(~taxa)
  if(showModels) p1 <- p1 + geom_line(data=mod_melt, aes(y=contam), color="red", linetype="solid")
  if(showModels) p1 <- p1 + geom_line(data=mod_melt, aes(y=non.contam), color="black", linetype="dashed")
  p1 + geom_point()
}
