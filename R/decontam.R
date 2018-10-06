#' Identify contaminant sequences.
#'
#' The frequency of each sequence (or OTU) in the input feature table as a function of the concentration of
#' amplified DNA in each sample is used to identify contaminant sequences.
#'
#' @param seqtab (Required). \code{Integer matrix} or \code{phyloseq} object.
#' A feature table recording the observed abundances of each sequence variant (or OTU) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#' If a phyloseq object is provided, the otu-table component will be extracted.
#'
#' @param conc (Optional). \code{numeric}. Required if performing frequency-based testing.
#' A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
#' All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
#' If \code{seqtab} was prodivded as a phyloseq object, the name of the appropriate sample-variable in that
#' phyloseq object can be provided.
#'
#' @param neg (Optional). \code{logical}. Required if performing prevalence-based testing.
#' TRUE if sample is a negative control, and FALSE if not (NA entries are not included in the testing).
#' Extraction controls give the best results.
#' If \code{seqtab} was provided as a phyloseq object, the name of the appropriate sample-variable in that
#' phyloseq object can be provided.
#'
#' @param method (Optional). \code{character}. The method used to test for contaminants.
#' \describe{
#'   \item{auto}{(Default). frequency, prevalence or combined will be automatically selected based on whether
#'               just \code{conc}, just \code{neg}, or both were provided.}
#'   \item{frequency}{Contaminants are identified by frequency that varies inversely with sample DNA concentration.}
#'   \item{prevalence}{Contaminants are identified by increased prevalence in negative controls.}
#'   \item{combined}{The frequency and prevalence probabilities are combined with Fisher's method and used to identify contaminants.}
#'   \item{minimum}{The minimum of the frequency and prevalence probabilities is used to identify contaminants.}
#'   \item{either}{Contaminants are called if identified by either the frequency or prevalance methods.}
#'   \item{both}{Contaminants are called if identified by both the frequency and prevalance methods.}
#' }
#'
#' @param batch (Optional). \code{factor}, or any type coercible to a \code{factor}. Default NULL.
#' If provided, should be a vector of length equal to the number of input samples which specifies which batch
#' each sample belongs to (eg. sequencing run). Contaminants identification will be performed independently
#' within each batch.
#' If \code{seqtab} was provided as a phyloseq object, the name of the appropriate sample-variable in that
#' phyloseq object can be provided.
#'
#' @param batch.combine (Optional). Default "minimum".
#' For each input sequence variant (or OTU) the probabilities calculated in each batch are combined into a
#' single probability that is compared to `code{threshold}` to classify contaminants.
#' Valid values: "minimum", "product", "fisher".
#'
#' @param threshold (Optional). Default \code{0.1}.
#' The probability threshold below which (strictly less than) the null-hypothesis (not a contaminant) should be rejected in favor of the
#' alternate hypothesis (contaminant). A length-two vector can be provided when using the \code{either} or \code{both} methods:
#' the first value is the threshold for the frequency test and the second for the prevalence test.
#'
#' @param normalize (Optional). Default TRUE.
#' If TRUE, the input \code{seqtab} is normalized so that each row sums to 1 (converted to frequency).
#' If FALSE, no normalization is performed (the data should already be frequencies or counts from equal-depth samples).
#'
#' @param detailed (Optional). Default TRUE.
#' If TRUE, the return value is a \code{data.frame} containing diagnostic information on the contaminant decision.
#' If FALSE, the return value is a \code{logical} vector containing the binary contaminant classifications.
#'
#' @return
#' If \code{detailed=TRUE} a \code{data.frame} with classification information.
#' If \code{detailed=FALSE} a \code{logical} vector is returned, with TRUE indicating contaminants.
#'
#' @importFrom methods as
#' @importFrom methods is
#' @importFrom stats pchisq
#'
#' @export
#'
#' @examples
#' st <- readRDS(system.file("extdata", "st.rds", package="decontam"))
#' # conc should be positive and non-zero
#' conc <- c(6413, 3581.0, 5375, 4107, 4291, 4260, 4171, 2765, 33, 48)
#' neg <- c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
#' # Use frequency or frequency and prevalence to identify contaminants
#' isContaminant(st, conc=conc, method="frequency", threshold=0.2)
#' isContaminant(st, conc=conc, neg=neg, method="both", threshold=c(0.1,0.5))
#'
isContaminant <- function(seqtab, conc=NULL, neg=NULL,
                          method=c("auto", "frequency", "prevalence", "combined", "minimum", "either", "both"),
                          batch=NULL, batch.combine=c("minimum", "product", "fisher"),
                          threshold = 0.1, normalize=TRUE, detailed=TRUE) {
  # Validate input
  if(is(seqtab, "phyloseq")) {
    ps <- seqtab
    seqtab <- as(ps@otu_table, "matrix")
    if(ps@otu_table@taxa_are_rows) { seqtab <- t(seqtab) }
    if(is.character(conc) && length(conc)==1) { conc <- getFromPS(ps, conc) }
    if(is.character(neg) && length(neg)==1) { neg <- getFromPS(ps, neg) }
    if(is.character(batch) && length(batch)==1) { batch <- getFromPS(ps, batch) }
  }
  if(!(is(seqtab, "matrix") && is.numeric(seqtab))) stop("seqtab must be a numeric matrix.")
  if(any(rowSums(seqtab) == 0)) { # Catch and remove zero-count samples
    zero.count <- rowSums(seqtab) == 0
    seqtab <- seqtab[!zero.count,]
    if(!is.null(conc)) conc <- conc[!zero.count]
    if(!is.null(neg)) neg <- neg[!zero.count]
    if(!is.null(batch)) batch <- batch[!zero.count]
    warning("Removed ", sum(zero.count), " samples with zero total counts (or frequency).")
  }
  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  method <- match.arg(method)
  if(method == "auto") {
    if(!missing(conc) && missing(neg)) method <- "frequency"
    else if(missing(conc) && !missing(neg)) method <- "prevalence"
    else method <- "combined"
  }
  do.freq <- FALSE; do.prev <- FALSE; p.freq <- NA; p.prev <- NA
  if(method %in% c("frequency", "minimum", "combined", "minimum", "either", "both")) do.freq <- TRUE
  if(method %in% c("prevalence", "combined", "minimum", "either", "both")) do.prev <- TRUE
  if(do.prev) {
    if(missing(neg)) stop("neg must be provided to perform prevalence-based contaminant identification.")
  }
  if(do.freq) {
    if(missing(conc)) stop("conc must be provided to perform frequency-based contaminant identification.")
    if(!(is.numeric(conc) && all(conc>0))) stop("conc must be positive numeric.")
    if(nrow(seqtab) != length(conc)) stop("The length of conc must match the number of samples (the rows of seqtab).")
    if(missing(neg)) neg <- rep(FALSE, length(conc)) # Don't ignore any samples
  }
  if(is.numeric(threshold) && all(threshold >= 0) && all(threshold <= 1)) {
    if(method %in% c("either", "both")) {
      if(length(threshold) == 1) {
        message("Using same threshold value for the frequency and prevalence contaminant identification.")
        threshold <- c(threshold, threshold)
      }
    } else if(length(threshold) != 1) {
      stop("threshold should be a single value.")
    }
  } else {
    stop("threshold must be a numeric value from 0 to 1 (inclusive).")
  }
  if(missing(batch) || is.null(batch)) {
    batch <- factor(rep(1, nrow(seqtab)))
  }
  if(nrow(seqtab) != length(batch)) stop("The length of batch must match the number of samples (the rows of seqtab).")
  batch.combine <- match.arg(batch.combine)
  batch <- factor(batch)
  # Loop over batches
  p.freqs <- matrix(NA, nrow=nlevels(batch), ncol=ncol(seqtab))
  rownames(p.freqs) <- levels(batch)
  p.prevs <- matrix(NA, nrow=nlevels(batch), ncol=ncol(seqtab))
  rownames(p.prevs) <- levels(batch)
  for(bat in levels(batch)) {
    # Calculate frequency p-value
    if(do.freq) {
      p.freqs[bat,] <- apply(seqtab[batch==bat & !neg,], 2, isContaminantFrequency, conc=conc[batch==bat & !neg])
    }
    # Calculate prevalence p-value
    if(do.prev) {
      p.prevs[bat,] <- apply(seqtab[batch==bat,], 2, isContaminantPrevalence, neg=neg[batch==bat])
    }
  }
  # Combine batch p-values
  if(batch.combine == "minimum") {
    if(do.freq) {
      suppressWarnings(p.freq <- apply(p.freqs, 2, min, na.rm=TRUE))
      p.freq[is.infinite(p.freq)] <- NA # If NA in all batches, min sets to infinite
    }
    if(do.prev) {
      suppressWarnings(p.prev <- apply(p.prevs, 2, min, na.rm=TRUE))
      p.prev[is.infinite(p.prev)] <- NA # If NA in all batches, min sets to infinite
    }
  } else if(batch.combine == "product") {
    if(do.freq) {
      suppressWarnings(p.freq <- apply(p.freqs, 2, prod, na.rm=TRUE))
    }
    if(do.prev) {
      suppressWarnings(p.prev <- apply(p.prevs, 2, prod, na.rm=TRUE))
    }
  } else if(batch.combine == "fisher") {
    if(do.freq) {
      p.freq <- apply(p.freqs, 2, fish.combine, na.replace=0.5)
    }
    if(do.prev) {
      p.prev <- apply(p.prevs, 2, fish.combine, na.replace=0.5)
    }
  } else {
    stop("Invalid batch.combine value.")
  }
  # Calculate overall p-value
  if(method=="frequency") { pval <- p.freq }
  else if(method=="prevalence") { pval <- p.prev }
  else if(method=="minimum") { pval <- pmin(p.freq, p.prev) }
  else if(method=="combined") { pval <- pchisq(-2*log(p.freq * p.prev), df=4, lower.tail=FALSE) }
  else if(method %in% c("either", "both")) { pval <- rep(NA, length(p.freq)) }
  else { stop("Invalid method specified.") }

  if(method=="either") { # Two tests
    isC <- (p.freq < threshold[[1]]) | (p.prev < threshold[[2]])
  } else if(method =="both") {
    isC <- (p.freq < threshold[[1]]) & (p.prev < threshold[[2]])
  } else { # One test
    isC <- (pval < threshold)
  }
  isC[is.na(isC)] <- FALSE # NA pvals are not called contaminants
  # Make return value
  if(detailed) {
    rval <- data.frame(freq=apply(seqtab,2,mean), prev=apply(seqtab>0,2,sum), p.freq=p.freq, p.prev=p.prev, p=pval, contaminant=isC)
  } else {
    rval <- isC
  }
  return(rval)
}
#' @importFrom stats lm
#' @importFrom stats pf
#'
#' @keywords internal
isContaminantFrequency <- function(freq, conc) {
  df <- data.frame(logc=log(conc), logf=log(freq))
  df <- df[!is.na(freq) & freq>0,]
  if(nrow(df)>1) {
    lm1 <- lm(logf~offset(-1*logc), data=df)
    SS1 <- sum(lm1$residuals^2)
    lm0 <- lm(logf~1, data=df)
    SS0 <- sum(lm0$residuals^2)
    dof <- sum(freq>0)-1
    pval <- pf(SS1/SS0,dof,dof)
  } else {
    pval <- NA
  }
  return(pval)
}
#' @importFrom stats fisher.test
#' @importFrom stats prop.test
#'
#' @keywords internal
isContaminantPrevalence <- function(freq, neg, method="auto") {
  fisher.pval <- function(tab, alternative) {
    excess <- fisher.test(tab, alternative="greater")$p.value + fisher.test(tab, alternative="less")$p.value - 1
    pval <- fisher.test(tab, alternative=alternative)$p.value
    pval <- pval - excess/2
    pval
  }
  if(sum(freq>0)>1 && sum(neg,na.rm=TRUE) > 0 && sum(neg,na.rm=TRUE) < sum(!is.na(neg))) {
    tab <- table(factor(neg, levels=c(TRUE, FALSE)), factor(freq>0, levels=c(TRUE, FALSE)))
    # First entry (1,1) is the neg prevalence, so alternative is "greater"
    if((tab[1,2] + tab[2,2]) == 0) { # Present in all samples
      pval <- 0.5
    } else if(method == "fisher") {
      pval <- fisher.pval(tab, alternative="greater")
    } else if(method == "chisq") {
      pval <- prop.test(tab, alternative="greater")$p.value
    } else {
      pval <- tryCatch(prop.test(tab, alternative="greater")$p.value, warning=function(w) fisher.pval(tab, alternative="greater"))
    }
    if(is.na(pval)) {
      warning("NA probability calculated.")
    }
  } else {
    pval <- NA
  }
  return(pval)
}
# fisher.test(matrix(c(1,10,40,40), nrow=2), alternative="greater")
# contingency table, test is whether the first entry is less than expected under fixed marginals
# so, is there a lower fraction of the first column in row 1 than row 2
# prop.test(matrix(c(1,10,40,40), nrow=2), alternative="greater")
# Same test but using chisq approx, which fails at low numbers
# Warns for low numbers (conditionally use prop.test based on that?)
# tab <- table(factor(df$neg, levels=c(TRUE, FALSE)), factor(df$present, levels=c(TRUE, FALSE)))
# tryCatch(prop.test(tab, alternative="less"), warning=function(w) fisher.test(tab, alternative="less"))

#' Identify non-contaminant sequences.
#'
#' The prevalence of each sequence (or OTU) in the input feature table across samples and negative controls
#' is used to identify non-contaminant sequences. Note that the null hypothesis
#' here is that sequences **are** contaminants. This function is intended for use on low-biomass samples
#' in which a large proportion of the sequences are likely to be contaminants.
#'
#' @param seqtab (Required). Integer matrix.
#' A feature table recording the observed abundances of each sequence (or OTU) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#'
#' @param neg (Required). \code{logical}
#' The negative control samples. Extraction controls give the best results.
#'
#' @param method (Optional). Default "prevalence".
#' The method used to test for contaminants. Currently the only method supported is prevalence.
#' prevalence: Contaminants are identified by increased prevalence in negative controls.
#'
#' @param threshold (Optional). Default \code{0.5}.
#' The probability threshold below which (strictly less than) the null-hypothesis (a contaminant) should be rejected in favor of the
#' alternate hypothesis (not a contaminant).
#'
#' @param normalize (Optional). Default TRUE.
#' If TRUE, the input \code{seqtab} is normalized so that each row sums to 1 (converted to frequency).
#' If FALSE, no normalization is performed (the data should already be frequencies or counts from equal-depth samples).
#'
#' @param detailed (Optional). Default FALSE.
#' If TRUE, the return value is a \code{data.frame} containing diagnostic information on the non-contaminant decision.
#' If FALSE, the return value is a \code{logical} vector containing the non-contaminant decisions.
#'
#' @return
#' If \code{detailed=FALSE} a \code{logical} vector is returned, with TRUE indicating non-contaminants.
#' If \code{detailed=TRUE} a \code{data.frame} is returned instead.
#'
#' @export
#'
#' @examples
#' st <- readRDS(system.file("extdata", "st.rds", package="decontam"))
#' samdf <- readRDS(system.file("extdata", "samdf.rds", package="decontam"))
#' isNotContaminant(st, samdf$quant_reading, threshold=0.05)
#'
isNotContaminant <- function(seqtab, neg=NULL, method="prevalence", threshold = 0.5, normalize=TRUE, detailed=FALSE) {
  if(!method %in% c("prevalence")) stop("isNotContaminant only supports the following methods: prevalence")
  df <- isContaminant(seqtab, conc=NULL, neg=neg, method=method, threshold=threshold, normalize=normalize, detailed=TRUE)
  df$p.freq <- 1-df$p.freq
  df$p.prev <- 1-df$p.prev
  # Calculate overall p-value
  if(method=="prevalence") { pval <- df$p.prev }
  # Make contaminant calls
  isNotC <- (pval < threshold)
  isNotC[is.na(isNotC)] <- FALSE # NA pvals are not called not-contaminants
  df$p <- pval
  df$contaminant <- NULL
  df$not.contaminant <- isNotC
  # Make return value
  if(detailed) {
    rval <- df
  } else {
    rval <- isNotC
  }
  return(rval)
}

#' @keywords internal
list_along <- function(nm) {
  if(!is.character(nm)) stop("list_along requires character input.")
  rval <- vector("list", length(nm))
  names(rval) <- nm
}

#' @keywords internal
#' @importFrom stats pchisq
fish.combine <- function(vec, na.replace=NA) {
  vec[is.na(vec)] <- na.replace
  vec <- vec[!is.na(vec)]
  if(any(vec<0 | vec>1)) stop("fish.combine expects values between 0 and 1.")
  p <- prod(vec)
  pchisq(-2*log(p), df=2*length(vec), lower.tail=FALSE)
}

#' @keywords internal
getFromPS <- function(ps, nm) {
  i <- match(nm, ps@sam_data@names)
  if(is.na(i)) stop(paste(nm, "is not a valid sample-variable in the provided phyloseq object."))
  ps@sam_data@.Data[[i]]
}
