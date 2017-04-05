#' Identify contaminant sequences.
#'
#' The frequency of each sequence (or OTU) in the input feature table as a function of the concentration of
#' amplified DNA in each sample is used to identify contaminant sequences.
#'
#' @param seqtab (Required). Integer matrix.
#' A feature table recording the observed abundances of each sequence (or OTU) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#'
#' @param conc (Required). \code{numeric}.
#' A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
#' All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
#'
#' @param threshold (Optional). Default \code{1e-3}.
#' The p-value threshold at which the null-hypothesis (not a contaminant) should be rejected in favor of the
#' alternate hypothesis (contaminant).
#'
#' @param normalize (Optional). Default TRUE.
#' If TRUE, the input \code{seqtab} is normalized so that each row sums to 1 (converted to frequency).
#' If FALSE, no normalization is performed (the data should already be frequencies or counts from equal-depth samples).
#'
#' @param detailed (Optional). Default FALSE.
#' If TRUE, the return value is a \code{data.frame} containing diagnostic information on the contaminant decision.
#' If FALSE, the return value is a \code{logical} vector containing the contaminant decisions.
#'
#' @return
#' If \code{detailed=FALSE} a \code{logical} vector is returned, with TRUE indicating contaminants.
#' If \code{detailed=TRUE} a \code{data.frame} is returned instead.
#'
#' @importFrom methods as
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' filterAndTrim(testFastqs, filtFastqs, truncQ=2, truncLen=200, rm.phix=TRUE)
#'
isContaminant <- function(seqtab, conc, threshold = 1e-3, normalize=TRUE, detailed=FALSE) {
  # Validate input
  if(!(is(seqtab, "matrix") && is.numeric(seqtab))) stop("seqtab must be a numeric matrix.")
  if(!(is.numeric(conc) && all(conc>0))) stop("conc must be a positive numeric vector.")
  if(nrow(seqtab) != length(conc)) stop("The length of conc must match the number of samples (the rows of seqtab).")
  # Test each input feature (column)
  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  out <- t(apply(seqtab, 2, isContaminantSingle, conc=conc))
  # Make return value
  if(detailed) {
    rval <- data.frame(out)
    colnames(rval) <- c("SS.0", "SS.1", "dof", "pval")
  } else {
    rval <- (out[,4] < threshold)
  }
  return(rval)
}
#' importFrom stats lm
#' @keywords internal
isContaminantSingle <- function(freq, conc) {
  df <- data.frame(logc=log(conc), logf=log(freq))
  df <- df[freq>0,]
  if(sum(freq>0)>1) {
    lm1 <- lm(logf~offset(-1*logc), data=df)
    SS1 <- sum(lm1$residuals^2)
    lm0 <- lm(logf~1, data=df)
    SS0 <- sum(lm0$residuals^2)
    dof <- sum(freq>0)-1
    pval <- pf(SS1/SS0,dof,dof)
  } else {
    SS1 <- NA
    SS0 <- NA
    dof <- NA
    pval <- 1
  }
  return(c(SS0, SS1, dof, pval))
}

