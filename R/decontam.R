#' Identify contaminant sequences.
#'
#' The frequency of each sequence (or OTU) in the input feature table as a function of the concentration of
#' amplified DNA in each sample is used to identify contaminant sequences.
#'
#' @param seqtab (Required). Integer matrix.
#' A feature table recording the observed abundances of each sequence (or OTU) in each sample.
#' Rows should correspond to samples, and columns to sequences (or OTUs).
#'
#' @param conc (Optional). \code{numeric}.
#' A quantitative measure of the concentration of amplified DNA in each sample prior to sequencing.
#' All values must be greater than zero. Zero is assumed to represent the complete absence of DNA.
#' REQUIRED if performing frequency-based testing.
#'
#' @param neg (Optional). \code{logical}
#' The negative control samples. Extraction controls give the best results.
#' REQUIRED if performing prevalence-based testing.
#'
#' @param method(Optional). Default "frequency".
#' The method used to test for contaminants.
#' frequency: Contaminants are identified by increased frequency in lower biomass samples.
#' prevalence: Contaminants are identified by increased prevalence in negative controls.
#' combined: Both frequency and prevalence are used to identify contaminants.
#'
#' @param threshold (Optional). Default \code{1e-2}.
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
#' isContaminant(st, conc, threshold=0.05)
#'
isContaminant <- function(seqtab, conc=NULL, neg=NULL, method="frequency", threshold = 1e-2, normalize=TRUE, detailed=FALSE) {
  # Validate input
  if(!(is(seqtab, "matrix") && is.numeric(seqtab))) stop("seqtab must be a numeric matrix.")
  if(normalize) seqtab <- sweep(seqtab, 1, rowSums(seqtab), "/")
  if(!method %in% c("frequency", "prevalence", "combined")) stop("Valid method arguments: frequency, prevalence, combined")
  do.freq <- FALSE; do.prev <- FALSE; p.freq <- NA; p.prev <- NA
  if(method %in% c("frequency", "combined")) do.freq <- TRUE
  if(method %in% c("prevalence", "combined")) do.prev <- TRUE
  # Calculate frequency p-value
  if(do.freq) {
    if(!(is.numeric(conc) && all(conc>0))) stop("conc must be a positive numeric vector.")
    if(nrow(seqtab) != length(conc)) stop("The length of conc must match the number of samples (the rows of seqtab).")
    p.freq <- apply(seqtab, 2, isContaminantFrequency, conc=conc)
  }
  # Calculate prevalence p-value
  if(do.prev) {
    p.prev <- apply(seqtab, 2, isContaminantPrevalence, neg=neg)
  }
  # Calculate overall p-value
  if(method=="frequency") { pval <- p.freq }
  else if(method=="prevalence") { pval <- p.prev }
  else { pval <- pchisq(-2*log(p.freq * p.prev), df=4, lower.tail=FALSE) }
  # Make contaminant calls ( ALLOW LENGTH 2 THRESH FOR COMBINED METHOD? )
  isC <- (pval < threshold)
  # Make return value
  if(detailed) {
    rval <- data.frame(freq=apply(seqtab,2,mean), prev=apply(seqtab>0,2,sum),p.freq=p.freq, p.prev=p.prev, pval=pval, contaminant=isC)
  } else {
    rval <- isC
  }
  return(rval)
}
#' importFrom stats lm
#'
#' @export
#'
### @keywords internal
isContaminantFrequency <- function(freq, conc) {
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
    pval <- 1
  }
  return(pval)
}
#' importFrom stats chisq.test
#'
#' @export
#'
### @keywords internal
isContaminantPrevalence <- function(freq, neg) {
  tab <- table(factor(freq>0, levels=c(FALSE, TRUE)), factor(neg, levels=c(FALSE, TRUE)))
  chi1 <- suppressWarnings(chisq.test(tab, simulate.p.value=TRUE, B=10000))
  pval <- chi1$p.value
  if(tab[2,2]/(tab[1,2]+tab[2,2]) < tab[2,1]/(tab[1,1]+tab[2,1])) { # more freq in positives
    pval <- 1
  } else { # In the expected direction
    pval <- pval/2
  }
  if(is.na(pval)) pval <- 1
  return(pval)
}

