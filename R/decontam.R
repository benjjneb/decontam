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
#' @param threshold (Optional). Default \code{0.1}.
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
#' If \code{detailed=TRUE} a \code{data.frame} with additional information (such as the p-value) is returned.
#'
#' @importFrom methods as
#' @importFrom methods is
#'
#' @export
#'
#' @examples
#' isContaminant(st, conc, threshold=0.2)
#'
isContaminant <- function(seqtab, conc=NULL, neg=NULL, method="frequency", threshold = 0.1, normalize=TRUE, detailed=FALSE) {
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
  isC[is.na(isC)] <- FALSE # NA pvals are not called contaminants
  # Make return value
  if(detailed) {
    rval <- data.frame(freq=apply(seqtab,2,mean), prev=apply(seqtab>0,2,sum), p.freq=p.freq, p.prev=p.prev, pval=pval, contaminant=isC)
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
    pval <- NA
  }
  return(pval)
}
#' importFrom stats chisq.test
#' importFrom stats fisher.test
#'
#' @export
#'
### @keywords internal
isContaminantPrevalence <- function(freq, neg, method="auto") {
  fisher.pval <- function(tab, alternative) {
    excess <- fisher.test(tab, alternative="greater")$p.value + fisher.test(tab, alternative="less")$p.value - 1
    pval <- fisher.test(tab, alternative=alternative)$p.value
    pval <- pval - excess/2
    pval
  }
  if(sum(freq>0)>1 && sum(neg) > 0 && sum(neg) < length(neg)) {
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
      warning("NA p-value calculated.")
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

#' Identify non-contaminant sequences. This is the inverse function of \code{\link{isContaminant}}.
#'
#' The frequency of each sequence (or OTU) in the input feature table as a function of the concentration of
#' amplified DNA in each sample is used to identify non-contaminant sequences. Note that the null hypothesis
#' here is that sequences **are** contaminants. This function is intended for use on low-biomass samples
#' in which a large proportion of the sequences are likely to be contaminants.
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
#' @param threshold (Optional). Default \code{0.5}.
#' The p-value threshold at which the null-hypothesis (a contaminant) should be rejected in favor of the
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
#' isNotContaminant(st, conc, threshold=0.05)
#'
isNotContaminant <- function(seqtab, conc=NULL, neg=NULL, method="frequency", threshold = 0.5, normalize=TRUE, detailed=FALSE) {
  df <- isContaminant(seqtab, conc=conc, neg=neg, method=method, threshold=threshold, normalize=normalize, detailed=TRUE)
  df$p.freq <- 1-df$p.freq
  df$p.prev <- 1-df$p.prev
  # Calculate overall p-value
  if(method=="frequency") { pval <- df$p.freq }
  else if(method=="prevalence") { pval <- df$p.prev }
  else { pval <- pchisq(-2*log(df$p.freq * df$p.prev), df=4, lower.tail=FALSE) }
  # Make contaminant calls ( ALLOW LENGTH 2 THRESH FOR COMBINED METHOD? )
  isNotC <- (pval < threshold)
  isNotC[is.na(isNotC)] <- FALSE # NA pvals are not called not-contaminants
  df$pval <- pval
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
