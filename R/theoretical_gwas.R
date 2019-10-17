
#' Expected se given beta, maf, n and vy
#'
#' se = sqrt(sigma_e^2 / ss(x))
#'
#' @param beta array of effect sizes
#' @param maf array of allele frequencies
#' @param n sample size
#' @param vy variance of y
#'
#' @export
#' @return array of standard errors
expected_se <- function(beta, maf, n, vy)
{
	sqrt((vy - beta^2 * 2 * maf * (1-maf)) / ((2 * maf * (1-maf)) * n))
}

#' Get the expected se for a gwas given n, h2, beta maf
#'
#'
#' @param n sample size
#' @param h2 heritability
#' @param beta array of effect sizes
#' @param maf array of allele frequencies
#'
#' @export
#' @return array of standard errors
gwas_se <- function(n, h2, beta, maf)
{
	nsnp <- length(beta)
	stopifnot(length(maf) == nsnp)
	vg <- sum(beta^2 * 2 * maf * (1-maf))
	ve <- (vg - h2 * vg) / h2
	vy <- vg + ve
	se <- expected_se(beta, maf, n, vy)
	return(se)
}

#' Sample beta values given standard errors
#'
#' @param beta array of beta values
#' @param se array of se values
#'
#' @export
#' @return array of beta hats
sample_beta <- function(beta, se)
{
	rnorm(length(beta), beta, se)
}

#' Create a theoretical GWAS dataset
#'
#' Choose nsnp, nid, h2
#' Create true effects
#' Create expected se
#' Create sampled effects
#'
#' @param beta Array of beta values
#' @param maf array of maf values
#' @param h2 h2 of trait
#' @param nid sample size
#' @param minmaf minimum allowed maf. default=0.01 to prevent instability
#'
#' @export
#' @return list of data frames
theoretical_gwas <- function(beta, maf, h2, nid, minmaf=0.01)
{
	stopifnot(length(beta) == length(maf))
	stopifnot(all(maf <= 0.5))
	maf <- pmax(minmaf, maf)
	dat <- dplyr::tibble(
		snp = 1:length(beta),
		beta = beta,
		maf = maf,
		se = gwas_se(nid, h2, beta, maf),
		betahat = sample_beta(beta, se),
		pval = pnorm(abs(betahat / se), low=FALSE) * 2
	)
	return(dat)
}


