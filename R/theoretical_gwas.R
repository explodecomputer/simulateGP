#' Calculate expected MSE
#'
#' @param beta array of effect sizes
#' @param af array of allele frequencies
#' @param vy variance of y
#'
#' @export
#' @return Numeric
expected_mse <- function(beta, af, vy)
{
	vy - beta^2 * 2 * af * (1-af)
}

#' Calculate expected SSX
#'
#' @param af array of allele frequencies
#' @param n sample size
#'
#' @export
#' @return Numeric
expected_ssx <- function(af, n)
{
	(2 * af * (1-af)) * (n - 1)
}

#' Expected se given beta, af, n and vy
#'
#' se = sqrt(sigma_e^2 / ss(x))
#'
#' @param beta array of effect sizes
#' @param af array of allele frequencies
#' @param n sample size
#' @param vy variance of y
#'
#' @export
#' @return array of standard errors
expected_se <- function(beta, af, n, vy)
{
	sqrt(expected_mse(beta, af, vy) / expected_ssx(af, n))
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
	stats::rnorm(length(beta), beta, se)
}

#' Create a GWAS summary dataset
#'
#'
#' @param beta Array of true beta values
#' @param af Array of effect allele frequency values
#' @param nid sample size
#' @param vy Variance of trait
#' @param minmaf minimum allowed maf. default=0.01 to prevent instability
#'
#' @export
#' @return list of data frames
generate_gwas_ss <- function(beta, af, nid, vy=1, minmaf=0.01)
{
	stopifnot(length(beta) == length(af))
	stopifnot(all(af > 0 & af < 1))
	af <- pmax(minmaf, af)
	af <- pmin(1-minmaf, 1-af)
	dat <- dplyr::tibble(
		snp = 1:length(beta),
		b = beta,
		maf = af,
		se = expected_se(b, af, nid, vy),
		bhat = sample_beta(b, se),
		fval = (bhat/se)^2,
		n = nid,
		pval = pf(fval, df1=1, df2=nid-1, lower.tail=FALSE) * 2
	)
	return(dat)
}


#' Generate SNP effects given MAF, h2 and selection
#'
#' @param af Vector of effect allele frequencies, one for each SNP
#' @param h2 Variance explained by all SNPs
#' @param S Selection coefficient on trait. Default = 0
#'
#' @export
#' @return data frame
generate_gwas_params <- function(af, h2, S=0)
{
	nsnp <- length(af)
	if(h2 == 0)
	{
		return(dplyr::tibble(beta=0, af=af))
	}
	beta <- stats::rnorm(nsnp, mean=0, sd = sqrt((af * 2 * (1-af))^S))
	vg <- sum(af * 2 * (1-af) * beta^2)
	ve <- (vg - h2 * vg) / h2
	vy <- vg + ve
	beta <- beta / sqrt(vy)
	return(dplyr::tibble(beta=beta, af=af))
}
