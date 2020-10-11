#' Calculate expected MSE
#'
#' @param beta <what param does>
#' @param maf <what param does>
#' @param vy <what param does>
#'
#' @export
#' @return Numeric
expected_mse <- function(beta, maf, vy)
{
	vy - beta^2 * 2 * maf * (1-maf)
}

#' Calculate expected SSX
#'
#' @param maf <what param does>
#' @param n <what param does>
#'
#' @export
#' @return Numeric
expected_ssx <- function(maf, n)
{
	(2 * maf * (1-maf)) * (n - 1)
}

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
	sqrt(expected_mse(beta, maf, vy) / expected_ssx(maf, n))
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

#' Create a GWAS summary dataset
#'
#'
#' @param beta Array of true beta values
#' @param maf array of maf values
#' @param nid sample size
#' @param vy Variance of trait
#' @param minmaf minimum allowed maf. default=0.01 to prevent instability
#' @param reference If present then matches the LD to the 
#'
#' @export
#' @return list of data frames
generate_gwas_ss <- function(beta, maf, nid, vy=1, minmaf=0.01, reference=NULL)
{
	stopifnot(length(beta) == length(maf))
	stopifnot(all(maf > 0 & maf < 1))
	maf <- pmax(minmaf, maf)
	dat <- dplyr::tibble(
		snp = 1:length(beta),
		b = beta,
		maf = maf,
		se = expected_se(b, maf, nid, vy),
		bhat = sample_beta(b, se),
		pval = pnorm(abs(bhat / se), low=FALSE) * 2
	)
	return(dat)
}


#' Generate SNP effects given MAF, h2 and selection
#'
#' @param maf Vector of allele frequencies, one for each SNP
#' @param h2 Variance explained by all SNPs
#' @param S Selection coefficient on trait. Default = 0
#'
#' @export
#' @return data frame of maf and beta
generate_gwas_params <- function(maf, h2, S=0)
{
	nsnp <- length(maf)
	if(h2 == 0)
	{
		return(dplyr::tibble(beta=0, maf=maf))
	}
	beta <- rnorm(nsnp, mean=0, sd = sqrt((maf * 2 * (1-maf))^S))
	vg <- sum(maf * 2 * (1-maf) * beta^2)
	ve <- (vg - h2 * vg) / h2
	vy <- vg + ve
	beta <- beta / sqrt(vy)
	return(dplyr::tibble(beta=beta, maf=maf))
}
