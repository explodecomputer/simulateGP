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


#' Generate SNP effects given MAF, h2 and selection
#'
#' @param af Vector of effect allele frequencies, one for each SNP
#' @param h2 Variance explained by all SNPs
#' @param S Selection coefficient on trait. Default = 0
#' @param snp Variant identifiers. Default NULL means 1:nsnp
#'
#' @export
#' @return data frame
generate_gwas_params <- function(af, h2, S=0, snp=NULL)
{
	nsnp <- length(af)
	if(is.null(snp))
	{
		snp <- 1:nsnp
	} else {
		stopifnot(length(snp) == nsnp)
	}
	if(h2 == 0)
	{
		return(dplyr::tibble(snp=snp, beta=0, af=af))
	}
	beta <- stats::rnorm(nsnp, mean=0, sd = sqrt((af * 2 * (1-af))^S))
	vg <- sum(af * 2 * (1-af) * beta^2)
	ve <- (vg - h2 * vg) / h2
	vy <- vg + ve
	beta <- beta / sqrt(vy)
	dat <- dplyr::tibble(snp=snp, beta=beta, af=af)
	return(dat)
}


#' Modify SNP effects to account for LD
#'
#' After generating the set of causal effects at each SNP, use an LD correlation matrix to transform the effects to reflect the correlation structure at the SNPs. Note if running many repeats, only need to generate the LD-modified params once and then can repeatedly re-sample using generate_gwas_ss
#'
#' @param params Output from \code{generate_gwas_params}
#' @param ldobj LD object (e.g. see test_ldobj)
#' @param nthreads Number of threads (can be slow for complete GWAS and large LD regions)
#'
#' @export
#' @return Updated params
add_ld_to_params <- function(params, ldobj, nthreads=1)
{
	stopifnot(is.list(ldobj))
	stopifnot(all(c("map", "ld") %in% names(ldobj[[1]])))
	nchunk <- 1:length(ldobj)
	l <- parallel::mclapply(1:length(ldobj), function(i)
	{
		parami <- dplyr::inner_join(ldobj[[i]][["map"]], params, by="snp")
		parami$beta <- parami$beta %*% ldobj[[i]][["ld"]] %>% drop()
		return(parami)
	}, mc.cores=nthreads) %>%
		dplyr::bind_rows() %>%
		dplyr::arrange(chr, pos)
	return(l)
}


#' Create a GWAS summary dataset
#'
#' Determines SE and generates effect estimates given input parameters
#'
#' @param params Output from \code{generate_gwas_params}
#' @param nid sample size
#' @param vy Variance of trait
#' @param minmaf minimum allowed maf. default=0.01 to prevent instability
#'
#' @export
#' @return list of data frames
generate_gwas_ss <- function(params, nid, vy=1, minmaf=0.01)
{
	nsnp <- nrow(params)
	stopifnot(all(params$af > 0 & params$af < 1))

	params <- params %>%
		dplyr::mutate(
			af = pmax(minmaf, af) %>% pmin(1-minmaf, 1-af),
			se = expected_se(beta, af, nid, vy),
			bhat = sample_beta(beta, se),
			fval = (bhat / se)^2,
			n = nid,
			pval = pf(fval, df1=1, df2=nid-1, lower.tail=FALSE) * 2
		) %>%
		dplyr::select(-beta)
	return(params)
}



#' Create test LD object
#'
#' @param nsnp Number of SNPs
#' @param chunksize Chunksize for splitting
#'
#' @export
#' @return list of chunks, which each contain map and LD matrix
test_ldobj <- function(nsnp, chunksize)
{
	snp <- 1:nsnp
	nchunk <- ceiling(nsnp/chunksize)
	start <- 0:(nchunk-1) * chunksize + 1
	end <- pmin(1:nchunk * chunksize, nsnp)

	ldobj <- lapply(1:nchunk, function(i){
		n <- length(start[i]:end[i])
		p <- qr.Q(qr(matrix(rnorm(n^2), n)))
		Sigma <- crossprod(p, p*(5:1))
		denom <- sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma)))
		rho <- Sigma / denom
		map <- dplyr::tibble(
			snp=start[i]:end[i],
			chr=i,
			pos=1:n,
			ref="A",
			alt="C"
		)
		rownames(rho) <- colnames(rho) <- map$snp
		return(list(map=map, ld=rho))
	})
	return(ldobj)
}



