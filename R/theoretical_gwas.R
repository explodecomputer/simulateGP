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
#' @param map Data frame containing at least `af` allele frequency and `snp` SNP columns. SNPs must be unique, `af` must be between 0 and 1. Optionally also include the chr, pos, ref, alt columns if using LD-aware simulations
#' @param h2 Variance explained by all SNPs
#' @param S Selection coefficient on trait. Default = 0
#' @param Pi Proportion of variants that have an effect - sampled randomly. Default=1
#'
#' @export
#' @return data frame
generate_gwas_params <- function(map, h2, S=0, Pi=1)
{
	stopifnot(all(c("snp", "af") %in% names(map)))
	stopifnot(all(map$af > 0 & map$af < 1))
	stopifnot(!any(duplicated(map$snp)))
	nsnp <- nrow(map)
	if(h2 == 0)
	{
		map$beta <- 0
		return(map)
	}

	index <- sample(1:nsnp, ceiling(nsnp * Pi), replace=FALSE)

	map$beta <- 0
	map$beta[index] <- stats::rnorm(length(index), mean=0, sd = sqrt((map$af[index] * 2 * (1-map$af[index]))^S))
	vg <- sum(map$af * 2 * (1-map$af) * map$beta^2)
	ve <- (vg - h2 * vg) / h2
	vy <- vg + ve
	map$beta <- map$beta / sqrt(vy)
	return(map)
}


#' Modify SNP effects to account for LD
#'
#' After generating the set of causal effects at each SNP, use an LD correlation matrix to transform the effects to reflect the correlation structure at the SNPs. Note if running many repeats, only need to generate the LD-modified params once and then can repeatedly re-sample using generate_gwas_ss
#'
#' @param params Output from \code{generate_gwas_params}
#' @param ldobj LD objects (e.g. see test_ldobj)
#' @param ldobjlist List of LD objects 
#' @param ldobjfiles Array of filenames containing LD object files (e.g. see \code{generate_ldobj})
#' @param ldobjdir Directory containing output from \code{generate_ldobj}
#' @param nthreads Number of threads (can be slow for complete GWAS and large LD regions)
#'
#' @export
#' @return Updated params
add_ld_to_params <- function(params, ldobj=NULL, ldobjlist=NULL, ldobjfiles=NULL, ldobjdir=NULL, nthreads=1)
{
	params <- subset(params, !duplicated(snp))
	stopifnot(c("chr", "pos", "snp", "beta") %in% names(params))
	if(!is.null(ldobj))
	{
		stopifnot(is.list(ldobj))
		stopifnot(all(c("map", "ld") %in% names(ldobj)))
		params <- subset(params, snp %in% ldobj$map$snp)
		if(nrow(params) > 0)
		{
			params$beta <- params$beta %*% ldobj[["ld"]] %>% drop()
		}
		return(params)
	}

	if(!is.null(ldobjlist))
	{
		stopifnot(is.list(ldobjlist))
		stopifnot(all(c("map", "ld") %in% names(ldobjlist[[1]])))
		nchunk <- 1:length(ldobjlist)
		l <- pbapply::pblapply(1:length(ldobjlist), function(i)
		{
			parami <- subset(params, snp %in% ldobjlist[[i]][["map"]]$snp)
			if(nrow(parami) > 0)
			{
				parami$beta <- parami$beta %*% ldobjlist[[i]][["ld"]] %>% drop()
			}
			return(parami)
		}, cl=nthreads) %>%
			dplyr::bind_rows() %>%
			dplyr::arrange(chr, pos)
		return(l)
	}

	if(!is.null(ldobjfiles))
	{
		stopifnot(is.character(ldobjfiles))
		stopifnot(all(file.exists(ldobjfiles)))
		nchunk <- 1:length(ldobjfiles)
		l <- pbapply::pblapply(1:length(ldobjfiles), function(i)
		{
			ldobj <- readRDS(ldobjfiles[i])
			parami <- subset(params, snp %in% ldobj[["map"]]$snp)
			if(nrow(parami) > 0)
			{
				parami$beta <- parami$beta %*% ldobj[["ld"]] %>% drop()
			}
			return(parami)
		}, cl=nthreads) %>%
			dplyr::bind_rows() %>%
			dplyr::arrange(chr, pos)
		return(l)
	}

	if(!is.null(ldobjdir))
	{
		fn <- list.files(ldobjdir, full.names=TRUE) %>% grep("ldobj_chr", .data, value=TRUE)
		message("Found ", length(fn), " region files")
		message("Splitting map by region")
		code <- fn %>%
			basename() %>%
			gsub("ldobj_chr", "", .data) %>%
			gsub("\\.rds", "", .data)

		l <- pbapply::pblapply(1:length(code), function(i)
		{
			ldobj <- readRDS(fn[i])
			parami <- subset(params, region == code[i])
			parami <- subset(parami, snp %in% ldobj[["map"]]$snp)
			if(nrow(parami) > 0)
			{
				parami$beta <- parami$beta %*% ldobj[["ld"]] %>% drop()
			}
			return(parami)
		}, cl=nthreads) %>%
			dplyr::bind_rows() %>%
			dplyr::arrange(chr, pos)
		return(l)
	}

	stop("At least one of ldobj, ldobjlist, ldobjfiles, ldobjdir should be specified")
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


