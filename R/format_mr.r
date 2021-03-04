#' Organise outputs from \code{gwas} into harmonised dat format
#'
#' @param gwasx Output from \code{gwas}
#' @param gwasy Output from \code{gwas}
#' @param xname exposure name
#' @param yname outcome name
#'
#' @export
#' @return data frame
merge_exp_out <- function(gwasx, gwasy, xname="X", yname="Y")
{
	d <- dplyr::inner_join(gwasx, gwasy, by='snp')
	dat <- dplyr::tibble(
		SNP = d$snp,
		exposure=xname,
		id.exposure=xname,
		outcome=yname,
		id.outcome=yname,
		beta.exposure=d$bhat.x,
		beta.outcome=d$bhat.y,
		se.exposure=d$se.x,
		se.outcome=d$se.y,
		pval.exposure=d$pval.x,
		pval.outcome=d$pval.y,
		samplesize.exposure=d$n.x,
		samplesize.outcome=d$n.y,
		units.exposure = "SD",
		units.outcome = "SD",
		rsq.exposure = d$fval.x / (d$fval.x + d$n.x - 2),
		rsq.outcome = d$fval.y / (d$fval.y + d$n.y - 2),
		mr_keep=TRUE
	)
	return(dat)	
}

#' Simple recoding to have every effect on x positive
#'
#' @param dat Output from get_effs
#'
#' @export
#' @return Data frame
recode_dat_simple <- function(dat)
{
	.Deprecated('recode_dat')
	sign0 <- function(x) {
		x[x == 0] <- 1
		return(sign(x))
	}
	index <- sign0(dat$beta.exposure) == -1
	dat$beta.exposure <- abs(dat$beta.exposure)
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
	return(dat)
}

#' Intercept recoding to have every effect on x positive
#'
#' Tries to avoid issue of recoding by finding intercept and pivoting negative g-x associations around intercept
#'
#' @param dat Output from get_effs
#'
#' @export
#' @return Data frame
recode_dat_intercept <- function(dat)
{
	.Deprecated('recode_dat')
	a <- lm(beta.outcome ~ beta.exposure, dat)$coefficients[1]
	index <- dat$beta.exposure < 0
	dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1 + 2 * a
	dat$index <- index
	return(dat)		
}

#' Recode data to make every effect on x positive
#'
#' Can use simple method or by pivoting around intercept
#'
#' @param dat Output from get_effs
#' @param method Default 'intercept'. Alternatively can specify 'simple'
#'
#' @export
#' @return Data frame
recode_dat <- function(dat, method='intercept')
{
	if(method == 'intercept')
	{
		a <- lm(beta.outcome ~ beta.exposure, dat)$coefficients[1]
		index <- dat$beta.exposure < 0
		dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
		dat$beta.outcome[index] <- dat$beta.outcome[index] * -1 + 2 * a
		dat$index <- index
		return(dat)		
	} else if(method == 'simple') {
		sign0 <- function(x) {
			x[x == 0] <- 1
			return(sign(x))
		}
		index <- sign0(dat$beta.exposure) == -1
		dat$beta.exposure <- abs(dat$beta.exposure)
		dat$beta.outcome[index] <- dat$beta.outcome[index] * -1
		return(dat)
	} else {
		stop('method must be intercept or simple')
	}
}


#' Take several exposures and one outcome and make the data required for multivariable MR
#'
#' @param exposures List of exposure vectors
#' @param y Vector of outcomes
#' @param g Matrix of genotypes
#'
#' @export
#' @return mv_harmonise_data output
make_mvdat <- function(exposures, y, g)
{
	stopifnot(is.list(exposures))
	message("There are ", length(exposures), " exposures")
	if(is.null(names(exposures)))
	{
		names(exposures) <- paste0("x", 1:length(exposures))
	}
	exposure_dat <- lapply(exposures, function(x) gwas(x, g))
	# exposure_dat1 <- gwas(x1, g)
	# exposure_dat2 <- gwas(x2, g)
	af <- colSums(g) / (nrow(g) * 2)
	out <- gwas(y, g)
	mvexp <- data.frame(
		SNP=1:ncol(g),
		effect_allele.exposure="A",
		other_allele.exposure="G",
		eaf.exposure=rep(af, times=length(exposures)),
		exposure=rep(names(exposures), each=ncol(g)),
		id.exposure=rep(names(exposures), each=ncol(g)),
		beta.exposure = lapply(exposure_dat, function(x) x$bhat) %>% unlist,
		se.exposure = lapply(exposure_dat, function(x) x$se) %>% unlist,
		pval.exposure = lapply(exposure_dat, function(x) x$pval) %>% unlist
	)
	outcome_dat <- data.frame(
		SNP=rep(1:ncol(g), times=length(exposures)),
		outcome="y",
		id.outcome="y",
		effect_allele.outcome="A",
		other_allele.outcome="G",
		eaf.outcome=rep(af, times=length(exposures)),
		beta.outcome = rep(out$bhat, times=length(exposures)),
		se.outcome = rep(out$se, times=length(exposures)),
		pval.outcome = rep(out$pval, times=length(exposures))
	)
	mvdat <- TwoSampleMR::mv_harmonise_data(mvexp, outcome_dat)
	return(mvdat)
}
