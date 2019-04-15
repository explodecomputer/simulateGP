#' Get summary statistics in simple linear regression
#'
#' @param y Vector of dependent variable
#' @param x Vector of independent variable
#'
#' @export
#' @return List
fast_assoc <- function(y, x)
{
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lower.tail=FALSE)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p, n=n
	))
}


#' For known effect size, sample size and variances, what is the expected association?
#'
#' @param eff Vector of effect sizes
#' @param n Vector of sample sizes
#' @param vx Vector of variances of x
#' @param vy vector of varinces of y
#'
#' @export
#' @return Data frame of associations
expected_gwas <- function(eff, n, vx, vy)
{
	rsq <- eff^2 * vx / vy
	varexp_xy <- rsq * vy
	ve <- varexp_xy / rsq - varexp_xy
	se <- sqrt(ve / ((n-2) * vx))
	tval <- eff / se
	p <- pt(abs(tval), n-1, lower.tail=FALSE)
	dat <- tibble::data_frame(bhat=eff, se=se, fval=tval^2, pval=p, n=n)
	return(dat)
}


#' Perform association of many SNPs against phenotype
#'
#' @param y Vector of phenotypes
#' @param g Matrix of genotypes
#'
#' @export
#' @return Data frame
gwas <- function(y, g)
{
	out <- matrix(0, ncol(g), 6)
	for(i in 1:ncol(g))
	{
		o <- fast_assoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- tibble::as_data_frame(out)
	names(out) <- names(o)
	out$snp <- 1:ncol(g)
	return(out)
}

#' Get effs for two traits and make dat format
#'
#' @param x Vector of exposure phenotype
#' @param y Vector of outcome phenotype
#' @param g Matrix of genotypes
#' @param xname
#' @param yname
#'
#' @export
#' @return Data frame
get_effs <- function(x, y, g, xname="X", yname="Y")
{
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)
	return(make_dat(gwasx, gwasy))
}

#' Organise outputs from \code{gwas} into harmonised dat format
#'
#' @param gwasx Output from \code{gwas}
#' @param gwasy Output from \code{gwas}
#' @param xname
#' @param yname
#'
#' @export
#' @return data frame
make_dat <- function(gwasx, gwasy, xname="X", yname="Y")
{
	d <- dplyr::inner_join(gwasx, gwasy, by='snp')
	dat <- tibble::data_frame(
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
		exposure=rep(paste0("x", 1:length(exposures)), each=ncol(g)),
		id.exposure=rep(paste0("x", 1:length(exposures)), each=ncol(g)),
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
