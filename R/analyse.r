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
	dat <- data.frame(bhat=eff, se=se, fval=tval^2, pval=p, n=n)
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
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

#' Get effs for two traits and make dat format
#'
#' @param x Vector of exposure phenotype
#' @param y Vector of outcome phenotype
#' @param g Matrix of genotypes
#'
#' @export
#' @return Data frame
get_effs <- function(x, y, g)
{
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)

	dat <- data.frame(
		exposure="X",
		id.exposure="X",
		outcome="Y",
		id.outcome="Y",
		beta.exposure=gwasx$bhat,
		beta.outcome=gwasy$bhat,
		se.exposure=gwasx$se,
		se.outcome=gwasy$se,
		pval.exposure=gwasx$pval,
		pval.outcome=gwasy$pval,
		samplesize.exposure=gwasx$n,
		samplesize.outcome=gwasy$n,
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
	a <- lm(beta.outcome ~ beta.exposure, dat)$coefficients[1]
	index <- dat$beta.exposure < 0
	dat$beta.exposure[index] <- dat$beta.exposure[index] * -1
	dat$beta.outcome[index] <- dat$beta.outcome[index] * -1 + 2 * a
	dat$index <- index
	return(dat)
}

