
#' Simulate variable based on the influences of other variables
#'
#' I need a function that takes some variables, and an effect size for each of them, and this creates a new variable where y = Xb + e, where X is the matrix of independent variables, b is their effect sizes, and e is a noise term. By default it's easier if b relates to the variance explained by each vector in X, and e has variance of 1 - sum(b).
#'
#' @param effs An array of effect sizes that the variables in indep have on the variable that you are simulating
#' @param indep A matrix of variables, rows = Samples and columns = variables that have an influence on the variable that you are simulating
#' @param vy What variance the output should have. Default = 1, meaning effs relate to the signed rsq of the influence of the indep variables on the outcome
#' @param vx What variance the indep variables should have. Default is to set to 1, meaning that the effects are the signed variance explained
#' @param my mean value of y to be output
#'
#' @export
#' @return Vector of y values
make_phen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0)
{
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	sc <- sum(cors^2)
	if(sc >= 1)
	{
		print(sc)
		stop("effects explain more than 100% of variance")
	}
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, stats::rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

#' Create genotype matrix
#'
#' @param nid Number of samples
#' @param nsnp Number of SNPs
#' @param af Allele frequency of all SNPs (all SNPs have the same allele frequency) 
#'
#' @export
#' @return Matrix of genotypes, rows = individuals, columns = snps
make_geno <- function(nid, nsnp, af)
{
	return(matrix(rbinom(nid * nsnp, 2, af), nid, nsnp))
}

#' Get vector of effects that explain some amount of variance
#'
#' @param nsnp Number of SNPs
#' @param totvar Total variance explained by all SNPs
#' @param sqrt Output effect sizes in terms of standard deviations. Default=TRUE
#' @param mua Constant term to be added to effects. Default = 0
#'
#' @export
#' @return Vector of effects
choose_effects <- function(nsnp, totvar, sqrt=TRUE, mua=0)
{
	eff <- stats::rnorm(nsnp)
	eff <- sign(eff) * eff^2
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out + mua)
}

#' Convert continuous trait to binary 
#'
#' @param y Phenotype vector
#' @param prevalence Disease prevalence. Default = NULL
#' @param threshold Disease threshold Default = NULL
#'
#' @export
#' @return Vector of binary trait
y_to_binary <- function(y, prevalence=NULL, threshold=NULL)
{
	if(is.null(prevalence) & is.null(threshold)) stop("Prevalence or threshold needs to be non-null")
	if(!is.null(prevalence))
	{
		d <- y
		t <- quantile(d, 1-prevalence)
		d[y >= t] <- 1
		d[y < t] <- 0
		return(d)
	}
	if(!is.null(threshold))
	{
		d <- y
		d[y >= threshold] <- 1
		d[y < threshold] <- 0
		return(d)
	}
}


#' Ascertain some proportion of cases and controls from binary phenotype
#'
#' @param d Vector of 1/0
#' @param prop_cases Proportion of 1s to retain
#'
#' @export
#' @return Array of IDs
ascertain_samples <- function(d, prop_cases)
{
	d <- d[!is.na(d)]
	d <- d[d %in% c(0,1)]
	x1 <- sum(d==1)
	x0 <- sum(d==0)
	exp_cases <- x1
	exp_controls <- (x1 - prop_cases * x1) / prop_cases
	if(round(exp_controls) > x0)
	{
		exp_controls <- x0
		exp_cases <- (x0 - (1 - prop_cases) * x0) / (1 - prop_cases)
	}
	i0 <- which(d == 0)
	i0 <- sample(i0, exp_controls, replace=FALSE)
	i1 <- which(d == 1)
	i1 <- sample(i1, exp_cases, replace=FALSE)

	return(sort(c(i0, i1)))
}

