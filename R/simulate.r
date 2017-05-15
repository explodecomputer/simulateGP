#' Create phenotype from vector of effects and matrix of factors
#'
#' @param effs Vector of effects. If vy and vx are set to default then effect sizes are equivalent to sqrt(variance explained) of each factor
#' @param indep Matrix of factors. Must be same number of columns as length of effs
#' @param vy=1 Variance of phenotype to be outputed
#' @param vx=rep(1, length(effs)) Vector of variances of each factor in indep
#' @param my=0 Mean of phenotype to be outputted
#'
#' @export
#' @return Phenotype vector
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
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
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
#' @param sqrt=TRUE Output effect sizes in terms of standard deviations
#'
#' @export
#' @return Vector of effects
choose_effects <- function(nsnp, totvar, sqrt=TRUE)
{
	eff <- rnorm(nsnp)
	eff <- sign(eff) * eff^2
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out)
}

#' Convert continuous trait to binary 
#'
#' @param y Phenotype vector
#' @param prevalence=NULL Disease prevalence
#' @param threshold=NULL Disease threshold
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

