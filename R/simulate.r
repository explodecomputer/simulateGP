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

