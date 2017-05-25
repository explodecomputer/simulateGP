#' Obtain 2x2 contingency table from marginal parameters and odds ratio
#'
#' Columns are the case and control frequencies
#' Rows are the frequencies for allele 1 and allele 2
#'
#' @param af Allele frequency of effect allele 
#' @param prop Proportion of cases
#' @param odds_ratio Odds ratio
#' @param eps=1e-15 tolerance
#'
#' @export
#' @return 2x2 contingency table as matrix
contingency <- function(af, prop, odds_ratio, eps=1e-15)
{
	a <- odds_ratio-1
	b <- (af+prop)*(1-odds_ratio)-1
	c_ <- odds_ratio*af*prop

	if (abs(a) < eps)
	{
		z <- -c_ / b
	} else {
		d <- b^2 - 4*a*c_
		if (d < eps*eps) 
		{
			s <- 0
		} else {
			s <- c(-1,1)
		}
		z <- (-b + s*sqrt(max(0, d))) / (2*a)
	}
	y <- vapply(z, function(a) zapsmall(matrix(c(a, prop-a, af-a, 1+a-af-prop), 2, 2)), matrix(0.0, 2, 2))
	i <- apply(y, 3, function(u) all(u >= 0))
	return(y[,,i])
}

#' Estimate allele frequency from SNP
#'
#' @param g Vector of 0/1/2
#'
#' @export
#' @return Allele frequency 
allele_frequency <- function(g)
{
	(sum(g == 1) + 2 * sum(g == 2)) / (2 * sum(!is.na(g)))
}


#' Estimate the allele frequency in population from case/control summary data
#'
#' @param af Effect allele frequency (or MAF)
#' @param prop Proportion of samples that are cases
#' @param odds_ratio Odds ratio
#' @param prevalence Population disease prevalence
#'
#' @export
#' @return Population allele frequency
get_population_allele_frequency <- function(af, prop, odds_ratio, prevalence)
{
	co <- contingency(af, prop, odds_ratio)
	af_controls <- co[1,2] / (co[1,2] + co[2,2])
	af_cases <- co[1,1] / (co[1,1] + co[2,1])
	af <- af_controls * (1 - prevalence) + af_cases * prevalence
	return(af)
}


#' Estimate proportion of variance of liability explained by SNP in general population
#'
#' This uses equation 10 in Genetic Epidemiology 36 : 214–224 (2012)
#' 
#' @param b Log odds ratio
#' @param af 
#' @param ncase <what param does>
#' @param ncontrol <what param does>
#' @param prevalence <what param does>
#' @param model Is the effect size estiamted in "logit" (default) or "probit" model
#'
#' @export
#' @return Rsq
lor_to_rsq <- function(b, af, ncase, ncontrol, prevalence, model="logit")
{
	if(model == "logit")
	{
		ve <- pi^2/3
	} else if(model == "probit") {
		ve <- 1
	} else {
		stop("Model must be probit or logit")
	}
	af <- get_population_allele_frequency(af, ncase / (ncase + ncontrol), exp(b), prevalence)
	vg <- b^2 * af * (1-af)
	return(vg / (vg + ve) / 0.58)
}

