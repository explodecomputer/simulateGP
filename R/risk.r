#' Scale variable to have range between 0 and 1
#'
#' @param x Vector
#'
#' @export
#' @return Vector
range01 <- function(x)
{
	(x-min(x))/(max(x)-min(x))
}

#' Translate risk from liability to probability scale
#'
#' @param Genetic score on liability scale 
#' @param h2x Disease heritability on liability scale
#' @param prev Prevalence
#'
#' @export
#' @return Vector of disease probabilities
gx_to_gp <- function(gx, h2x, prev)
{
	x_prime <- qnorm(prev, 0, 1)
	p <- pnorm(x_prime, mean=gx, sd = sqrt(1 - h2x), lower.tail=FALSE)
	return(p)
}

#' Plot liability vs probability disease risk
#'
#' @param x Disease risk on liability scale
#' @param o Disease risk on probability scale
#' @param xlab="Values (low to high)" Xlab
#' @param ylab="" Ylab
#' @param title="" Title
#' @param xname="GRS" Name of liability
#' @param oname="Disease" Name of disease
#'
#' @export
#' @return ggplot
risk_cross_plot <- function(x, o, xlab="Values (low to high)", ylab="", title="", xname="GRS", oname="Disease")
{
	require(ggplot2)
	d <- data.frame(
		value = c(range01(o), range01(x)),
		key = c(rep(oname, length(o)), rep(xname, length(x))),
		gr = rep(1:length(x), times=2)
	)
	d$key <- factor(d$key, levels=c(xname, oname))
	ggplot(d, aes(x=value, y=key)) +
	geom_line(aes(group=gr), alpha=0.1) +
	geom_point(aes(colour=key)) +
	labs(x=xlab,y=ylab,title=title) +
	scale_colour_discrete(guide=FALSE) +
	theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())
}

#' Make simulation to compare disease and liability scales
#'
#' Compares the liability and probability of disease under two scenarios - where all SNPs are known, and where only some SNPs are known
#'
#' @param G Matrix of genotypes
#' @param eff SNP effects on liability scale
#' @param prevalence Disease prevalence
#' @param prop_discovered Proportion of SNPs discovered
#'
#' @export
#' @return Data frame
risk_simulation <- function(G, eff, prevalence, prop_discovered)
{
	nid <- nrow(G)
	nsnp <- ncol(G)
	h2x <- var(eff)
	gx_true <- as.numeric(scale(G %*% eff)) * sqrt(h2x)
	prob_disease <- gx_to_gp(gx_true, h2x, 1-prevalence)
	disease <- rbinom(nid, 1, prob_disease)
	eff_pred <- eff
	eff_pred[sample(1:nsnp, nsnp * (1-prop_discovered))] <- 0
	gx_pred <- as.numeric(G %*% eff_pred / sqrt(nsnp))
	dat <- data.frame(gx_true=gx_true, gx_pred=gx_pred, prob_disease=prob_disease, disease=disease)
	return(dat)
}

