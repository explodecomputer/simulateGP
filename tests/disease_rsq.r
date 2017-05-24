library(tidyverse)
library(devtools)
load_all()

contingency <- function(maf, prop, odds_ratio, eps=1e-15)
{
	a <- odds_ratio-1
	b <- (maf+prop)*(1-odds_ratio)-1
	c_ <- odds_ratio*maf*prop

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
	y <- vapply(z, function(a) zapsmall(matrix(c(a, prop-a, maf-a, 1+a-maf-prop), 2, 2)), matrix(0.0, 2, 2))
	i <- apply(y, 3, function(u) all(u >= 0))
	return(y[,,i])
}

allele_frequency <- function(g)
{
	(sum(g == 1) + 2 * sum(g == 2) / (2 * sum(!is.na(g)))
}

allele_frequency(g)
allele_frequency(ga)



get_population_allele_frequency <- function(maf, prop, odds_ratio, prevalence)
{
	co <- contingency(maf, prop, odds_ratio)
	af_controls <- co[1,2] / (co[1,2] + co[2,2])
	af_cases <- co[1,1] / (co[1,1] + co[2,1])
	af <- af_controls * (1 - prevalence) + af_cases * prevalence
	return(af)
}


lor_to_rsq <- function(b, af, ncase, ncontrol, prevalence)
{
	af <- get_population_allele_frequency(af, ncase / (ncase + ncontrol), exp(b), prevalence)
	vg <- b^2 * af * (1-af)
	return(vg / (vg + pi^2/3) / 0.58)
}


param <- data_frame(
	va = runif(100, 0.001, 0.0001),
	n = 100000,
	prevalence = runif(100, 0.01, 0.1),
	proportion = runif(100, 0.2, 0.8),
	maf = runif(100, 0.05, 0.5),
	rsq = NA
)


for(i in 1:nrow(param))
{
	message(i)
	va <- param$va[i]
	maf <- param$maf[i]
	prevalence <- 0.1
	nid <- param$n[i]
	proportion <- param$proportion[i]
	g <- make_geno(nid, 1, maf)
	x <- make_phen(sqrt(va), g)
	bin <- y_to_binary(x, prevalence=prevalence)
	asc <- ascertain_samples(bin, proportion)
	ga <- g[asc, ]
	xa <- x[asc]
	bina <- bin[asc]
	a <- fast_assoc(xa, ga)
	mod <- summary(glm(bina ~ ga, family="binomial"))
	eff <- coefficients(mod)[2,1]
	pval <- coefficients(mod)[2,4]
	param$rsq[i] <- lor_to_rsq(eff, allele_frequency(ga), sum(bina==1), sum(bina==0), prevalence)
	param$paf[i] <- get_population_allele_frequency(allele_frequency(ga), sum(bina==1) / length(bina), exp(eff), prevalence)
	param$pval[i] <- pval
}

plot(rsq ~ va, param)
summary(lm(rsq ~ va + prevalence + proportion + maf, param))
plot(-log10(param$pval))

ggplot(param, aes(x=va, y=rsq)) +
geom_point(aes(colour=-log10(pval))
)
get_population_allele_frequency()

paf <- get_population_allele_frequency(allele_frequency(ga), sum(bina==1) / length(bina), exp(eff), prevalence)
pvg <- eff^2 * paf * (1-paf)
saf <- allele_frequency(ga)
svg <- eff^2 * saf * (1-saf)
pvg / (pvg + pi^2/3)
svg / (svg + pi^2/3)

summary(lm(rsq ~ va, param))

vg <- eff^2 * 2 * allele_frequency(g) * (1 - allele_frequency(g))
vg / (vg + pi^2/3)

vg <- eff^2 * 2 * allele_frequency(ga) * (1 - allele_frequency(ga))
vg / (vg + pi^2/3)


co <- contingency(allele_frequency(ga), sum(bina==1)/length(bina), exp(eff))
allele_frequency(ga)
sum(co[1,])
co
sum(co[,1])

allele_frequency(g[bin==0])





contingency(0.1, 0.2, exp(coefficients(mod)[2]))






freq <- sum(ga) / (2 * length(ga))

rsq <- mod$coefficients[2]^2 * freq * (1 - freq) * 2 / var_cc(0.1, 0.9)
rsq


	h2l


mod <- lm(xa ~ ga)
rsq <- mod$coefficients[2]^2 * freq * (1 - freq) * 2 / var_cc(0.8, 0.9)
rsq



summary(lm(x ~ g))
summary(glm(bin ~ g, family="binomial"))
var_cc(prevalence, prevalence)


res <- array(0, 50)
proportion <- seq(0.1, 0.9, length.out=50)
for(i in 1:50)
{
	message(i)
	va <- 0.1
	maf <- 0.5
	prevalence <- 0.1
	nid <- 100000

	g <- make_geno(nid, 1, maf)
	x <- make_phen(sqrt(va), g)

	d <- risk_simulation(g, eff=0.1, prevalence=0.1, prop_discovered=1)
	bin <- d$disease

	asc <- ascertain_samples(bin, proportion[i])
	ga <- g[asc, ]
	bina <- bin[asc]
	mod <- glm(bina ~ ga, family="binomial")
	res[i] <- mod$coefficients[2]
}


plot(res ~ proportion)
summary(glm(bin ~ g, family="binomial"))




summary(glm(bin ~ g, family="binomial"))



0.9 * 0.1
table(bin)
var(bin)





var_cc <- function(proportion, prevalence)
{
	tval <- qnorm(prevalence, lower.tail=FALSE)
	zval <- dnorm(tval)
	i <- zval / prevalence
	lambda <- (proportion - prevalence) / (1 - prevalence)
	vp <- 1 + i * lambda * (tval - i * lambda)
	return(vp)
}





