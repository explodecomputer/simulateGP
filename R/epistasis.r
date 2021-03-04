epistasis_problem.simulate_system <- function(nid, r, p1, p2, p3, rsq)
{
	cis <- simulate_geno(nid, r, p1, p2)
	trans <- rbinom(nid, 2, p3)
	y <- scale(cis[,1]) * sqrt(rsq) + stats::rnorm(nid, sd=sqrt(1-rsq))
	dat <- data.frame(
		cis=cis[,1],
		cist=cis[,2],
		trans=trans,
		y=y
	)
	return(dat)
}

epistasis_problem.run1 <- function(param, i)
{
	set.seed(param$seed[i])
	a <- with(param[i, ], simulate_system(nid, r, p1, p2, p3, rsq))
	moda1 <- lm(y ~ as.factor(cis) + as.factor(trans), a)
	moda2 <- lm(y ~ as.factor(cis) * as.factor(trans), a)
	moda <- anova(moda1, moda2)
	modb1 <- lm(y ~ as.factor(cist) + as.factor(trans), a)
	modb2 <- lm(y ~ as.factor(cist) * as.factor(trans), a)
	modb <- anova(modb1, modb2)
	param$pint[i] <- moda$P[2]
	param$pintt[i] <- modb$P[2]
	return(param)
}

epistasis_problem <- function()
{
	a <- simulate_haplotypes(563, 0.9, 0.5, 0.5)
	b <- simulate_geno(1000, 0.9, 0.5, 0.5)
	a <- correlated_binomial(10000, 0.5, 0.5, 0.1)
	a <- simulate_system(1000, 0.8, 0.5, 0.5, 0.5, 0.3)

	param <- expand.grid(
		sim=1:300,
		nid=1000,
		r=c(sqrt(0.1), sqrt(0.5), sqrt(0.9)),
		p1=0.5,
		p2=0.5,
		p3=0.5,
		rsq=c(0.01, 0.1, 0.5),
		seed=NA,
		pint=NA,
		pintt=NA
	)
	param$seed <- 1:nrow(param)


	for(i in 1:nrow(param))
	{
		message(i)
		param <- run1(param, i)
	}

	graphics::hist(param$pintt)
	min(param$pintt)

	graphics::hist(param$pint, breaks=20)
	min(param$pint)

}
