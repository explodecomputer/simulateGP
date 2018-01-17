# Also see https://github.com/explodecomputer/pic_haps/

correlated_binomial <- function(nid, p1, p2, rho, n=2, round=TRUE, print=FALSE)
{
	# from https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
	p <- p1
	q <- p2
	a <- function(rho, p, q) {
		rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
	}

	a.0 <- a(rho, p, q)
	prob <- c(`(0,0)`=a.0, `(1,0)`=1-q-a.0, `(0,1)`=1-p-a.0, `(1,1)`=a.0+p+q-1)
	if (min(prob) < 0) {
		print(prob)
		stop("Error: a probability is negative.")
	}
	#
	# Illustrate generation of correlated Binomial variables.
	#
	n.sim <- nid
	u <- sample.int(4, n.sim * n, replace=TRUE, prob=prob)
	y <- floor((u-1)/2)
	x <- 1 - u %% 2
	x <- colSums(matrix(x, nrow=n)) # Sum in groups of `n`
	y <- colSums(matrix(y, nrow=n)) # Sum in groups of `n`

	if(round)
	{
		x <- round(x)
		y <- round(y)
	}

	if(print)
	{
		print(table(x, y))
		print(cor(x, y))
	}
	return(cbind(x, y))

	#
	# Plot the empirical bivariate distribution.
	#
	# plot(x+rnorm(length(x), sd=1/8), y+rnorm(length(y), sd=1/8),
	#      pch=19, cex=1/2, col="#00000010",
	#      xlab="X", ylab="Y",
	#      main=paste("Correlation is", signif(cor(x,y), 3)))
	# abline(v=mean(x), h=mean(y), col="Red")
	# abline(lm(y ~ x), lwd=2, lty=3)
}

hap_freqs <- function(r, p1, p2)
{
	# d = pAB - p1p2
	# d = p1q2 - pAb
	# d = q1p2 - paB
	# d = pab - q1q2
	# r = d / denom
	denom <- sqrt(p1 * (1-p1) * p2 * (1-p2))
	p <- c(
		`(1,1)` = p1 * p2 + r * denom, 
		`(0,0)` = (1-p1) * (1-p2) + r * denom, 
		`(1,0)` = p1 * (1-p2) - r * denom, 
		`(0,1)` = (1-p1) * p2 - r * denom
	)
	stopifnot(all(p >= 0))
	return(p)
}

test_hap_freqs <- function(r, p1, p2)
{
	a <- try(hap_freqs(r, p1, p2), silent=TRUE)
	return(ifelse(class(a)=='try-error', FALSE, TRUE))
}

simulate_haplotypes <- function(nid, r, p1, p2)
{
	p <- hap_freqs(r, p1, p2)
	n <- round(p * nid)
	diff <- sum(n) - nid
	i <- 1
	while(diff != 0)
	{
		if(diff > 0)
		{
			n[i] <- n[i] - 1
		} else {
			n[i] <- n[i] + 1
		}
		diff <- sum(n) - nid
		i <- ifelse(i == 4, 1, i + 1)
	}
	mat <- rbind(
		cbind(rep(1, n[1]), rep(1, n[1])),
		cbind(rep(0, n[2]), rep(0, n[2])),
		cbind(rep(1, n[3]), rep(0, n[3])),
		cbind(rep(0, n[4]), rep(1, n[4]))
	)
	ind1 <- sample(1:nid, nid, replace=FALSE)
	ind2 <- sample(1:nid, nid, replace=FALSE)
	haps <- cbind(
		A1 = mat[ind1, 1],
		A2 = mat[ind2, 1],
		B1 = mat[ind1, 2],
		B2 = mat[ind2, 2]
	)
	return(haps)
}

simulate_geno <- function(nid, r, p1, p2)
{
	a <- simulate_haplotypes(nid, r, p1, p2)
	b <- cbind(a[,1]+a[,2], a[,3]+a[,4])
	return(b)
}

epistasis_problem.simulate_system <- function(nid, r, p1, p2, p3, rsq)
{
	cis <- simulate_geno(nid, r, p1, p2)
	trans <- rbinom(nid, 2, p3)
	y <- scale(cis[,1]) * sqrt(rsq) + rnorm(nid, sd=sqrt(1-rsq))
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

	hist(param$pintt)
	min(param$pintt)

	hist(param$pint, breaks=20)
	min(param$pint)

}
