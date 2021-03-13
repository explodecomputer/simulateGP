context("fastassoc")
library(simulateGP)

test_that("fastassoc", {
	set.seed(1234)
	n <- 10000
	x <- rnorm(n)
	y <- x + rnorm(n)
	a <- fast_assoc(y, x)
	b <- summary(lm(y ~ x))
	expect_equal(a$bhat, b$coefficients[2,1])
	expect_equal(a$se, b$coefficients[2,2])
	expect_equal(a$pval, b$coefficients[2,4])
})


test_that("gwas", {
	set.seed(1234)
	nsnp <- 3
	nid <- 10000
	rsq_gx <- 0.05
	u <- rnorm(nid)
	g <- make_geno(nid=nid, nsnp=nsnp, af=0.5)
	effs <- choose_effects(nsnp=nsnp, totvar=rsq_gx)
	x <- make_phen(effs=c(effs, 0.3), indep=cbind(g, u))
	res <- gwas(x, g)	
	expect_equal(nrow(res), nsnp)
	mod <- summary(lm(x ~ g[,1]))
	expect_equal(res$bhat[1], mod$coefficients[2,1])
})


test_that("get_effs", {

	set.seed(1234)
	beta_xy <- -0.3
	nsnp <- 3
	nid <- 10000
	rsq_gx <- 0.05
	u <- rnorm(nid)
	g <- make_geno(nid=nid, nsnp=nsnp, af=0.5)
	effs <- choose_effects(nsnp=nsnp, totvar=rsq_gx)
	x <- make_phen(effs=c(effs, 0.3), indep=cbind(g, u))
	y <- make_phen(effs=c(beta_xy, 0.3), cbind(x, u))
	res <- gwas(x, g)
	expect_equal(nrow(res), nsnp)
	mod <- summary(lm(x ~ g[,1]))
	expect_equal(res$bhat[1], mod$coefficients[2,1])

	res <- get_effs(x, y, g)
	expect_equal(nrow(res), nsnp)
})


test_that("logistic", {

	set.seed(1234)
	beta_xy <- -0.3
	nsnp <- 3
	nid <- 10000
	rsq_gx <- 0.05
	u <- rnorm(nid)
	g <- make_geno(nid=nid, nsnp=nsnp, af=0.5)
	effs <- choose_effects(nsnp=nsnp, totvar=rsq_gx)
	x <- make_phen(effs=c(effs, 0.3), indep=cbind(g, u))
	z <- rbinom(nid, 1, plogis(x))
	res <- gwas(z, g, logistic=TRUE)

	expect_equal(nrow(res), nsnp)
	mod <- summary(glm(z ~ g[,1], family="binomial"))
	expect_equal(res$bhat[1], mod$coefficients[2,1])
})
