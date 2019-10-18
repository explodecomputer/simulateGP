context("Theoretical GWAS")
library(simulateGP)

test_that("expected se works", {

	set.seed(1)
	h2 <- 0.3
	nid <- 1000
	nsnp <- 1
	eff <- rnorm(nsnp)
	g <- make_geno(1000, 1, 0.5)
	grs <- g %*% eff
	vare <- var(g) * (1 - h2) / h2
	y <- grs + rnorm(nid, 0, sqrt(vare))
	empirical_se <- summary(lm(y ~ g))$coefficients[2,2]
	ese <- expected_se(eff, mean(g)/2, nid, var(y))

	# 5% error maximum?
	expect_true(abs(empirical_se - ese)/ese < 0.05)

})



