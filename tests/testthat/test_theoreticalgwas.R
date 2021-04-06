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


test_that("gwas summary data", {
	set.seed(1345)
	nsnp <- 1000
	ldobjlist <- test_ldobj(nsnp, 100)
	map <- lapply(ldobjlist, function(x) x$map) %>%
		dplyr::bind_rows()
	out <- map %>%
		generate_gwas_params(h2=0.3, S=0.3, Pi=3/nsnp) %>%
		generate_gwas_ss(100000, ldobjlist=ldobjlist)
	expect_equal(nrow(out), nsnp)
	expect_true(sum(out$pval < 5e-8) > 3)
	# ggplot(out, aes(pos, -log10(pval))) + geom_point() + facet_grid(. ~ chr)

	ldobj <- ldobjlist[[1]]
	nsnp <- nrow(ldobj$map)
	out <- ldobj$map %>%
		generate_gwas_params(h2=0.3, S=0.3, Pi=3/nsnp) %>%
		generate_gwas_ss(100000, ldobj=ldobj)
	expect_equal(nrow(out), nsnp)
	expect_true(sum(out$pval < 5e-8) > 3)

	out <- ldobj$map %>%
		generate_gwas_params(h2=0.3, S=0.3, Pi=3/nsnp) %>%
		generate_gwas_ss(100000, ld=ldobj$ld)
	expect_equal(nrow(out), nsnp)
	expect_true(sum(out$pval < 5e-8) > 3)


})

