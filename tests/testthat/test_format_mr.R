context("format MR")
library(simulateGP)


test_that("mvmr", {

	# Simulate 100 genotypes
	g <- make_geno(10000, 80, 0.5)

	# Choose effect sizes for instruments for each trait
	effs1 <- choose_effects(50, 0.05)
	effs2 <- choose_effects(50, 0.05)

	# Create X1 and X2, where they overlap some variants
	x1 <- make_phen(effs1, g[,1:50])
	x2 <- make_phen(effs2, g[,31:80])

	# Create Y - x1 has a -0.3 influence on it and x2 has a +0.3 influence on it
	y <- make_phen(c(-0.3, 0.3), cbind(x1, x2))

	# Perform separate MR on each
	dat1 <- get_effs(x1, y, g)
	dat2 <- get_effs(x2, y, g)

	# Do multivariable MR
	# First get the effects for x1, x2 and y, and put them in mv format
	mvdat <- make_mvdat(list(x1, x2), y, g)

	expect_true(is.list(mvdat))

})


test_that("merge_exp_out", {

	# Simulate 100 genotypes
	g <- make_geno(10000, 80, 0.5)

	# Choose effect sizes for instruments for each trait
	effs <- choose_effects(80, 0.05)

	# Create X1 and X2, where they overlap some variants
	x <- make_phen(effs, g)

	# Create Y - x1 has a -0.3 influence on it and x2 has a +0.3 influence on it
	y <- make_phen(-0.3, x)

	gx <- gwas(x, g)
	gy <- gwas(y, g)

	dat <- merge_exp_out(gx, gy)
	expect_true(nrow(dat) == 80)

	dat2 <- recode_dat(dat)
	expect_true(nrow(dat2) == 80)
})

