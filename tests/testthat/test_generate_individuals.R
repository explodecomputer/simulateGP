context("generate individuals")
library(simulateGP)

test_that("make_geno", {

	g <- make_geno(10000, 80, 0.5)

	expect_true(nrow(g) == 10000)
	expect_true(ncol(g) == 80)
})

test_that("make_phen", {

	set.seed(1234)
	g <- make_geno(10000, 50, 0.5)
	effs <- choose_effects(50, 0.05)
	expect_equal(sum(effs^2), 0.05)

	x <- make_phen(effs, g)
	prs <- g %*% effs
	cor(prs, x)^2
	expect_equal(cor(prs, x)[1,1]^2, sum(effs^2), tolerance=0.015)
})


test_that("y_to_binary", {

	set.seed(1234)
	g <- make_geno(10000, 50, 0.5)
	effs <- choose_effects(50, 0.05)
	expect_equal(sum(effs^2), 0.05)

	x <- make_phen(effs, g)
	y <- y_to_binary(x, prevalence=0.01)	

	expect_equal(sum(y)/length(y), 0.01)

	a <- ascertain_samples(y, 0.5)
	expect_equal(sum(y[a])/length(a), 0.5)
})
