context("mr system")
library(simulateGP)

test_that("create system", {

	dat <- create_system(nidx=1000, nidy=1000, nidu=0, nu=0, na=0, nb=0, var_x.y=0.1, nsnp_x=100, var_gx.x=0.1, var_gx.y=0, mu_gx.y=0, prop_gx.y=1, nsnp_y=0, var_gy.y=0, var_gy.x=0, mu_gy.x=0, prop_gy.x=1)
	expect_true(is.list(dat))
})

test_that("init_parameters", {

	dat <- init_parameters(var_x.y=0.1, nsnp_x=100, var_gx.x=0.1, var_gx.y=0, mu_gx.y=0, prop_gx.y=1, nsnp_y=0, var_gy.y=0, var_gy.x=0, mu_gy.x=0, prop_gy.x=1)
	expect_true(is.list(dat))

	dat <- add_u(dat, nsnp_u=100, var_u.x=0.1, var_u.y=0.1, var_gu.u=0.1)
	expect_true(is.list(dat))

	dat2 <- sample_system_effects(dat)
	expect_true(is.list(dat2))

	dat3 <- simulate_population(dat2, 1000)
	expect_true(is.list(dat3))

	dat4 <- estimate_system_effects(dat3)
	expect_true(is.list(dat4))
})


test_that("test_system", {

	skip("need random forest")
	ss <- create_system(nidx=10000, nidy=10000, nidu=0, nu=0, na=0, nb=0, var_x.y=0.1, nsnp_x=10, var_gx.x=0.1, var_gx.y=0, mu_gx.y=0, prop_gx.y=1, nsnp_y=10, var_gy.y=0.1, var_gy.x=0, mu_gy.x=0, prop_gy.x=1)
	res <- test_system(ss)
	expect_true(is.list(res))
})
