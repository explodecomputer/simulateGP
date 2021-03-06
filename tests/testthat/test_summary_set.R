context("GWAS summary data simulations")
library(simulateGP)


test_that("no sample overlap", {

	set.seed(1234)
	nsnp <- 100
	map <- dplyr::tibble(
		af = runif(nsnp, 0.01, 0.99),
		snp = 1:nsnp
	)
	beta_gx <- generate_gwas_params(map, h2=0.5)
	beta_gy <- beta_gx
	beta_gy$beta <- beta_gy$beta * 0.3

	dat <- summary_set(
		beta_gx = beta_gx$beta,
		beta_gy = beta_gy$beta,
		af = beta_gx$af,
		n_gx = 10000,
		n_gy = 10000,
		n_overlap = 0,
		cor_xy = 0,
		prev_y=NA,
		sigma_x=1,
		sigma_y=1
	)

	res <- dat %>%
		dplyr::filter(pval.exposure < 5e-8) %>%
		TwoSampleMR::mr(., method_list="mr_ivw")

	expect_equal(0.3, res$b, tolerance=res$se * 1.96)
})



test_that("complete sample overlap linear linear", {

	set.seed(1234)
	nsnp <- 100
	map <- dplyr::tibble(
		af = runif(nsnp, 0.01, 0.99),
		snp = 1:nsnp
	)	
	beta_gx <- generate_gwas_params(map, h2=0.5)
	beta_gy <- beta_gx
	beta_gy$beta <- beta_gy$beta * 0.3

	dat <- summary_set(
		beta_gx = beta_gx$beta,
		beta_gy = beta_gy$beta,
		af = beta_gx$af,
		n_gx = 10000,
		n_gy = 10000,
		n_overlap = 10000,
		cor_xy = 0.6,
		prev_y=NA,
		sigma_x=1,
		sigma_y=1
	)

	res <- dat %>%
		dplyr::filter(pval.exposure < 5e-8) %>%
		TwoSampleMR::mr(., method_list="mr_ivw")

	expect_equal(0.3, res$b, tolerance=res$se * 1.96)
})



test_that("partial sample overlap linear linear", {

	set.seed(123456)
	nsnp <- 100
	map <- dplyr::tibble(
		af = runif(nsnp, 0.01, 0.99),
		snp = 1:nsnp
	)
	beta_gx <- generate_gwas_params(map, h2=0.7)
	beta_gy <- beta_gx
	beta_gy$beta <- beta_gy$beta * 0.3

	dat <- summary_set(
		beta_gx = beta_gx$beta,
		beta_gy = beta_gy$beta,
		af = beta_gx$af,
		n_gx = 10000,
		n_gy = 10000,
		n_overlap = 5000,
		cor_xy = 0.1,
		prev_y=NA,
		sigma_x=1,
		sigma_y=1
	)

	res <- dat %>%
		dplyr::filter(pval.exposure < 5e-8) %>%
		TwoSampleMR::mr(., method_list="mr_ivw")

	expect_equal(0.3, res$b, tolerance=res$se * 1.96)
})


test_that("partial sample overlap linear logistic", {

	set.seed(1234567)
	nsnp <- 100
	map <- dplyr::tibble(
		af = runif(nsnp, 0.01, 0.99),
		snp = 1:nsnp
	)	
	beta_gx <- generate_gwas_params(map, h2=0.5)
	beta_gy <- beta_gx
	beta_gy$beta <- beta_gy$beta * 0.3

	dat <- summary_set(
		beta_gx = beta_gx$beta,
		beta_gy = beta_gy$beta,
		af = beta_gx$af,
		n_gx = 10000,
		n_gy = 10000,
		n_overlap = 5000,
		cor_xy = 0.6,
		prev_y=0.1,
		sigma_x=1,
		sigma_y=1
	)

	res <- dat %>%
		dplyr::filter(pval.exposure < 5e-8) %>%
		TwoSampleMR::mr(., method_list="mr_ivw")

	expect_equal(0.3, res$b, tolerance=res$se * 1.96)
})


