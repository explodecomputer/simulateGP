context("LD")
library(simulateGP)

test_that("get ld matrix", {

	skip("needs gwas files")
	chr <- 1
	from <- 1892607
	to <- 3582736
	bfile <- "/Users/gh13047/repo/mr-base-api/app/ld_files/EUR"
	plink_bin <- genetics.binaRies::get_plink_binary()
	a <- get_ld(chr, from, to, bfile, plink_bin)	
	expect_true(length(a) == 2)
	expect_true(nrow(a[[1]]) == nrow(a[[2]]))
	expect_true(nrow(a[[1]]) == ncol(a[[1]]))
})
