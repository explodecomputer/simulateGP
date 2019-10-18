library(devtools)
load_all()
library(TwoSampleMR)


ss <- try(create_system(
	nidx=sample(20000:500000, 1),
	nidy=sample(20000:500000, 1),
	nidu=0,
	nu=sample(0:10, 1),
	na=0,
	nb=0,
	var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
	nsnp_x=sample(1:200, 1),
	nsnp_y=sample(1:200, 1),
	var_gx.x=runif(1, 0.01, 0.1),
	var_gy.y=runif(1, 0.01, 0.1),
	var_gx.y=runif(1, 0.001, 0.01),
	mu_gx.y=runif(1, -0.005, 0.005),
	prop_gx.y=runif(1, 0, 1),
	var_gy.x=runif(1, 0.001, 0.01),
	mu_gy.x=runif(1, -0.005, 0.005),
	prop_gy.x=runif(1, 0, 1)
))





nop <- try(create_system(
	nidx=100000,
	nidy=100000,
	nidu=0,
	nu=0,
	na=0,
	nb=0,
	var_x.y=0,
	nsnp_x=100,
	nsnp_y=100,
	var_gx.x=0.15,
	var_gy.y=0.15,
	var_gx.y=0,
	mu_gx.y=0,
	prop_gx.y=0,
	var_gy.x=0,
	mu_gy.x=0,
	prop_gy.x=0
))




withu <- try(create_system(
	nidx=10000,
	nidy=10000,
	nidu=10000,
	nu=2,
	na=2,
	nb=2,
	var_x.y=0,
	nsnp_x=100,
	nsnp_y=100,
	var_gx.x=0.15,
	var_gy.y=0.15,
	var_gx.y=0,
	mu_gx.y=0,
	prop_gx.y=0,
	var_gy.x=0,
	mu_gy.x=0,
	prop_gy.x=0
))

a <- make_dat(withu$x$x, withu$y$y) %>% filter(pval.exposure < 1e-8)
b <- make_dat(nop$x$x, nop$y$y) %>% filter(grepl("x", SNP), pval.exposure < 1e-8)

mr_scatter_plot(mr(a), a)
mr_scatter_plot(mr(b), b)
mr_heterogeneity(a)
mr_heterogeneity(b)



nop1 <- try(create_system(
	nidx=10000,
	nidy=10000,
	nidu=0,
	nu=0,
	na=0,
	nb=0,
	var_x.y=0.1,
	nsnp_x=100,
	nsnp_y=100,
	var_gx.x=0.15,
	var_gy.y=0.15,
	var_gx.y=0,
	mu_gx.y=0,
	prop_gx.y=0,
	var_gy.x=0,
	mu_gy.x=0,
	prop_gy.x=0
))




withu1 <- try(create_system(
	nidx=10000,
	nidy=10000,
	nidu=10000,
	nu=2,
	na=2,
	nb=2,
	var_x.y=0.1,
	nsnp_x=100,
	nsnp_y=100,
	var_gx.x=0.15,
	var_gy.y=0.15,
	var_gx.y=0,
	mu_gx.y=0,
	prop_gx.y=0,
	var_gy.x=0,
	mu_gy.x=0,
	prop_gy.x=0
))


a <- make_dat(withu1$x$x, withu1$y$y) %>% filter(pval.exposure < 1e-8)
b <- make_dat(nop1$x$x, nop1$y$y) %>% filter(grepl("x", SNP), pval.exposure < 1e-8)

mr_scatter_plot(mr(a), a)
mr_scatter_plot(mr(b), b)
mr_heterogeneity(a)
mr_heterogeneity(b)

