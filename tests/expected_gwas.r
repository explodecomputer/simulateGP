library(devtools)
load_all()
ss <- try(simulateGP::create_system(
    nidx=400000,
    nidy=400000,
    nidu=0,
    nu=3,
    na=0,
    nb=0,
    var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
    nsnp_x=80,
    nsnp_y=80,
    var_gx.x=runif(1, 0.01, 0.1),
    var_gy.y=runif(1, 0.01, 0.1),
    var_gx.y=runif(1, 0.001, 0.01),
    mu_gx.y=runif(1, -0.005, 0.005),
    prop_gx.y=runif(1, 0, 1),
    var_gy.x=runif(1, 0.001, 0.01),
    mu_gy.x=runif(1, -0.005, 0.005),
    prop_gy.x=runif(1, 0, 1)
))


o <- test_system(ss)




theoretical approach to simulating ss


g -> x
g -> y
g -> u -> x
g -> u -> y
g -> x -> y


a <- init_parameters(
var_x.y=sample(c(0, runif(5, 0.001, 0.1)), 1),
nsnp_x=80,
nsnp_y=80,
var_gx.x=runif(1, 0.01, 0.1),
var_gy.y=runif(1, 0.01, 0.1),
var_gx.y=runif(1, 0.001, 0.01),
mu_gx.y=runif(1, -0.005, 0.005),
prop_gx.y=runif(1, 0, 1),
var_gy.x=runif(1, 0.001, 0.01),
mu_gy.x=runif(1, -0.005, 0.005),
prop_gy.x=runif(1, 0, 1)
)
a <- sample_system_effects(a)

u <- rnorm(1000)
x <- rnorm(1000) + u*2
y <- rnorm(1000) + u + x*4

lm(x ~ u)
lm(y ~ u)

expected_gwas(9, 10000, 5, 99)

expected_gwas(0, 10000, 1, 2)

x <- rnorm(10000)
y <- rnorm(10000) + x
summary(lm(y ~ x))$coeff[2,3]^2

(cor(x,y)^2 * 10000-2) / (1-cor(x,y)^2)
F = 

pop <- simulate_population(a, 500000)


# direct effects
emp <- estimate_system_effects(pop)
qw


expx <- expected_gwas(a$eff_gx.x/sqrt(0.5), 500000, 1, 1)
plot(x$x$bhat[x$x$inst == 'x'], expx$bhat)
abline(0,1)

expyx <- expected_gwas(a$eff_gy.x/sqrt(0.5), 500000, 1, 1)
plot(x$x$bhat[x$x$inst == 'y'], expyx$bhat)
abline(0,1)

expxy <- expected_gwas(a$eff_gx.x/sqrt(0.5) * a$eff_x.y + a$eff_gx.y/sqrt(0.5), 500000, 1, 1)
plot(x$y$bhat[x$y$inst == 'x'], expxy$bhat)
abline(0,1)

expy <- expected_gwas(a$eff_gy.y/sqrt(0.5), 500000, 1, 1)
plot(x$y$bhat[x$y$inst == 'y'], expy$bhat)
abline(0,1)



expx <- expected_gwas(a$eff_gx.x/sqrt(0.5), 500000, 1, 1)
plot(emp$x$se[emp$x$inst == 'x'], expx$se)
abline(0,1)

expyx <- expected_gwas(a$eff_gy.x/sqrt(0.5), 500000, 1, 1)
plot(emp$x$se[emp$x$inst == 'y'], expyx$se)
abline(0,1)

expxy <- expected_gwas(a$eff_gx.x/sqrt(0.5) * a$eff_x.y + a$eff_gx.y/sqrt(0.5), 500000, 1, 1)
plot(emp$y$se[emp$y$inst == 'x'], expxy$se)
abline(0,1)

expy <- expected_gwas(a$eff_gy.y/sqrt(0.5), 500000, 1, 1)
plot(emp$y$se[emp$y$inst == 'y'], expy$se)
abline(0,1)

F = rsq / ()

alternative_expected_gwas <- function(eff, n, vx, vy)
{
	rsq <- eff^2 * vx / vy
	fval = (rsq * (n-2)) / (1 - rsq)
	tval = sqrt(fval)
	se = abs(eff / tval)
	p <- pt(abs(tval), n-1, lower.tail=FALSE)
	dat <- tibble::data_frame(bhat=eff, se=se, fval=tval^2, pval=p, n=n)
	return(dat)
}

# problem - not getting correct standard error
# fval seems to be only correct if half sample size is assumed??


