library(devtools)
load_all()
n <- 10000
maf <- 0.3
g <- make_geno(n, 1, maf)
table(g)
e <- rnorm(n)
beta <- 2
y <- g * beta + e

summary(lm(y ~ g))

sigma <- 1 / (n-2) * sum(e^2)
denom <- sum((g - mean(g))^2)
sqrt(sigma / denom)

vy <- var(y)
vg <- var(g)

r2 <- 2^2 * 2 * 0.3 * 0.7 / vy

r2 <- vg / vy = vg / (vg + ve)
r2*vg + r2*ve = vg
ve = (vg - r2*vg)/r2

ve = vg*(1-r2)/r2

vy

tsigma <- vy - vy * r2



y <- make_phen(sqrt(0.4), g)


nsnp <- 10000
nid <- 10000
maf <- runif(nsnp)
g <- make_geno(nid, nsnp, maf)
eff <- rnorm(nsnp, sd=sqrt(0.2))
y <- g %*% eff + rnorm(nid, sd=sqrt(0.8))
# out <- gwas(y, g)
# summary(lm(y ~ g))
sum(eff^2 * 2 * maf * (1-maf))
mean(abs(eff))^2 * 2 * mean(maf) * (1-mean(maf)) * length(maf)

sqrt(0.2) * sqrt(2)/sqrt(pi)





library(simulateGP)

beta <- c(100)
hist(beta, breaks=100)
maf <- rbeta(length(beta), 0.8, 1)/2
hist(maf)
dat <- theoretical_gwas(beta, maf, 0.2, c(300000, 150000), minmaf=0.01)

library(ggplot2)

o <- subset(dat[[1]], pval < 5e-8)
ggplot(o, aes(x=abs(beta), y=abs(betahat))) +
geom_point() +
geom_smooth(method=lm) +
geom_abline(slope=1,intercept=0)

summary(lm(abs(betahat) ~ abs(beta), o))


