load_all()


va <- 0.1
maf <- 0.5
prevalence <- 0.1
nid <- 100000
proportion <- 0.5


g <- make_geno(nid, 1, maf)
x <- make_phen(sqrt(va), g)
bin <- y_to_binary(x, prevalence=prevalence)








asc <- ascertain_samples(bin, proportion)
ga <- g[asc, ]
xa <- x[asc]
bina <- bin[asc]

a <- fast_assoc(xa, ga)
mod <- glm(bina ~ ga, family="binomial")

freq <- sum(ga) / (2 * length(ga))

rsq <- mod$coefficients[2]^2 * freq * (1 - freq) * 2 / var_cc(0.1, 0.9)
rsq


	h2l


mod <- lm(xa ~ ga)
rsq <- mod$coefficients[2]^2 * freq * (1 - freq) * 2 / var_cc(0.8, 0.9)
rsq



summary(lm(x ~ g))
summary(glm(bin ~ g, family="binomial"))
var_cc(prevalence, prevalence)


res <- array(0, 50)
proportion <- seq(0.1, 0.9, length.out=50)
for(i in 1:50)
{
	message(i)
	va <- 0.1
	maf <- 0.5
	prevalence <- 0.1
	nid <- 100000

	g <- make_geno(nid, 1, maf)
	x <- make_phen(sqrt(va), g)

	d <- risk_simulation(g, eff=0.1, prevalence=0.1, prop_discovered=1)
	bin <- d$disease

	asc <- ascertain_samples(bin, proportion[i])
	ga <- g[asc, ]
	bina <- bin[asc]
	mod <- glm(bina ~ ga, family="binomial")
	res[i] <- mod$coefficients[2]
}


plot(res ~ proportion)
summary(glm(bin ~ g, family="binomial"))




summary(glm(bin ~ g, family="binomial"))



0.9 * 0.1
table(bin)
var(bin)





var_cc <- function(proportion, prevalence)
{
	tval <- qnorm(prevalence, lower.tail=FALSE)
	zval <- dnorm(tval)
	i <- zval / prevalence
	lambda <- (proportion - prevalence) / (1 - prevalence)
	vp <- 1 + i * lambda * (tval - i * lambda)
	return(vp)
}





contingency <- function(maf, prop, odds_ratio, eps=1e-15) {
  a <- odds_ratio-1
  b <- (maf+prop)*(1-odds_ratio)-1
  c_ <- odds_ratio*maf*prop

  if (abs(a) < eps) {
    z <- -c_ / b
  } else {
    d <- b^2 - 4*a*c_
    if (d < eps*eps) s <- 0 else s <- c(-1,1)
    z <- (-b + s*sqrt(max(0, d))) / (2*a)
  }
  y <- vapply(z, function(a) zapsmall(matrix(c(a, prop-a, maf-a, 1+a-maf-prop), 2, 2)),
              matrix(0.0, 2, 2))
  i <- apply(y, 3, function(u) all(u >= 0))
  return(y[,,i])
}





