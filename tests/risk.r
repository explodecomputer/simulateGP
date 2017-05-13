library(pROC)
library(ggplot2)
library(dplyr)
library(tidyr)

# h2l <- h2o * k^2 * (1-k)^2 / (p * (1 - p) * dnorm(threshold)^2)

# Calculating disease risk on observed scale
# Requires some 
# Check that gx to go works

d <- expand.grid(
	g = seq(-4,4,by=0.1),
	h2x = c(0.35, 0.5, 1),
	prev = c(0.004, 0.2, 0.5)
)
d$gr <- 1:nrow(d)
d$x_prime <- qnorm(d$prev, 0, 1)
d$e2x <- 1 - d$h2x
d$z <- dnorm(d$x_prime, mean=d$g, sd = sqrt(d$e2x))
d$p <- pnorm(d$x_prime, mean=d$g, sd = sqrt(d$e2x), lower.tail=FALSE)
d$p1 <- gx_to_gp(d$g, d$h2x, d$prev)

d <- group_by(d, h2x, prev) %>%
	mutate(p = range01(p), g = range01(g)) 


d1 <- gather(d, key, value, g, p)


ggplot(d1, aes(x=value, y=key)) +
geom_point(aes(colour=key)) +
geom_line(aes(group=gr)) +
facet_grid(h2x ~ prev)


# Model 1
# Common variant-common disease

nid <- 1000
nsnp <- 1000
h2x <- 0.3
prev <- 0.5
G_cdcv <- scale(make_geno(nid, nsnp, 0.5))
eff_cdcv <- rnorm(nsnp, sd=sqrt(h2x))
dat_cdcv <- risk_simulation(
	G=G_cdcv, 
	eff=eff_cdcv, 
	prevalence=prev,
	prop_discovered=0.1
)
plot(roc(disease ~ gx_true, dat_cdcv))



table(dat_cdcv$disease)
var(dat_cdcv$gx_true)
plot(roc(disease ~ gx_true, dat_cdcv))
plot(roc(disease ~ gx_pred, dat_cdcv))
risk_cross_plot(o=dat_cdcv$prob_disease, x=dat_cdcv$gx_true)
risk_cross_plot(o=dat_cdcv$disease, x=dat_cdcv$gx_true, title="Genetic values mapped to disease")
risk_cross_plot(o=dat_cdcv$disease, x=dat_cdcv$gx_pred, title="Genetic predictor of disease")


nid <- 100
nsnp <- 1000
h2x <- 0.8
prev <- 0.5
G_cdcv <- scale(make_geno(nid, nsnp, 0.5))
eff_cdcv <- rnorm(nsnp, sd=sqrt(h2x))
dat_cdcv <- risk_simulation(
	G=G_cdcv, 
	eff=eff_cdcv, 
	prevalence=prev,
	prop_discovered=0.1
)
pdf(file="roc_0.8_0.5.pdf")
plot(roc(disease ~ gx_true, dat_cdcv))
dev.off()

pdf(file="roc_0.8_0.5_0.1.pdf")
plot(roc(disease ~ gx_pred, dat_cdcv))
dev.off()

risk_cross_plot(o=dat_cdcv$disease, x=dat_cdcv$gx_true, title="Genetic values mapped to disease")
ggsave(file="cp_0.8_0.5.pdf", width=8,height=3)
risk_cross_plot(o=dat_cdcv$disease, x=dat_cdcv$gx_pred, title="Genetic predictor of disease using only 10% of causal variants")
ggsave(file="cp_0.8_0.5_0.1.pdf", width=8,height=3)



nid <- 1000
nsnp <- 1000
h2x <- 0.3
prev <- 0.5
G_cdcv <- scale(make_geno(nid, nsnp, 0.5))
eff_cdcv <- rnorm(nsnp, sd=sqrt(h2x))
dat_cdcv <- risk_simulation(
	G=G_cdcv, 
	eff=eff_cdcv, 
	prevalence=prev,
	prop_discovered=0.1
)
pdf(file="roc_0.3_0.5.pdf")
plot(roc(disease ~ gx_true, dat_cdcv))
dev.off()

pdf(file="roc_0.3_0.5_0.1.pdf")
plot(roc(disease ~ gx_pred, dat_cdcv))
dev.off()

risk_cross_plot(o=dat_cdcv$disease, x=dat_cdcv$gx_pred, title="Genetic predictor of disease using only 10% of causal variants")


nid <- 10000
nsnp <- 1000
h2x <- 0.3
prev <- 0.01
G_cdcv <- scale(make_geno(nid, nsnp, 0.5))
eff_cdcv <- rnorm(nsnp, sd=sqrt(h2x))
dat_cdcv <- risk_simulation(
	G=G_cdcv, 
	eff=eff_cdcv, 
	prevalence=prev,
	prop_discovered=0.1
)
pdf(file="roc_0.3_0.01.pdf")
plot(roc(disease ~ gx_true, dat_cdcv))
dev.off()


# Model 2
# Every case has a specific mutation
# Rare variant-common disease

nid <- 1000
nsnp <- 1000
h2x <- 0.3
prev <- 0.5
G_cdrv <- diag(nid)
diag(G_cdrv)[1:(nid/2)] <- 0
eff_cdrv <- scale(rnorm(nid, diag(G_cdrv), 0.001)) * sqrt(h2x)
dat_cdrv <- risk_simulation(
	G=G_cdrv, 
	eff=eff_cdrv, 
	prevalence=prev,
	prop_discovered=0.1
)

table(dat_cdrv$disease)
var(dat_cdrv$gx_true)

pdf("roc_0.3_0.5_rare.pdf")
plot(roc(disease ~ gx_true, dat_cdrv))
dev.off()

pdf("roc_0.3_0.5_0.1_rare.pdf")
plot(roc(disease ~ gx_pred, dat_cdrv))
dev.off()

risk_cross_plot(o=dat_cdrv$prob_disease, x=dat_cdrv$gx_true)
risk_cross_plot(o=dat_cdrv$disease, x=dat_cdrv$gx_true, title="Genetic values mapped to disease")
risk_cross_plot(o=dat_cdrv$disease, x=dat_cdrv$gx_pred, title="Genetic predictor of disease")

