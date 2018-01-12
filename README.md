# Simulate stuff

A collection of instruments for simulating various things.

To install:

```
devtools::install_github("explodecomputer/simulateGP")
```

## MR

Here we want to achieve the following:

1. Simulate some genetic effects for the exposure
2. Simulate the effect of the exposure on the outcome
3. Get the SNP-exposure and SNP-outcome effects and standard errors for these simulations

### Simple example

Here is how to do 1 and 2:

```r
# Genotypes for 20 SNPs and 10000 individuals, where the SNPs have MAF = 0.5:
g <- make_geno(10000, 20, 0.5)

# These SNPs instrument some exposure, and together explain 5% of the variance
effs <- choose_effects(20, 0.05)

# Create X
x <- make_phen(effs, g)

# Check that the SNPs explain 5% of the variance in x
sum(cor(x, g)^2)

# Create Y, which is negatively influenced by x, explaining 10% of the variance in Y
y <- make_phen(-sqrt(0.1), x)
```


We now have an X and Y, and the genotypes. To perform 2-sample MR on this we need to obtain the summary data.

```r
# Get the summary data
dat <- get_effs(x, y, g)

# Perform MR
library(TwoSampleMR)
mr(dat)
```

The effect size is approx -0.31, as we simulated.


### More complex example

Let's simulate an instance where there are two traits that have different effects on the outcome, but they share some of the same instruments

```r
# Simulate 100 genotypes
g <- make_geno(100000, 80, 0.5)

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
library(TwoSampleMR)
mr(subset(dat1, pval.exposure < 5e-8))
mr(subset(dat2, pval.exposure < 5e-8))

# Do multivariable MR
# First get the effects for x1, x2 and y, and put them in mv format
mvdat <- make_mvdat(list(x1, x2), y, g)

# Perform MV MR
mv_multiple(mvdat)
```


