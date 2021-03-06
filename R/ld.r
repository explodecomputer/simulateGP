#' Create test LD object
#'
#' @param nsnp Number of SNPs
#' @param chunksize Chunksize for splitting
#'
#' @export
#' @return list of chunks, which each contain map and LD matrix
test_ldobj <- function(nsnp, chunksize)
{
	snp <- 1:nsnp
	nchunk <- ceiling(nsnp/chunksize)
	start <- 0:(nchunk-1) * chunksize + 1
	end <- pmin(1:nchunk * chunksize, nsnp)

	ldobj <- lapply(1:nchunk, function(i){
		n <- length(start[i]:end[i])
		p <- qr.Q(qr(matrix(rnorm(n^2), n)))
		Sigma <- crossprod(p, p*(5:1))
		denom <- sqrt(diag(Sigma)) %*% t(sqrt(diag(Sigma)))
		rho <- Sigma / denom
		map <- dplyr::tibble(
			snp=start[i]:end[i],
			chr=i,
			pos=1:n,
			ref="A",
			alt="C"
		)
		rownames(rho) <- colnames(rho) <- map$snp
		return(list(map=map, ld=rho))
	})
	return(ldobj)
}


#' Get LD matrix for a specified region from bfile reference panel
#'
#' @param chr Chromosome
#' @param from from bp
#' @param to to bp
#' @param bfile LD reference panel
#' @param snplist list of rsids to retain. Default is NULL (i.e. retain all)
#' @param idlist list of sample ids to retain. Default is NULL (i.e. retain all)
#' @param plink_bin Plink binary default=genetics.binaRies::get_plink_binary()
#'
#' @export
#' @return List of LD matrix and map info including MAF
get_ld <- function(chr, from, to, bfile, plink_bin=genetics.binaRies::get_plink_binary())
{
	# Make textfile
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	fn <- tempfile()

	fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --chr ", chr,
		" --from-bp ", from, 
		" --to-bp ", to,
		" --r square ", 
		" --keep-allele-order ",
		" --make-just-bim ",
		" --freq ",
		" --out ", shQuote(fn, type=shell)
	)

	stat <- system(fun1, ignore.stdout=TRUE)

	if(stat==0)
	{
		x <- data.table::fread(paste0(fn, ".ld")) %>% as.matrix()
		y <- data.table::fread(paste0(fn, ".bim")) %>% dplyr::as_tibble()
		z <- data.table::fread(paste0(fn, ".frq")) %>% dplyr::as_tibble()
		names(y) <- c("chr", "snp", "gp", "pos", "alt", "ref")
		y <- dplyr::select(y, -gp)
		y$af <- z$MAF
		out <- list(ld=x, map=y)
	} else {
		out <- NULL
	}
	unlink(paste0(fn, c(".ld", ".bim", ".frq")))
	return(out)
}
# Also see https://github.com/explodecomputer/pic_haps/


#' <brief desc>
#'
#' <full description>
#'
#' @param bfile <what param does>
#' @param regionfile <what param does>
#' @param plink_bin=genetics.binaRies::get_plink_binary() <what param does>
#' @param nthreads=1 <what param does>
#'
#' @export
#' @return
generate_ldobj <- function(outdir, bfile, regions, plink_bin=genetics.binaRies::get_plink_binary(), nthreads=1)
{
	codes <- paste0(regions$chr, "_", regions$start, "_", regions$stop)
	map <- parallel::mclapply(1:nrow(regions), function(i)
	{
		message(i, " of ", nrow(regions))
		out <- get_ld(
			chr=regions$chr[i] %>% gsub("chr", "", .),
			from=regions$start[i],
			to=regions$stop[i],
			bfile=bfile,
			plink_bin=plink_bin
		)
		if(!is.null(out))
		{
			fn <- file.path(outdir, paste0("ldobj_", codes[i], ".rds"))
			saveRDS(out, file=fn)
			return(out$map)
		} else {
			return(NULL)
		}
	}, mc.cores=nthreads) %>%
		dplyr::bind_rows()
	saveRDS(map, file.path(outdir, "map.rds"))
	return(map)
}


#' Simulate two correlated binomial variables
#'
#' @param nid Number of samples
#' @param p1 Frequency 1
#' @param p2 Frequency 2
#' @param rho Target correlation
#' @param n Binomial parameter, should be 2 (default) for genotypes
#' @param round Round or not Default=TRUE
#' @param print Print or not Default=FALSE
#'
#' @export
#' @return Matrix
correlated_binomial <- function(nid, p1, p2, rho, n=2, round=TRUE, print=FALSE)
{
	# from https://stats.stackexchange.com/questions/284996/generating-correlated-binomial-random-variables
	p <- p1
	q <- p2
	a <- function(rho, p, q) {
		rho * sqrt(p*q*(1-p)*(1-q)) + (1-p)*(1-q)
	}

	a.0 <- a(rho, p, q)
	prob <- c(`(0,0)`=a.0, `(1,0)`=1-q-a.0, `(0,1)`=1-p-a.0, `(1,1)`=a.0+p+q-1)
	if (min(prob) < 0) {
		print(prob)
		stop("Error: a probability is negative.")
	}
	#
	# Illustrate generation of correlated Binomial variables.
	#
	n.sim <- nid
	u <- sample.int(4, n.sim * n, replace=TRUE, prob=prob)
	y <- floor((u-1)/2)
	x <- 1 - u %% 2
	x <- colSums(matrix(x, nrow=n)) # Sum in groups of `n`
	y <- colSums(matrix(y, nrow=n)) # Sum in groups of `n`

	if(round)
	{
		x <- round(x)
		y <- round(y)
	}

	if(print)
	{
		print(table(x, y))
		print(stats::cor(x, y))
	}
	return(cbind(x, y))

	#
	# Plot the empirical bivariate distribution.
	#
	# plot(x+rnorm(length(x), sd=1/8), y+rnorm(length(y), sd=1/8),
	#      pch=19, cex=1/2, col="#00000010",
	#      xlab="X", ylab="Y",
	#      main=paste("Correlation is", signif(cor(x,y), 3)))
	# abline(v=mean(x), h=mean(y), col="Red")
	# abline(lm(y ~ x), lwd=2, lty=3)
}

#' Estimate haplotype frequencies for two loci
#'
#' @param r Required LD r
#' @param p1 Freq 1
#' @param p2 Freq 2
#'
#' @export
#' @return vector
hap_freqs <- function(r, p1, p2)
{
	# d = pAB - p1p2
	# d = p1q2 - pAb
	# d = q1p2 - paB
	# d = pab - q1q2
	# r = d / denom
	denom <- sqrt(p1 * (1-p1) * p2 * (1-p2))
	p <- c(
		`(1,1)` = p1 * p2 + r * denom, 
		`(0,0)` = (1-p1) * (1-p2) + r * denom, 
		`(1,0)` = p1 * (1-p2) - r * denom, 
		`(0,1)` = (1-p1) * p2 - r * denom
	)
	stopifnot(all(p >= 0))
	return(p)
}

test_hap_freqs <- function(r, p1, p2)
{
	a <- try(hap_freqs(r, p1, p2), silent=TRUE)
	return(ifelse(class(a)=='try-error', FALSE, TRUE))
}

#' Simulate haplotypes of two loci
#'
#' @param nid Number of samples
#' @param r Desired LD r
#' @param p1 Freq 1
#' @param p2 Freq 2
#'
#' @export
#' @return Matrix
simulate_haplotypes <- function(nid, r, p1, p2)
{
	p <- hap_freqs(r, p1, p2)
	n <- round(p * nid)
	diff <- sum(n) - nid
	i <- 1
	while(diff != 0)
	{
		if(diff > 0)
		{
			n[i] <- n[i] - 1
		} else {
			n[i] <- n[i] + 1
		}
		diff <- sum(n) - nid
		i <- ifelse(i == 4, 1, i + 1)
	}
	mat <- rbind(
		cbind(rep(1, n[1]), rep(1, n[1])),
		cbind(rep(0, n[2]), rep(0, n[2])),
		cbind(rep(1, n[3]), rep(0, n[3])),
		cbind(rep(0, n[4]), rep(1, n[4]))
	)
	ind1 <- sample(1:nid, nid, replace=FALSE)
	ind2 <- sample(1:nid, nid, replace=FALSE)
	haps <- cbind(
		A1 = mat[ind1, 1],
		A2 = mat[ind2, 1],
		B1 = mat[ind1, 2],
		B2 = mat[ind2, 2]
	)
	return(haps)
}

#' Simulate genotypes from haplotypes
#'
#' @param nid Number of samples
#' @param r Desired LD r
#' @param p1 Freq 1
#' @param p2 Freq 2
#'
#' @export
#' @return Matrix
simulate_geno <- function(nid, r, p1, p2)
{
	a <- simulate_haplotypes(nid, r, p1, p2)
	b <- cbind(a[,1]+a[,2], a[,3]+a[,4])
	return(b)
}

