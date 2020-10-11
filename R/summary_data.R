variant_reference <- function(bfile, plink_bin=genetics.binaRies::get_plink_binary(), fn=tempfile())
{
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")

	fun1 <- paste0(
		shQuote(plink_bin, type=shell),
		" --bfile ", shQuote(bfile, type=shell),
		" --freq ",
		" --out ", shQuote(fn, type=shell)
	)
	message("Generating MAF")
	system(fun1, ignore.stdout=TRUE)

	message("Reading variant info")
	frq <- data.table::fread(paste0(fn, ".frq"))
	bim <- data.table::fread(paste0(bfile, ".bim"))

	names(bim) <- c("CHR", "RSID", "GP", "POS", "EA", "NEA")
	bim$EAF <- frq$MAF
	bim <- subset(bim, !is.na(EAF))
	return(bim)
}


get_regions <- function(pop="ASN")
{
	system.file(paste0("extdata/ldetect/", pop, ".bed"), package="simulateGP") %>%
		data.table::fread(., header=TRUE) %>%
		dplyr::mutate(
			chr=as.numeric(gsub("chr", "", chr)),
			start=as.numeric(start),
			stop=as.numeric(stop)
		) %>% 
		dplyr::as_tibble()
}


generate_ld_matrices <- function(regions, varref, bfile, plink_bin=genetics.binaRies::get_plink_binary(), fn=tempfile())
{
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	message("Calculating LD for ", nrow(regions), " regions")
	l <- list()
	for(i in 1:nrow(regions))
	{
		message("Region ", i, " of ", nrow(regions))
		variants <- subset(varref, CHR == regions$chr[i] & POS > regions$start[i] & POS < regions$stop[i])
		l[[i]] <- list(
			info=variants,
			ld=ieugwasr::ld_matrix(variants$RSID, bfile=bfile, plink_bin=plink_bin, with_alleles=FALSE)
		)
	}
}

# LD from 
generate_ld_matrices_slow <- function(regions, varref, bfile, plink_bin=genetics.binaRies::get_plink_binary(), fn=tempfile())
{
	shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
	message("Calculating LD for ", nrow(regions), " regions")
	l <- list()
	for(i in 1:nrow(regions))
	{
		message("Region ", i, " of ", nrow(regions))
		variants <- subset(varref, CHR == regions$chr[i] & POS > regions$start[i] & POS < regions$stop[i])
		write.table(data.frame(variants$RSID), file=fn, row.names=F, col.names=F, quote=F)
		fun2 <- paste0(
			shQuote(plink_bin, type=shell),
			" --bfile ", shQuote(bfile, type=shell),
			" --extract ", shQuote(fn, type=shell), 
			" --recode A ", 
			" --out ", shQuote(fn, type=shell)
		)
		system(fun2, ignore.stdout=TRUE)
		x <- data.table::fread(paste0(fn, ".raw")) %>% {.[,-c(1:6)]} %>% as.matrix()
		l[[i]] <- list(
			info=variants,
			ld=cor(x, use="pair")
		)
		unlink(paste0(fn, ".raw"))
	}
}


greedy_remove <- function(r, threshold=0.99)
{
	diag(r) <- 0
	flag <- 1
	rem <- c()
	nom <- colnames(r)
	while(flag == 1)
	{
		message("iteration")
		count <- apply(r, 2, function(x) sum(x >= threshold))
		if(any(count > 0))
		{
			worst <- which.max(count)[1]
			rem <- c(rem, names(worst))
			r <- r[-worst,-worst]
		} else {
			flag <- 0
		}
	}
	return(which(nom %in% rem))
}


ld_multiplier <- function(varref, ld)
{
	varref$beta_rho <- NA
	varref$b_rho <- NA
	for(i in 1:length(ld))
	{
		message("Region ", i, " of ", length(ld))
		m1 <- match(ld[[i]]$info$RSID, varref$RSID) %>% na.exclude %>% as.numeric
		variants_dat <- varref$RSID[m1]
		if(length(variants_dat) > 0)
		{
			m2 <- match(variants_dat, rownames(ld[[i]]$ld))
			r <- ld[[i]]$ld[m2,m2]
			varref$beta_rho[m1] <- varref$beta[m1] %*% r
			varref$b_rho[m1] <- varref$b[m1] %*% r
		}
	}
	varref$beta_rho_se <- expected_se(varref$beta_rho, varref$EAF, varref$N, varref$vy)
	varref$b_rho_se <- expected_se(varref$b_rho, varref$EAF, varref$N, varref$vy)
	return(varref)
}


draw_betas_multi_sample <- function(b, se, N, pcor, Nrep=1)
{
	stopifnot(length(b) == length(se))
	stopifnot(length(b) == nrow(N))
	stopifnot(nrow(N) == nrow(N))
	stopifnot(all(dim(N) == dim(pcor)))
	stopifnot(all(diag(N) == 1))
	stopifnot(all(diag(pcor) == 1))
	ses <- diag(se)
	covmat <- ses %*% (pcor * N) %*% ses
	mvrnorm(Nrep, b, covmat)
}


calculate_overlap_2 <- function(n1, n2, n_overlap)
{
	n_overlap / sqrt(n1 * n2)
}


calculate_overlap <- function(n, N_overlap)
{
	stopifnot(length(n) == nrow(N_overlap))
	stopifnot(nrow(N_overlap) == ncol(N_overlap))
	N_overlap / sqrt(n %*% t(n))
}


