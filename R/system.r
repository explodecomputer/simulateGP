library(dplyr)


#' Choose initial parameters for direct effects on X and Y
#'
#' @param nsnp_x
#' @param var_gx.x
#' @param var_x.y
#' @param var_gx.y=0
#' @param nsnp_y=0
#' @param var_gy.y=0
#' @param mu_gx.y=0
#' @param var_gy.x=0
#' @param mu_gy.x=0
#' @param prop_gy.x=1
#' @param prop_gx.y=1
#'
#' @export
#' @return List of model parameters
init_parameters <- function(nsnp_x, var_gx.x, var_x.y, var_gx.y=0, nsnp_y=0, var_gy.y=0, mu_gx.y=0, var_gy.x=0, mu_gy.x=0, prop_gy.x=1, prop_gx.y=1)
{
	parameters <- list(
		# Causal effect
		var_x.y = var_x.y,

		# Direct effects on x
		nsnp_x = nsnp_x,
		var_gx.x = var_gx.x,
		var_gx.y = var_gx.y,
		mu_gx.y = mu_gx.y,
		prop_gx.y = prop_gx.y,

		nsnp_y = nsnp_y,
		var_gy.y = var_gy.y,
		var_gy.x = var_gy.x,
		mu_gy.x = mu_gy.x,
		prop_gy.x = prop_gy.x,
		u = list()
	)
	return(parameters)
}

#' Add confounder variables and their instruments
#'
#' @param parameters Output from \code{init_parameters}
#' @param nsnp_u
#' @param var_u.x
#' @param var_u.y
#' @param var_gu.u
#'
#' @export
#' @return List of model parameters
add_u <- function(parameters, nsnp_u, var_u.x, var_u.y, var_gu.u)
{
	nom <- "u"
	if(var_u.x == 0 & var_u.y != 0) nom <- "b"
	if(var_u.x != 0 & var_u.y == 0) nom <- "a"
	i <- length(parameters$u) + 1
	parameters$u[[i]] <- list(
		nsnp_u = nsnp_u,
		var_u.x = var_u.x,
		var_u.y = var_u.y,
		var_gu.u = var_gu.u,
		name_u = nom
	)
	return(parameters)
}

#' Sample the actual effects based on initial parameters
#'
#' @param parameters Output from \code{init_parameters} or \code{add_u}
#'
#' @export
#' @return List of effect parameters
sample_system_effects <- function(parameters)
{
	nu <- length(parameters$u)
	if(nu > 0)
	{
		for(i in 1:nu)
		{
			parameters$u[[i]]$eff_gu.u <- choose_effects(parameters$u[[i]]$nsnp_u, parameters$u[[i]]$var_gu.u)
			parameters$u[[i]]$eff_u.x <- choose_effects(1, parameters$u[[i]]$var_u.x)
			parameters$u[[i]]$eff_u.y <- choose_effects(1, parameters$u[[i]]$var_u.y)
		}
	}

	# Make genetic effects for x instruments
	parameters$eff_gx.x <- choose_effects(parameters$nsnp_x, parameters$var_gx.x)

	# Make genetic effects for y instruments
	parameters$eff_gy.y <- choose_effects(parameters$nsnp_y, parameters$var_gy.y)


	# Create pleiotropic effects for some proportion of the effects
	# X-Y
	nchoose <- round(parameters$nsnp_x * parameters$prop_gx.y)
	ind <- sample(1:parameters$nsnp_x, nchoose, replace=FALSE)
	parameters$eff_gx.y <- rep(0, parameters$nsnp_x)
	if(nchoose > 0)
	{
		parameters$eff_gx.y[ind] <- choose_effects(nchoose, parameters$var_gx.y, mua=parameters$mu_gx.y)
	}

	# Create pleiotropic effects for some proportion of the effects
	# Y-X
	nchoose <- round(parameters$nsnp_y * parameters$prop_gy.x)
	ind <- sample(1:parameters$nsnp_y, nchoose, replace=FALSE)
	parameters$eff_gy.x <- rep(0, parameters$nsnp_y)
	if(nchoose > 0)
	{
		parameters$eff_gy.x[ind] <- choose_effects(nchoose, parameters$var_gy.x, mua=parameters$mu_gy.x)
	}

	parameters$eff_u.x <- sapply(parameters$u, function(x) choose_effects(1, x$var_u.x))
	parameters$eff_u.y <- sapply(parameters$u, function(x) choose_effects(1, x$var_u.y))
	parameters$eff_x.y <- sqrt(parameters$var_x.y)

	return(parameters)
}


#' Simulate individual level data from initial parameters
#'
#' @param parameters Output from \code{init_parameters} or \code{add_u}
#' @param nid Sample size to generate
#'
#' @export
#' @return List of matrices and vectors that represent individual level data
simulate_population <- function(parameters, nid)
{
	require(dplyr)

	Gx <- matrix(rbinom(parameters$nsnp_x * nid, 2, 0.5), nid, parameters$nsnp_x)

	Gy <- matrix(rbinom(parameters$nsnp_y * nid, 2, 0.5), nid, parameters$nsnp_y)

	U <- lapply(parameters$u, function(param)
	{
		G <- matrix(rbinom(param$nsnp_u * nid, 2, 0.5), nid, param$nsnp_u)
		u <- make_phen(param$eff_gu.u, G)
		return(list(p=u, G=G, nom=param$name_u))
	})


	bx <- parameters$eff_gx.x
	by <- parameters$eff_gx.y

	Gt <- Gx
	if(parameters$nsnp_y > 0)
	{
		bx <- c(bx, parameters$eff_gy.x)
		by <- c(by, parameters$eff_gy.y)
		Gt <- cbind(Gt, Gy)
	}

	if(length(parameters$u) > 0)
	{
		bx <- c(bx, sapply(parameters$u, function(x) x$eff_u.x))
		by <- c(by, sapply(parameters$u, function(x) x$eff_u.y))
		Gt <- cbind(Gt, do.call(cbind, lapply(U, function(x) x$p)))
	}
	by <- c(by, parameters$eff_x.y)

	x <- make_phen(bx, Gt)
	y <- make_phen(by, cbind(Gt, x))

	return(list(
		y=y,
		x=x,
		Gx=Gx,
		Gy=Gy,
		U=U
	))

}


#' Estimate the effects of all SNPs on all phenotypes
#'
#' @param sim Output from \code{simulate_population}
#'
#' @export
#' @return Lists of SNP-trait effect estimates
estimate_system_effects <- function(sim)
{
	# Get effects of all X SNPs

	l <- list()

	gx.x <- gwas(sim$x, sim$Gx)
	gx.x$inst <- "x"
	gx.x$snp <- paste0(1:nrow(gx.x), "x")
	
	gx.y <- gwas(sim$y, sim$Gx)
	gx.y$inst <- "x"
	gx.y$snp <- paste0(1:nrow(gx.y), "x")

	gx <- gx.x
	gy <- gx.y

	# Get effects of all SNPs on Y

	if(ncol(sim$Gy) > 0)
	{
		gy.x <- gwas(sim$x, sim$Gy)
		gy.x$inst <- "y"
		gy.x$snp <- paste0(1:nrow(gy.x), "y")
		gx <- rbind(gx, gy.x)

		gy.y <- gwas(sim$y, sim$Gy)
		gy.y$inst <- "y"
		gy.y$snp <- paste0(1:nrow(gy.y), "y")
		gy <- rbind(gy, gy.y)
	}

	nconf <- length(sim$U)
	if(nconf > 0)
	{
		gu.x <- list()
		gu.y <- list()
		gx.u <- list()
		gy.u <- list()
		gu.u <- list()
		for(i in 1:nconf)
		{
			gu.u[[i]] <- gwas(sim$U[[i]]$p, sim$U[[i]]$G)
			gu.u[[i]]$inst <- paste0("u", i)
			gu.u[[i]]$snp <- paste0(1:nrow(gu.u[[i]]), gu.u[[i]]$inst)

			gu.x[[i]] <- gwas(sim$x, sim$U[[i]]$G)
			gu.x[[i]]$inst <- paste0("u", i)
			gu.x[[i]]$snp <- paste0(1:nrow(gu.x[[i]]), gu.x[[i]]$inst)

			gx.u[[i]] <- gwas(sim$U[[i]]$p, sim$Gx)
			gx.u[[i]]$inst <- "x"
			gx.u[[i]]$snp <- paste0(1:nrow(gx.u[[i]]), gx.u[[i]]$inst)

			gu.y[[i]] <- gwas(sim$y, sim$U[[i]]$G)
			gu.y[[i]]$inst <- paste0("u", i)
			gu.y[[i]]$snp <- paste0(1:nrow(gu.y[[i]]), gu.y[[i]]$inst)

			if(ncol(sim$Gy) > 0)
			{
				gy.u[[i]] <- gwas(sim$U[[i]]$p, sim$Gy)
				gy.u[[i]]$inst <- "y"
				gy.u[[i]]$snp <- paste0(1:nrow(gy.u[[i]]), gy.u[[i]]$inst)
			}
		}
		gx <- rbind(gx, bind_rows(gu.x))
		gy <- rbind(gy, bind_rows(gu.y))
		gu <- list()
		for(i in 1:nconf)
		{
			if(ncol(sim$Gy) > 0)
			{
				gu[[i]] <- rbind(gx.u[[i]], gy.u[[i]], gu.u[[i]])
			} else {
				gu[[i]] <- rbind(gx.u[[i]], gy.u[[i]])
			}
		}
		names(gu) <- paste0("u", 1:nconf)
		l <- gu
	}
	l$x <- gx
	l$y <- gy
	return(l)
}




#' Wrapper for simulation pipeline
#'
#' Based on the parameters specified this function will call \code{init_parameters}, \code{add_u}, \code{sample_system_effects}, \code{simulate_population} and \code{estimate_system_effects}. A separate population is generated for each phenotype (x, y and each of u) to allow 2SMR and PRS analyses
#'
#' @param nidx 
#' @param nidy 
#' @param nidu=0 If 0 then don't simulate separate populations for the u variables
#' @param nu=0 Number of confounders influencing x and y
#' @param na=0 Number of traits upstream of x
#' @param nb=0 Number of traits upstream of y
#' @param var_x.y 
#' @param nsnp_x 
#' @param var_gx.x 
#' @param var_gx.y=0 
#' @param mu_gx.y=0 
#' @param prop_gx.y=1 
#' @param nsnp_y=0 
#' @param var_gy.y=0 
#' @param var_gy.x=0 
#' @param mu_gy.x=0 
#' @param prop_gy.x=1 
#'
#' @export
#' @return List of effects across system
create_system <- function(nidx, nidy, nidu=0, nu=0, na=0, nb=0, var_x.y, nsnp_x, var_gx.x, var_gx.y=0, mu_gx.y=0, prop_gx.y=1, nsnp_y=0, var_gy.y=0, var_gy.x=0, mu_gy.x=0, prop_gy.x=1)
{
	parameters <- init_parameters(nsnp_x=nsnp_x, nsnp_y=nsnp_y, var_gx.x=var_gx.x, var_gy.y=var_gy.y, var_x.y=var_x.y, var_gx.y=var_gx.y, mu_gx.y=mu_gx.y, prop_gx.y=prop_gx.y, var_gy.x=var_gy.x, mu_gy.x=mu_gy.x, prop_gy.x=prop_gy.x)
	if(nu > 0)
	{
		for(i in 1:nu)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=runif(1, min=0.01, max=0.1), 
				var_u.y=runif(1, min=0.01, max=0.1), 
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}
	if(na > 0)
	{
		for(i in 1:na)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=runif(1, min=0.01, max=0.1), 
				var_u.y=0,
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}
	if(nb > 0)
	{
		for(i in 1:nb)
		{
			parameters <- add_u(
				parameters, 
				nsnp_u=sample(5:30, 1), 
				var_u.x=0,
				var_u.y=runif(1, min=0.01, max=0.1),
				var_gu.u=runif(1, min=0.02, 0.2)
			)
		}
	}

	parameters <- sample_system_effects(parameters)

	message("X")
	pop <- simulate_population(parameters, nidx)
	x <- estimate_system_effects(pop)

	message("Y")
	pop <- simulate_population(parameters, nidy)
	y <- estimate_system_effects(pop)

	out <- list(x=x, y=y, parameters=parameters)

	if(nidu > 0)
	{
		u <- list()
		if(nu > 0)
		{
			for(i in 1:nu)
			{
				message("U: ", i, " of ", nu)
				pop <- simulate_population(parameters, nidu)
				u[[i]] <- estimate_system_effects(pop)
			}
		}

		a <- list()
		if(na > 0)
		{
			for(i in 1:na)
			{
				message("A: ", i, " of ", na)
				pop <- simulate_population(parameters, nidu)
				a[[i]] <- estimate_system_effects(pop)
			}
		}

		b <- list()
		if(nb > 0)
		{
			for(i in 1:nb)
			{
				message("B: ", i, " of ", nb)
				pop <- simulate_population(parameters, nidu)
				b[[i]] <- estimate_system_effects(pop)
			}
		}
		u <- c(u, a, b)
		names(u) <- paste0("u", 1:length(u))
		out$u <- u
	}
	return(out)
}


#' Apply MR tests to system
#'
#'
#' @param ss Output from create_syste,
#' @param id string denoting simulation ID
#'
#' @export
#' @return List
test_system <- function(ss, id="test")
{
	out <- list()

	dx <- make_dat(ss$x$x, ss$y$y)
	dx$exposure <- paste0("X:", id)
	dx$outcome <- paste0("Y:", id)
	dx$id.exposure <- paste0("X:", id)
	dx$id.outcome <- paste0("Y:", id)

	dy <- make_dat(ss$y$y, ss$x$x)
	dy$exposure <- paste0("Y:", id)
	dy$outcome <- paste0("X:", id)
	dy$id.exposure <- paste0("Y:", id)
	dy$id.outcome <- paste0("X:", id)

	# Oracle
	ox <- subset(dx, grepl("x", SNP))
	out$ox <- try(TwoSampleMR::mr_wrapper(ox)[[1]])
	oy <- subset(dy, grepl("y", SNP))
	out$oy <- try(TwoSampleMR::mr_wrapper(oy)[[1]])

	# Empirical
	ex <- subset(dx, pval.exposure < 5e-8)
	out$ex <- try(TwoSampleMR::mr_wrapper(ex)[[1]])
	ey <- subset(dy, pval.exposure < 5e-8)
	out$ey <- try(TwoSampleMR::mr_wrapper(ey)[[1]])

	param <- expand.grid(
		hypothesis = c("x", "y"), 
		selection = c("e", "o"),
		type = c("x", "y", "u")
	)
	o <- list()
	for(i in 1:nrow(param))
	{
		o[[i]] <- data_frame(
			hypothesis = param$hypothesis[i],
			selection = param$selection[i],
			type = param$type[i],
			measure = c("nofilter", "outlier", "steiger", "either", "both"),
			counts = get_counts(param$type[i], get(paste0(param$selection[i], param$hypothesis[i])), out[[paste0(param$selection[i], param$hypothesis[i])]])
		)
	}

	out$instrument_validity <- dplyr::bind_rows(o)
	out$instrument_validity$id <- id

	# Best model
	out$ex$estimates <- best_model(out$ex$estimates, ss$parameters$eff_x.y)
	out$ox$estimates <- best_model(out$ox$estimates, 0)

	return(out)
}

get_counts <- function(node, dat, res)
{
	c(sum(grepl(node, subset(dat)$SNP)),
	sum(grepl(node, subset(dat, SNP %in% subset(res$snps_removed, !outlier)$SNP)$SNP)),
	sum(grepl(node, subset(dat, SNP %in% subset(res$snps_removed, !steiger)$SNP)$SNP)),
	sum(grepl(node, subset(dat, SNP %in% subset(res$snps_removed, !either)$SNP)$SNP)),
	sum(grepl(node, subset(dat, SNP %in% subset(res$snps_removed, !both)$SNP)$SNP)))
}

best_model <- function(res, bxy)
{
	res$beta_correct <- res$ci_low <= bxy & res$ci_upp >= bxy
	res$beta_best <- FALSE
	res$beta_best[which.min(abs(res$b - bxy))] <- TRUE
	res$pval_sig <- res$pval < 0.05
	res$pval_lowest <- FALSE
	res$pval_lowest[which.min(res$pval)] <- TRUE
	res$pval_highest <- FALSE
	res$pval_highest[which.max(res$pval)] <- TRUE
	return(res)
}
