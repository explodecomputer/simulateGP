#' Wrapper for generating a summary set 
#'
#' Allows arbitrary sample overlap
#'
#' @param beta_gx Array of true effects on x
#' @param beta_gy Array of true effects on y
#' @param af Array of effect allele frequencies
#' @param n_gx Sample size of g-x association
#' @param n_gy Sample size of g-y association
#' @param n_overlap Number of overlapping samples
#' @param cor_xy Observational correlation between x and y
#' @param prev_y Disease prevalence of y. Default = NA (in which case treated as continuous)
#' @param sigma_x SD of x. Default=1 
#' @param sigma_y SD of y. Default=1
#'
#' @export
#' @return Data frame of summary statistics for x and y
summary_set <- function(beta_gx, beta_gy, af, n_gx, n_gy, n_overlap, cor_xy, prev_y=NA, sigma_x=1, sigma_y=1)
{
  stopifnot(length(beta_gy) == length(beta_gx))
  stopifnot(length(af) == length(beta_gx))
  stopifnot(all(af < 1 & af > 0))
  stopifnot(cor_xy <= 1 & cor_xy >= -1)
  stopifnot(n_overlap <= min(n_gx, n_gy))

  nsnp <- length(beta_gx)

  if(n_overlap == 0)
  {
    bhat_gx <- generate_gwas_ss(beta=beta_gx, af=af, nid=n_gx, vy=sigma_x^2)
    bhat_gy <- generate_gwas_ss(beta=beta_gy, af=af, nid=n_gy, vy=sigma_y^2)
    dat <- merge_exp_out(bhat_gx, bhat_gy)
    return(dat)
  }

  G_prob <- cbind(af^2, 2*af*(1-af), (1-af)^2)

  if(!is.na(prev_y))
  {
    myfunc <- function(af, Gamma0, Gamma1, prev)
    {  
      af^2 * stats::plogis(Gamma0 + Gamma1 * 2) + 2 * af * (1 - af) * stats::plogis(Gamma0 + Gamma1) + (1 - af)^2 * stats::plogis(Gamma0) - prev 
    }

    Gamma0 <- numeric(nsnp)
    for(i in 1:nsnp)
    {
      Gamma0[i] <- stats::uniroot(myfunc, Gamma1=beta_gy[i], af=af[i], prev=prev_y, lower=-10, upper=10)$root
    }

    var_gy <- asymp_var_logistic(n_gy, G_prob, Gamma0, beta_gy)
    var_gx <- asymp_var_linear(n_gx, G_prob, sigma=sigma_x)
    cov_gx_gy <- asymp_cov_linear_logistic(n_overlap, n_gx, n_gy, G_prob, cor_xy=cor_xy, sigma=sigma_x, Gamma_0=Gamma0, Gamma_1=beta_gy, prev=prev_y)
  }

  var_gy <- asymp_var_linear(n_gy, G_prob, sigma=sigma_y)
  var_gx <- asymp_var_linear(n_gx, G_prob, sigma=sigma_x)
  cov_gx_gy <- asymp_cov_linear_linear(n_overlap, n_gx, n_gy, G_prob, sigma_x=sigma_x, sigma_y=sigma_y,cor_xy=cor_xy)

  cov_array <- array(dim=c(2, 2, nsnp))
  cov_array[1,1,] <- var_gx
  cov_array[2,1,] <- cov_gx_gy
  cov_array[1,2,] <- cov_array[2,1,]
  cov_array[2,2,] <- var_gy

  summary_stats <- apply(cov_array, 3, function(x)
  {
    MASS::mvrnorm(n=1, mu=c(0,0), Sigma=x)
  })

  summary_stats <- t(summary_stats + rbind(beta_gx, beta_gy))

  dat <- dplyr::tibble(
    SNP = 1:nsnp,
    id.exposure="X",
    id.outcome="Y",
    exposure="X",
    outcome="Y",
    beta.exposure = summary_stats[,1],
    beta.outcome = summary_stats[,2],
    se.exposure = sqrt(var_gx),
    se.outcome = sqrt(var_gy),
    fval.exposure = (beta.exposure/se.exposure)^2,
    fval.outcome = (beta.outcome/se.outcome)^2,
    pval.exposure = pf(fval.exposure, df1=1, df2=n_gx-1, lower.tail=FALSE),
    pval.outcome = pf(fval.outcome, df1=1, df2=n_gy-1, lower.tail=FALSE),
    eaf.exposure = af,
    eaf.outcome = af,
    mr_keep=TRUE
  )
  return(dat)
}



## calculate asymptotic version of (X^t W X) 
asymp_var_logistic <- function(n,G_prob,Gamma_0,Gamma_1){
  N <- length(Gamma_0)
  disease_probs <- stats::plogis(outer(Gamma_0,c(1,1,1)) +  outer(Gamma_1,c(0,1,2)))
  diag_weights <- disease_probs*(1-disease_probs)
    
  a <- n*apply(G_prob*diag_weights,1,sum)
  b <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b 
  d <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)                     
  
  ##  invert matrix and take element for Gamma1
  return(a/(a*d-b*c))

}

asymp_var_linear <- function(n,G_prob,sigma=1){
  N <- nrow(G_prob)
 
  a <- n*apply(G_prob,1,sum)
  b <- n*apply(G_prob*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b 
  d <- n*apply(G_prob*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)                     
  
  ##  invert matrix and take element for Gamma1
  return(a/(a*d-b*c)*sigma^2)
  
}


asymp_cov_linear_linear <- function(n_overlap,n_gx,n_gy,G_prob,sigma_x=1,sigma_y=1,cor_xy=.5){
  N <- nrow(G_prob)
    cov_xy <- cor_xy*sigma_x*sigma_y
  
  a <- apply(G_prob,1,sum)
  b <- apply(G_prob*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b 
  d <- apply(G_prob*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)                     
  
  ##  invert matrix and take element for Gamma1
  return(n_overlap*(a/(a*d-b*c))*cov_xy/(n_gx*n_gy))
  
}

asymp_cov_linear_logistic <- function(n_overlap,n_gx,n_gy,G_prob,cor_xy=0,sigma=1,Gamma_0,Gamma_1,prev=0.01){
  N <- length(Gamma_0)
  cov_xy <- cor_xy*sigma*sqrt(prev*(1-prev))
  disease_probs <- stats::plogis(outer(Gamma_0,c(1,1,1)) +  outer(Gamma_1,c(0,1,2)))
  diag_weights <- disease_probs*(1-disease_probs)
  
  a <- apply(G_prob*diag_weights,1,sum)
  b <- apply(G_prob*diag_weights*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b 
  d <- apply(G_prob*diag_weights*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)                     
  
  ##  invert matrix and take element for Gamma1
  return(n_overlap*(a/(a*d-b*c))*cov_xy/(n_gx*n_gy))
  
}

