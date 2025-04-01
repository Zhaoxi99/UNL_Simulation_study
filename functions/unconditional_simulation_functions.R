library(sn)
library(nor1mix)
library(future)
library(future.apply)
##simpson's rule
simpson_inte <- function(x.vec,f.vec) {
  n1=length(x.vec)
  n2=length(f.vec)
  if(n1!=n2){
    stop("length of x.vec is not equal to f.vec")
  }
  n=n1
  if(n%%2==0){
    stop("length of x.vec and f.vec must be odd to have an even number of subintervals")
  }
  x1.ind=seq(2,n-1,by=2)
  x2.ind=seq(3,n-2,by=2)
  
  #print(length(x1.ind))
  #print(length(x2.ind))
  
  h=(x.vec[n]-x.vec[1])/(n-1)
  f_start=f.vec[1]
  f_end=f.vec[n]
  f.vec1 <- f.vec[x1.ind]
  f.vec2 <- f.vec[x2.ind]
  h/3*(f_start + f_end + 4*sum(f.vec1) + 2*sum(f.vec2))   # return value
}

dpm <- function(y, prior, mcmc, standardise = FALSE) {
  
  multinom <- function(probs) {
    probs <- t(apply(probs,1,cumsum))
    res <- rowSums(probs - runif(nrow(probs)) < 0) + 1
    return(res) 
  }
  
  n <- length(y)
  yt <- y
  if(standardise) {
    yt <- (y-mean(y))/sd(y)
  }
  
  m0 <- prior$m0
  S0 <- prior$S0
  a <- prior$a
  b <- prior$b
  alpha <- prior$alpha
  L <- prior$L
  
  nburn <- mcmc$nburn
  nsave <- mcmc$nsave
  nskip <- mcmc$nskip
  nsim <- nsave*nskip+nburn
  
  p <- ns <- rep(0,L)
  v <- rep(1/L,L)
  v[L] <- 1
  
  prop <- matrix(NA_real_, nrow = n, ncol = L)
  
  z <- matrix(NA_real_, nrow = nsim, ncol = n)
  z_tmp <- vector(length = n)
  
  z[1,] <- rep(1,n)
  
  P <- Mu <- Sigma2 <- matrix(0, nrow = nsim, ncol = L)
  
  Mu_tmp <- Sigma2_tmp  <- vector(length = L)
  
  Mu[1,] <- rep(mean(yt), L)
  Sigma2[1,] <- rep(var(yt), L)
  
  Mu_tmp <- Mu[1,]
  Sigma2_tmp <- Sigma2[1,]
  
  for(i in 2:nsim) {
    Sigma2_tmp <- Sigma2[i-1,]
    
    cumv <- cumprod(1-v)
    p[1] <- v[1]
    p[2:L] <- v[2:L]*cumv[1:(L-1)]
    
    for(l in 1:L){
      prop[,l] <- p[l]*dnorm(yt, mean = Mu_tmp[l], sd = sqrt(Sigma2_tmp[l]))
    }
    prob <- prop/rowSums(prop)
    
    z_tmp <- multinom(prob)
    ns <- sapply(1:L, function(x, v) sum(v == x), v = z_tmp)
    yt_z_l <- sapply(1:L, function(x, v, y) sum(y[v == x]), v = z_tmp, y = yt)
    
    v[1:(L-1)] <- rbeta(L-1, 1 + ns[1:(L-1)], alpha + rev(cumsum(rev(ns[-1]))))
    
    varmu <- 1/((1/S0) + (ns/Sigma2_tmp))
    meanmu <- ((yt_z_l/Sigma2_tmp) + (m0/S0))/((1/S0) + (ns/Sigma2_tmp))
    Mu_tmp <- rnorm(L, mean = meanmu, sd = sqrt(varmu))
    yt_z_l_mu <- sapply(1:L, function(x, v, y, mu) sum((y[v == x] - mu[x])^2), v = z_tmp, y = yt, mu = Mu_tmp)
    
    Sigma2_tmp <- 1/rgamma(L, a + ns/2, b + 0.5*yt_z_l_mu)
    
    P[i,] <- p
    z[i,] <- z_tmp
    Mu[i,] <- Mu_tmp
    Sigma2[i,] <- Sigma2_tmp
  }
  if(standardise){
    Mu <- sd(y)*Mu + mean(y)
    Sigma2 <- var(y)*Sigma2
  }
  
  res <- list()
  res$z <- z[seq(nburn+1, nsim, by = nskip),]
  res$P <- P[seq(nburn+1, nsim, by = nskip),]
  res$Mu <- Mu[seq(nburn+1, nsim, by = nskip),]
  res$Sigma2 <- Sigma2[seq(nburn+1, nsim, by = nskip),]
  return(res)
}

true_unl_uncon<-function(case=1,lgrid=10001){
  if(case==1){
    ##########Scenario 1#############
    ##########Case 1################
    mu1 <- -3.25; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 3.25; sd3 <- 1
    
    lt <- min(qnorm(0.00001, mu1, sd1), qnorm(0.00001, mu2, sd2), qnorm(0.00001, mu3, sd3))
    ut <- max(qnorm(0.99999, mu1, sd1), qnorm(0.99999, mu2, sd2), qnorm(0.99999, mu3, sd3))
    d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sd1)
    d2t <- dnorm(seq(lt, ut, len = lgrid), mu2, sd2)
    d3t <- dnorm(seq(lt, ut, len = lgrid), mu3, sd3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
    
  }else if (case==2){
    ##########Case 2################
    mu1 <- -1.3; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 1.15; sd3 <- 1
    
    lt <- min(qnorm(0.00001, mu1, sd1), qnorm(0.00001, mu2, sd2), qnorm(0.00001, mu3, sd3))
    ut <- max(qnorm(0.99999, mu1, sd1), qnorm(0.99999, mu2, sd2), qnorm(0.99999, mu3, sd3))
    d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sd1)
    d2t <- dnorm(seq(lt, ut, len = lgrid), mu2, sd2)
    d3t <- dnorm(seq(lt, ut, len = lgrid), mu3, sd3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
    
  }else if (case==3){
    ##########Case 3################
    mu1 <- -0.2; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 0.15; sd3 <- 1
    
    lt <- min(qnorm(0.00001, mu1, sd1), qnorm(0.00001, mu2, sd2), qnorm(0.00001, mu3, sd3))
    ut <- max(qnorm(0.99999, mu1, sd1), qnorm(0.99999, mu2, sd2), qnorm(0.99999, mu3, sd3))
    d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sd1)
    d2t <- dnorm(seq(lt, ut, len = lgrid), mu2, sd2)
    d3t <- dnorm(seq(lt, ut, len = lgrid), mu3, sd3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
    
  }else if (case==4){
    ##########Scenario 2#############
    ##########Case 1################
    lt <- min(qgamma(0.00001, 3, 1), qsn(0.00001, 6, 2, 5), qsn(0.00001, 8, 2, 5))
    ut <- max(qgamma(0.99999, 3, 1), qsn(0.99999, 6, 2, 5), qsn(0.99999, 8, 2, 5))
    d1t <- dgamma(seq(lt, ut, len = lgrid), 3, 1)
    d2t <- dsn(seq(lt, ut, len = lgrid), 6, 2, 5)
    d3t <- dsn(seq(lt, ut, len = lgrid), 8, 2, 5)
    dmaxt <- pmax(d1t, d2t, d3t)
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)
    
  }else if(case==5){
    ##########Scenario 2#############
    ##########Case 2###############
    lt <- min(qgamma(0.00001, 3, 1), qsn(0.00001, 2, 2.5, 5), qsn(0.00001, 4.25, 2, 5))
    ut <- max(qgamma(0.99999, 3, 1), qsn(0.99999, 2, 2.5, 5), qsn(0.99999, 4.25, 2, 5))
    d1t <- dgamma(seq(lt, ut, len = lgrid), 3, 1)
    d2t <- dsn(seq(lt, ut, len = lgrid), 2, 2.5, 5)
    d3t <- dsn(seq(lt, ut, len = lgrid), 4.25, 2, 5)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
    
  }else if (case==6){
    ##########Scenario 2#############
    ##########Case 3################
    lt <- min(qgamma(0.00001, 1.5, 1), qsn(0.00001, 0.1, 2, 5), qsn(0.00001, 0.25, 2, 5))
    ut <- max(qgamma(0.99999, 1.5, 1), qsn(0.99999, 0.1, 2, 5), qsn(0.99999, 0.25, 2, 5))
    d1t <- dgamma(seq(lt, ut, len = lgrid), 1.5, 1)
    d2t <- dsn(seq(lt, ut, len = lgrid), 0.1, 2, 5)
    d3t <- dsn(seq(lt, ut, len = lgrid), 0.25, 2, 5)
    dmaxt <- pmax(d1t, d2t, d3t)
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
  }else if(case==7){
    ##########Scenario 3#############
    ##########Case 1################
    aux_1 <- norMix(mu = c(-6, -3), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_2 <- norMix(mu = c(0.5, 3.25), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_3 <- norMix(mu = c(3.5, 6.25), sigma = c(1, 1), w = c(0.5, 0.5))
    
    lt <- min(qnorMix(0.00001, aux_1), qnorMix(0.00001, aux_2), qnorMix(0.00001, aux_3))
    ut <- max(qnorMix(0.99999, aux_1), qnorMix(0.99999, aux_2), qnorMix(0.99999, aux_3))
    d1t <- dnorMix(seq(lt, ut, len = lgrid), aux_1)
    d2t <- dnorMix(seq(lt, ut, len = lgrid), aux_2)
    d3t <- dnorMix(seq(lt, ut, len = lgrid), aux_3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
  }else if(case==8){
    ##########Scenario 3#############
    ##########Case 2################
    aux_1 <- norMix(mu = c(-2.25, 0.5), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_2 <- norMix(mu = c(2.75, 5.5), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_3 <- norMix(mu = c(3, 5.75), sigma = c(1, 1), w = c(0.5, 0.5))
    
    lt <- min(qnorMix(0.00001, aux_1), qnorMix(0.00001, aux_2), qnorMix(0.00001, aux_3))
    ut <- max(qnorMix(0.99999, aux_1), qnorMix(0.99999, aux_2), qnorMix(0.99999, aux_3))
    d1t <- dnorMix(seq(lt, ut, len = lgrid), aux_1)
    d2t <- dnorMix(seq(lt, ut, len = lgrid), aux_2)
    d3t <- dnorMix(seq(lt, ut, len = lgrid), aux_3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
  }else if (case==9){
    ##########Scenario 3#############
    ##########Case 3################
    aux_1 <- norMix(mu = c(0.15, 2.75), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_2 <- norMix(mu = c(0.5, 3), sigma = c(1, 1), w = c(0.5, 0.5))
    aux_3 <- norMix(mu = c(0.85, 3.15), sigma = c(1, 1), w = c(0.5, 0.5))
    
    lt <- min(qnorMix(0.00001, aux_1), qnorMix(0.00001, aux_2), qnorMix(0.00001, aux_3))
    ut <- max(qnorMix(0.99999, aux_1), qnorMix(0.99999, aux_2), qnorMix(0.99999, aux_3))
    d1t <- dnorMix(seq(lt, ut, len = lgrid), aux_1)
    d2t <- dnorMix(seq(lt, ut, len = lgrid), aux_2)
    d3t <- dnorMix(seq(lt, ut, len = lgrid), aux_3)
    dmaxt <- pmax(d1t, d2t, d3t)
    
    unl <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)  
  }
  return(unl)
}

####
gen_data_uncon<-function(case,n1=200,n2=200,n3=200,nrep=100){
  y1 <- matrix(0, nrow = n1, ncol = nrep)
  y2 <- matrix(0, nrow = n2, ncol = nrep)
  y3 <- matrix(0, nrow = n3, ncol = nrep)
  set.seed(123)
  if(case==1){
    mu1 <- -3.25; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 3.25; sd3 <- 1
    for(l in 1:nrep){
      y1[, l] <- rnorm(n1, mu1, sd1)
      y2[, l] <- rnorm(n2, mu2, sd2)
      y3[, l] <- rnorm(n3, mu3, sd3)
    }
  }else if(case==2){
    mu1 <- -1.3; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 1.15; sd3 <- 1
    for(l in 1:nrep){
      y1[, l] <- rnorm(n1, mu1, sd1)
      y2[, l] <- rnorm(n2, mu2, sd2)
      y3[, l] <- rnorm(n3, mu3, sd3)
    }
  }else if(case==3){
    mu1 <- -0.2; sd1 <- 1
    mu2 <- 0; sd2 <- 1
    mu3 <- 0.15; sd3 <- 1
    for(l in 1:nrep){
      y1[, l] <- rnorm(n1, mu1, sd1)
      y2[, l] <- rnorm(n2, mu2, sd2)
      y3[, l] <- rnorm(n3, mu3, sd3)
    }
  }else if(case==4){
    alpha_1 <- 3
    beta_1 <- 1
    xi_2<-6
    omega_2<-2
    alpha_2<-5
    xi_3<-8
    omega_3<-2
    alpha_3<-5
    
    for(l in 1:nrep){
      y1[, l] <- rgamma(n1, shape=alpha_1, rate=beta_1)
      y2[, l] <- rsn(n2, xi=xi_2, omega = omega_2,alpha = alpha_2)
      y3[, l] <- rsn(n3, xi=xi_3, omega = omega_3,alpha = alpha_3)
    }
  }else if(case==5){
    #function of minimum density by denfinition
    
    alpha_1 <- 3
    beta_1 <- 1
    xi_2<-2
    omega_2<-2.5
    alpha_2<-5
    xi_3<-4.25
    omega_3<-2
    alpha_3<-5
    
    for(l in 1:nrep){
      y1[, l] <- rgamma(n1, shape=alpha_1, rate=beta_1)
      y2[, l] <- rsn(n2, xi=xi_2, omega = omega_2,alpha = alpha_2)
      y3[, l] <- rsn(n3, xi=xi_3, omega = omega_3,alpha = alpha_3)
    }
  }else if(case==6){
    #function of minimum density by denfinition
    
    alpha_1 <- 1.5
    beta_1 <- 1
    xi_2<-0.1
    omega_2<-2
    alpha_2<-5
    xi_3<-0.25
    omega_3<-2
    alpha_3<-5
    
    for(l in 1:nrep){
      y1[, l] <- rgamma(n1, shape=alpha_1, rate=beta_1)
      y2[, l] <- rsn(n2, xi=xi_2, omega = omega_2,alpha = alpha_2)
      y3[, l] <- rsn(n3, xi=xi_3, omega = omega_3,alpha = alpha_3)
    }
  }else if(case==7){
    #function of minimum density by denfinition
    mu1a <- -6
    sigma1a <- 1
    mu1b <- -3
    sigma1b <- 1
    
    mu2a <- 0.5
    sigma2a <- 1
    mu2b <- 3.25
    sigma2b <- 1
    
    mu3a <- 3.5
    sigma3a <- 1
    mu3b <- 6.25
    sigma3b <- 1
    
    aux1 <- norMix(mu = c(mu1a,mu1b), sigma = c(sigma1a,sigma1b), w = c(0.5,0.5))
    aux2 <- norMix(mu = c(mu2a,mu2b), sigma = c(sigma2a,sigma2b), w = c(0.5,0.5))
    aux3 <- norMix(mu = c(mu3a,mu3b), sigma = c(sigma3a,sigma3b), w = c(0.5,0.5))
    
    for(l in 1:nrep){
      y1[, l] <- rnorMix(n1, obj = aux1)
      y2[, l] <- rnorMix(n2, obj = aux2)
      y3[, l] <- rnorMix(n3, obj = aux3)
    }
    
  }else if(case==8){
    #function of minimum density by denfinition
    mu1a <- -2.25
    sigma1a <- 1
    mu1b <- 0.5
    sigma1b <- 1
    
    mu2a <- 2.75
    sigma2a <- 1
    mu2b <- 5.5
    sigma2b <- 1
    
    mu3a <- 3
    sigma3a <- 1
    mu3b <- 5.75
    sigma3b <- 1
    
    aux1 <- norMix(mu = c(mu1a,mu1b), sigma = c(sigma1a,sigma1b), w = c(0.5,0.5))
    aux2 <- norMix(mu = c(mu2a,mu2b), sigma = c(sigma2a,sigma2b), w = c(0.5,0.5))
    aux3 <- norMix(mu = c(mu3a,mu3b), sigma = c(sigma3a,sigma3b), w = c(0.5,0.5))
    
    for(l in 1:nrep){
      y1[, l] <- rnorMix(n1, obj = aux1)
      y2[, l] <- rnorMix(n2, obj = aux2)
      y3[, l] <- rnorMix(n3, obj = aux3)
    }
  }else if(case==9){
    #function of minimum density by denfinition
    mu1a <- 0.15
    sigma1a <- 1
    mu1b <- 2.75
    sigma1b <- 1
    
    mu2a <- 0.5
    sigma2a <- 1
    mu2b <- 3
    sigma2b <- 1
    
    mu3a <- 0.85
    sigma3a <- 1
    mu3b <- 3.15
    sigma3b <- 1
    
    aux1 <- norMix(mu = c(mu1a,mu1b), sigma = c(sigma1a,sigma1b), w = c(0.5,0.5))
    aux2 <- norMix(mu = c(mu2a,mu2b), sigma = c(sigma2a,sigma2b), w = c(0.5,0.5))
    aux3 <- norMix(mu = c(mu3a,mu3b), sigma = c(sigma3a,sigma3b), w = c(0.5,0.5))
    
    for(l in 1:nrep){
      y1[, l] <- rnorMix(n1, obj = aux1)
      y2[, l] <- rnorMix(n2, obj = aux2)
      y3[, l] <- rnorMix(n3, obj = aux3)
    }
  }
  return(list(y1=y1,y2=y2,y3=y3))
}


#case=1;n1=200;n2=500;n3=100;nrep=10
simulation_unl_uncon<-function(case,nlist,nrep,nsave=5000,nburn=2000){
  n1=nlist$n1;n2=nlist$n2;n3=nlist$n3
  set.seed(25)
  data_con=gen_data_uncon(case = case,n1=n1,n2=n2,n3=n3,nrep=nrep)
  
  prior <- list(m0 = 0, S0 = 10, a = 2, b = 0.5, alpha = 1, L = 20)
  mcmc <- list(nsave = nsave, nburn = nburn, nskip = 1)
  lgrid <- 501
  
  nrep_list <- list()
  for (i in 1:nrep) {
    nrep_list[[as.character(i)]] <- i
  }
  Est_UNL=future_lapply(X=nrep_list,FUN=fit_and_cal_unl,data_con=data_con,prior=prior,mcmc=mcmc,
                        lgrid=lgrid,future.seed = TRUE)
  return(Est_UNL)
}

#l=1
fit_and_cal_unl<-function(l,data_con,prior,mcmc,lgrid){
  data_fitting1=data.frame(y=data_con$y1[,l])
  data_fitting2=data.frame(y=data_con$y2[,l])
  data_fitting3=data.frame(y=data_con$y3[,l])
  
  y_merge=c(data_con$y1[,l],data_con$y2[,l],data_con$y3[,l])
  data_scaled1=data.frame(scale(data_fitting1))
  data_scaled2=data.frame(scale(data_fitting2))
  data_scaled3=data.frame(scale(data_fitting3))
  
  y_grid <- seq(from=min(y_merge)-1,to=max(y_merge)+1,length.out=lgrid)
  y_grid1=(y_grid-mean(data_fitting1$y))/sd(data_fitting1$y)
  y_grid2=(y_grid-mean(data_fitting2$y))/sd(data_fitting2$y)
  y_grid3=(y_grid-mean(data_fitting3$y))/sd(data_fitting3$y)
  
  sd_y1=sd(data_fitting1$y)
  sd_y2=sd(data_fitting2$y)
  sd_y3=sd(data_fitting3$y)
  
  set.seed(123)
  res1 <- dpm(y = data_scaled1$y, prior = prior, mcmc = mcmc, standardise = FALSE) 
  set.seed(123)
  res2 <- dpm(y = data_scaled2$y, prior = prior, mcmc = mcmc, standardise = FALSE) 
  set.seed(123)
  res3 <- dpm(y = data_scaled3$y, prior = prior, mcmc = mcmc, standardise = FALSE) 
  
  
  niter <- nrow(res1$Mu) 
  Est_UNL1 <- numeric(niter)
  Est_UNL2 <- numeric(niter)
  
  for(r in 1:niter){
    aux_dens_1 <- norMix(mu = res1$Mu[r,], sigma = sqrt(res1$Sigma2[r,]), w = res1$P[r,])
    aux_dens_2 <- norMix(mu = res2$Mu[r,], sigma = sqrt(res2$Sigma2[r,]), w = res2$P[r,])
    aux_dens_3 <- norMix(mu = res3$Mu[r,], sigma = sqrt(res3$Sigma2[r,]), w = res3$P[r,])
    
    dens_1 <- dnorMix(x = y_grid1, aux_dens_1)/sd_y1
    dens_2 <- dnorMix(x = y_grid2, aux_dens_2)/sd_y2
    dens_3 <- dnorMix(x = y_grid3, aux_dens_3)/sd_y3
    
    dens_max <- pmax(dens_1, dens_2, dens_3)
    Est_UNL1[r] <- simpson_inte(y_grid, dens_max)
    
    inte_unl=simpson_inte(y_grid, pmax(dens_1,dens_2,dens_3))
    Est_UNL2[r]=1+0.25*inte_unl
    
  }
  Est_UNL=cbind(Est_UNL1,Est_UNL2)
  return(Est_UNL)
}



tranform_list_simu<-function(simulation_object,nsave=5000,nrep=100){
  Est_UNL1=matrix(0,nrow=nsave,ncol=nrep)
  Est_UNL2=matrix(0,nrow=nsave,ncol=nrep)
  counter=1
  for (i in names(simulation_object)){
    Est_UNL1[,counter]=simulation_object[[i]][,1]
    Est_UNL2[,counter]=simulation_object[[i]][,2]
    counter=counter+1
  }
  return(list(Est_UNL1=Est_UNL1,Est_UNL2=Est_UNL2))
}