library(splines)
library(LSBP)
library(nor1mix)
library(future)
library(future.apply)



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

true_unl_con<-function(case=1,xgrid,lgrid=10001){
  if(case==1){
    unl=rep(0,length(xgrid))
    for(i in 1:length(xgrid)){
      mu1=0.25+xgrid[i];sd1=1
      mu2=1+1.5*xgrid[i];sd2=1.5
      mu3=2.5+4*xgrid[i];sd3=1.75
      lt <- min(qnorm(0.00001, mean=mu1,sd=sd1), qnorm(0.00001, mean=mu2,sd=sd2), 
                qnorm(0.00001, mean=mu3,sd=sd3))
      ut <- max(qnorm(0.99999, mean=mu1,sd=sd1), qnorm(0.99999, mean=mu2,sd=sd2), 
                qnorm(0.99999, mean=mu3,sd=sd3))
      d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sd1)
      d2t <- dnorm(seq(lt, ut, len = lgrid), mu2, sd2)
      d3t <- dnorm(seq(lt, ut, len = lgrid), mu3, sd3)
      dmaxt <- pmax(d1t, d2t, d3t)
      unl[i] <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)
    }
  }else if (case==2){
    unl=rep(0,length(xgrid))
    for(i in 1:length(xgrid)){
      mu1= -0.75 + sin(pi*xgrid[i] + 1.25);sd1=0.5
      mu2=0.75 + sin(pi*xgrid[i]);sd2=1.25 + xgrid[i]^2
      mu3=2.35 + xgrid[i]^2;sd3=1
      
      lt <- min(qnorm(0.00001, mean=mu1,sd=sd1), qnorm(0.00001, mean=mu2,sd=sd2), 
                qnorm(0.00001, mean=mu3,sd=sd3))
      ut <- max(qnorm(0.99999, mean=mu1,sd=sd1), qnorm(0.99999, mean=mu2,sd=sd2), 
                qnorm(0.99999, mean=mu3,sd=sd3))
      d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sd1)
      d2t <- dnorm(seq(lt, ut, len = lgrid), mu2, sd2)
      d3t <- dnorm(seq(lt, ut, len = lgrid), mu3, sd3)
      dmaxt <- pmax(d1t, d2t, d3t)
      unl[i] <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)
      
    }
  }else if(case==3){
    unl=rep(0,length(xgrid))
    for(i in 1:length(xgrid)){
      mu1 <- -0.75 + sin(pi*xgrid[i])
      #sigma1 <- sqrt(exp(xgrid[i]))
      sigma1 <- 1
      
      mu2_1 <- xgrid[i]
      mu2_2 <- xgrid[i]^2
      w2_1 <- exp(xgrid[i])/(1 + exp(xgrid[i]))
      w2_2 <- 1/(1 + exp(xgrid[i]))
      sigma2_1 <- 0.5
      sigma2_2 <- 0.75
      aux <- norMix(mu = c(mu2_1, mu2_2),
                    sigma = c(sigma2_1, sigma2_2),
                    w = c(w2_1, w2_2))
      
      shape3 <- 3 + xgrid[i]^2
      rate3 <- 0.5  + exp(xgrid[i])
      
      lt <- min(qnorm(0.0001, mu1, sigma1), qnorMix(0.0001, aux), qgamma(0.0001, shape = shape3, rate = rate3))
      ut <- max(qnorm(0.9999, mu1, sigma1), qnorMix(0.9999, aux), qgamma(0.9999, shape = shape3, rate = rate3))
      
      d1t <- dnorm(seq(lt, ut, len = lgrid), mu1, sigma1)
      d2t <- dnorMix(seq(lt, ut, len = lgrid), aux)
      d3t <- dgamma(seq(lt, ut, len = lgrid), shape = shape3, rate = rate3)
      dmaxt <- pmax(d1t, d2t, d3t)
      unl[i] <- simpson_inte(seq(lt, ut, len = lgrid), dmaxt)
    }
  }
  return(unl)
}

###real simulation
gen_data_con<-function(case=1,n1=200,n2=200,n3=200,nrep=100){
  set.seed(179)
  y1 <- matrix(0, nrow = n1, ncol = nrep)
  y2 <- matrix(0, nrow = n2, ncol = nrep)
  y3 <- matrix(0, nrow = n3, ncol = nrep)
  
  x1 <- matrix(0, nrow = n1, ncol = nrep)
  x2 <- matrix(0, nrow = n2, ncol = nrep)
  x3 <- matrix(0, nrow = n3, ncol = nrep)
  for(l in 1:nrep){
    x1[, l]=runif(n=n1,min=-1,max=1)
    x2[, l]=runif(n=n2,min=-1,max=1)
    x3[, l]=runif(n=n3,min=-1,max=1)
  }
  if(case==1){
    for(l in 1:nrep){
      mu1=0.25+x1[, l];sd1=1
      mu2=1+1.5*x2[, l];sd2=1.5
      mu3=2.5+4*x3[, l];sd3=1.75
      y1[, l]=rnorm(n1,mean=mu1,sd=sd1)
      y2[, l]=rnorm(n2,mean=mu2,sd=sd2)
      y3[, l]=rnorm(n3,mean=mu3,sd=sd3)
    }
  }else if(case==2){
    for(l in 1:nrep){
      mu1= -0.75 + sin(pi*x1[, l] + 1.25);sd1=0.5
      mu2=0.75 + sin(pi*x2[, l]);sd2=1.25 + x2[, l]^2
      mu3=2.35 + x3[, l]^2;sd3=1
      
      y1[, l]=rnorm(n1,mean=mu1,sd=sd1)
      y2[, l]=rnorm(n2,mean=mu2,sd=sd2)
      y3[, l]=rnorm(n3,mean=mu3,sd=sd3)
    }
  }else if(case==3){
    for(l in 1:nrep){
      
      mu1 <- -0.75 + sin(pi*x1[, l])
      sigma1 <- 1
      
      mu2_1 <- x2[, l]
      mu2_2 <- x2[, l]^2
      w2_1 <- exp(x2[, l])/(1 + exp(x2[, l]))
      w2_2 <- 1/(1 + exp(x2[, l]))
      sigma2_1 <- 0.5
      sigma2_2 <- 0.75
      
      shape3 <- 3 + x3[, l]^2
      rate3 <- 0.5  + exp(x3[, l])
      
      y1[, l]=rnorm(n1, mean=mu1, sd=sigma1)
      y3[, l]=rgamma(n3, shape = shape3, rate = rate3)
      for(i in 1:n2){
        aux2=norMix(mu = c(mu2_1[i], mu2_2[i]),
                    sigma = c(sigma2_1, sigma2_2),
                    w = c(w2_1[i], w2_2[i]))
        y2[i,l]=rnorMix(n=1, obj = aux2)
      }
    }
  }
  return(list(y1=y1,y2=y2,y3=y3,x1=x1,x2=x2,x3=x3))
}



cal_density<-function(fit_Gibbs1,fit_Gibbs2,fit_Gibbs3,
                      R,X1_1,X2_1,X1_2,X2_2,X1_3,X2_3,
                      y_grid1,y_grid2,y_grid3,
                      sd_y1,sd_y2,sd_y3){
  
  # Posterior density - Gibbs sampling
  #print("b")
  pred_Gibbs1 <- matrix(0,R,length(y_grid1))
  pred_Gibbs2 <- matrix(0,R,length(y_grid2))
  pred_Gibbs3 <- matrix(0,R,length(y_grid3))
  
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    #r=2
    for(i in 1:length(y_grid1)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs1[r,i] <- c(LSBP_density(y_grid1[i],t(matrix(X1_1[i,])),t(matrix(X2_1[i,])),
                                         fit_Gibbs1$param$beta_mixing[r,,],
                                         fit_Gibbs1$param$beta_kernel[r,,],
                                         fit_Gibbs1$param$tau[r,])) /sd_y1
    }
    
    for(i in 1:length(y_grid2)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs2[r,i] <- c(LSBP_density(y_grid2[i],t(matrix(X1_2[i,])),t(matrix(X2_2[i,])),
                                         fit_Gibbs2$param$beta_mixing[r,,],
                                         fit_Gibbs2$param$beta_kernel[r,,],
                                         fit_Gibbs2$param$tau[r,])) /sd_y2
    }
    
    for(i in 1:length(y_grid3)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs3[r,i] <- c(LSBP_density(y_grid3[i],t(matrix(X1_3[i,])),t(matrix(X2_3[i,])),
                                         fit_Gibbs3$param$beta_mixing[r,,],
                                         fit_Gibbs3$param$beta_kernel[r,,],
                                         fit_Gibbs3$param$tau[r,])) /sd_y3
    }
  }
  
  return(list(pred_Gibbs1=pred_Gibbs1,pred_Gibbs2=pred_Gibbs2,pred_Gibbs3=pred_Gibbs3))
}

cal_density_normal<-function(results1_mat,results2_mat,results3_mat,
                             R,X1_1,X1_2,X1_3,
                             y_grid1,y_grid2,y_grid3,
                             sd_y1,sd_y2,sd_y3){
  
  # Posterior density - Gibbs sampling
  # results1_mat=results1_HCI_mat;results2_mat=results2_HCI_mat;results3_mat=results3_HCI_mat
  # y_grid1=adni_1_scaled$HCI;y_grid2=adni_2_scaled$HCI;y_grid3=adni_3_scaled$HCI
  # sd_y1=sd_HCI1;sd_y2=sd_HCI2;sd_y3=sd_HCI3
  pred_Gibbs1 <- matrix(0,R,length(y_grid1))
  pred_Gibbs2 <- matrix(0,R,length(y_grid2))
  pred_Gibbs3 <- matrix(0,R,length(y_grid3))
  
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    #r=2
    for(i in 1:length(y_grid1)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs1[r,i] <- c(dnorm(y_grid1[i],
                                  mean=sum(X1_1[i,]*results1_mat[r,1:2]),
                                  sd=results1_mat[r,3])) /sd_y1
    }
    for(i in 1:length(y_grid2)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs2[r,i] <- c(dnorm(y_grid2[i],
                                  mean=sum(X1_2[i,]*results2_mat[r,1:2]),
                                  sd=results2_mat[r,3])) /sd_y2
    }
    for(i in 1:length(y_grid3)){  # Cycle over the GAD grid
      #i=4
      pred_Gibbs3[r,i] <- c(dnorm(y_grid3[i],
                                  mean=sum(X1_3[i,]*results3_mat[r,1:2]),
                                  sd=results3_mat[r,3])) /sd_y3
    }
  }
  return(list(pred_Gibbs1=pred_Gibbs1,pred_Gibbs2=pred_Gibbs2,pred_Gibbs3=pred_Gibbs3))
}

lpml_density<-function(dens_mat){
  #dens_mat is a matrix with nrow = niter, ncol = n
  niter=dim(dens_mat)[1];n=dim(dens_mat)[2]
  aux <- 1/dens_mat
  omegabari <- apply(aux, 2, mean)
  omegabari_1 <- sqrt(niter) * omegabari
  omegatilde <- matrix(0, nrow = niter, ncol = n)
  
  for(i in 1:n) {
    omegatilde[,i] <- pmin(aux[,i], omegabari_1[i])  
  }
  
  sum_omegatilde <- apply(omegatilde,2,sum)
  sum_term_omegatilde <- apply(dens_mat*omegatilde, 2, sum)
  cpo <- sum_term_omegatilde/sum_omegatilde
  
  lpml <- sum(log(cpo))
  return(lpml)
}

waic_density<-function(dens_mat){
  #dens_mat is a matrix with nrow = niter, ncol = n
  niter=dim(dens_mat)[1];n=dim(dens_mat)[2]
  logdens_mat <- log(dens_mat)
  lpd <- sum(log(apply(exp(logdens_mat),2,mean)))
  p2 <- sum(apply(logdens_mat,2,var))
  waic <- -2*(lpd-p2)
  return(waic)
}


fit_and_cal_unl_con_best<-function(l,data_con,R,burn_in,prior1,prior2,prior3,xlength,H,
                                   model_formula_y1,model_formula_y2,model_formula_y3,case=1){
  #l=3
  data_fitting1=data.frame(y=data_con$y1[,l],x=data_con$x1[,l])
  data_fitting2=data.frame(y=data_con$y2[,l],x=data_con$x2[,l])
  data_fitting3=data.frame(y=data_con$y3[,l],x=data_con$x3[,l])
  
  
  if(case==1){
    x_bs1 <- bs(data_con$x1[,l], intercept = FALSE, degree = 3, knots = c())
    x_bs2 <- bs(data_con$x2[,l], intercept = FALSE, degree = 3, knots = c())
    x_bs3 <- bs(data_con$x3[,l], intercept = FALSE, degree = 3, knots = c()) 
  }else if(case==2){
    x_bs1 <- bs(data_con$x1[,l], intercept = FALSE, degree = 3, knots = quantile(data_fitting1$x,probs = c(0.5)))
    x_bs2 <- bs(data_con$x2[,l], intercept = FALSE, degree = 3, knots = c())
    x_bs3 <- bs(data_con$x3[,l], intercept = FALSE, degree = 3, knots = c()) 
  }else if(case==3){
    x_bs1 <- bs(data_con$x1[,l], intercept = FALSE, degree = 3, knots = c())
    x_bs2 <- bs(data_con$x2[,l], intercept = FALSE, degree = 3, knots = c())
    x_bs3 <- bs(data_con$x3[,l], intercept = FALSE, degree = 3, knots = c()) 
  }
  
  data_fitting1=data.frame(data_fitting1,BS=x_bs1)
  data_fitting2=data.frame(data_fitting2,BS=x_bs2)
  data_fitting3=data.frame(data_fitting3,BS=x_bs3)
  
  y_merge=c(data_con$y1[,l],data_con$y2[,l],data_con$y3[,l])
  data_scaled1=data.frame(scale(data_fitting1))
  data_scaled2=data.frame(scale(data_fitting2))
  data_scaled3=data.frame(scale(data_fitting3))
  
  # Gibbs algorithm
  set.seed(23) # The seed is setted so that the Gibbs sampler is reproducible.
  fit_Gibbs1   <- LSBP_Gibbs(model_formula_y1, data=data_scaled1, H=H, prior=prior1, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs2   <- LSBP_Gibbs(model_formula_y2, data=data_scaled2, H=H, prior=prior2, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs3   <- LSBP_Gibbs(model_formula_y3, data=data_scaled3, H=H, prior=prior3, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  
  # Points for which we will evaluate
  x.points =seq(from=-1,to=1,length.out=xlength)
  x.points1=(x.points-mean(data_fitting1$x))/sd(data_fitting1$x)
  x.points2=(x.points-mean(data_fitting2$x))/sd(data_fitting2$x)
  x.points3=(x.points-mean(data_fitting3$x))/sd(data_fitting3$x)
  
  if(case==1){
    x_bs.points1=x.points1
    x_bs.points2=x.points2
    x_bs.points3=x.points3
  }else if(case==2){
    x1_bs.points <- as.matrix(bs(x.points, intercept = FALSE, degree = 3, 
                                 knots = quantile(data_fitting1$x,probs = c(0.5))))
    x_bs.points1=cbind((x1_bs.points[,1]-mean(data_fitting1$BS.1))/sd(data_fitting1$BS.1),
                       (x1_bs.points[,2]-mean(data_fitting1$BS.2))/sd(data_fitting1$BS.2),
                       (x1_bs.points[,3]-mean(data_fitting1$BS.3))/sd(data_fitting1$BS.3),
                       (x1_bs.points[,4]-mean(data_fitting1$BS.4))/sd(data_fitting1$BS.4))
    x_bs.points2=x.points2
    x_bs.points3=x.points3
  }else if(case==3){
    x1_bs.points <- as.matrix(bs(x.points, intercept = FALSE, degree = 3, 
                                 knots = c()))
    x_bs.points1=cbind((x1_bs.points[,1]-mean(data_fitting1$BS.1))/sd(data_fitting1$BS.1),
                       (x1_bs.points[,2]-mean(data_fitting1$BS.2))/sd(data_fitting1$BS.2),
                       (x1_bs.points[,3]-mean(data_fitting1$BS.3))/sd(data_fitting1$BS.3))
    x_bs.points2=x.points2
    x_bs.points3=x.points3
  }
  
  # And the correspondings design matrices
  
  X1_1=cbind(1,x_bs.points1)
  X2_1=cbind(1,x.points1)
  
  X1_2=cbind(1,x_bs.points2)
  X2_2=cbind(1,x.points2)
  
  X1_3=cbind(1,x_bs.points3)
  X2_3=cbind(1,x.points3)
  
  
  y_grid <- seq(from=min(y_merge)-0.5,to=max(y_merge)+0.5,length.out=501)
  y_grid1=(y_grid-mean(data_fitting1$y))/sd(data_fitting1$y)
  y_grid2=(y_grid-mean(data_fitting2$y))/sd(data_fitting2$y)
  y_grid3=(y_grid-mean(data_fitting3$y))/sd(data_fitting3$y)
  
  sd_y1=sd(data_fitting1$y)
  sd_y2=sd(data_fitting2$y)
  sd_y3=sd(data_fitting3$y)
  
  
  # Posterior density - Gibbs sampling
  #print("b")
  pred_Gibbs1 <- matrix(0,length(y_grid),dim(X1_1)[1])
  pred_Gibbs2 <- matrix(0,length(y_grid),dim(X1_2)[1])
  pred_Gibbs3 <- matrix(0,length(y_grid),dim(X1_3)[1])
  
  Est_UNL <- matrix(0,R,dim(X1_3)[1])
  #print("a")
  
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    #r=2
    for(i in 1:length(y_grid)){  # Cycle over the y grid
      #i=4
      pred_Gibbs1[i,] <- c(LSBP_density(y_grid1[i],X1_1,X2_1,
                                        fit_Gibbs1$param$beta_mixing[r,,],
                                        fit_Gibbs1$param$beta_kernel[r,,],
                                        fit_Gibbs1$param$tau[r,])) /sd_y1
      pred_Gibbs2[i,] <- c(LSBP_density(y_grid2[i],X1_2,X2_2,
                                        fit_Gibbs2$param$beta_mixing[r,,],
                                        fit_Gibbs2$param$beta_kernel[r,,],
                                        fit_Gibbs2$param$tau[r,])) /sd_y2
      pred_Gibbs3[i,] <- c(LSBP_density(y_grid3[i],X1_3,X2_3,
                                        fit_Gibbs3$param$beta_mixing[r,,],
                                        fit_Gibbs3$param$beta_kernel[r,,],
                                        fit_Gibbs3$param$tau[r,])) /sd_y3
      
    }
    for(j in 1:dim(X1_1)[1]){
      #j=3
      dens1=pred_Gibbs1[,j];dens2=pred_Gibbs2[,j];dens3=pred_Gibbs3[,j]
      dmax=pmax(dens1,dens2,dens3)
      inte_unl=simpson_inte(y_grid, dmax)
      Est_UNL[r,j]=inte_unl
    }
    
  }
  return(Est_UNL)
}


simulation_unl_con_best<-function(case,nlist,nrep=100){
  n1=nlist$n1;n2=nlist$n2;n3=nlist$n3
  #case=1;n1=200;n2=200;n3=200;nrep=100
  set.seed(25)
  data_con=gen_data_con(case = case,n1=n1,n2=n2,n3=n3,nrep=nrep)
  #n         <- nrow(data_con) # Number of observations
  p         <- 2         # Row and colums of the design for the kernel
  p_splines_bs <- 4         # Number of splines components 
  
  R         <- 5000     # Number of replications
  burn_in   <- 2000      # Burn-in period
  H         <- 20        # Number of mixture components
  
  prior_bsn1 <- prior_LSBP(p,p,
                           b_mixing = rep(0,p), B_mixing=diag(10,p), 
                           b_kernel = rep(0,p_splines_bs+1), B_kernel=diag(10,p_splines_bs+1), 
                           a_tau = 2, b_tau= 0.5)
  
  prior_bsn0 <- prior_LSBP(p,p,
                           b_mixing = rep(0,p), B_mixing=diag(10,p), 
                           b_kernel = rep(0,4), B_kernel=diag(10,4), 
                           a_tau = 2, b_tau= 0.5)
  
  prior <- prior_LSBP(p,p,
                      b_mixing = rep(0,p), B_mixing=diag(10,p), 
                      b_kernel = rep(0,p), B_kernel=diag(10,p), 
                      a_tau = 2, b_tau= 0.5)
  
  xlength=30
  #Est_OVL=array(0,dim=c(nrep,R,xlength))
  
  if(case==1){
    model_formula_y1 <- Formula::as.Formula(y ~ x | x)
    model_formula_y2 <- Formula::as.Formula(y ~ x | x)
    model_formula_y3 <- Formula::as.Formula(y ~ x | x)
    prior1=prior;prior2=prior;prior3=prior
  }else if(case==2){
    model_formula_y1 <- Formula::as.Formula(y ~ BS.1+BS.2+BS.3+BS.4 | x)
    model_formula_y2 <- Formula::as.Formula(y ~ x | x)
    model_formula_y3 <- Formula::as.Formula(y ~ x | x)
    prior1=prior_bsn1;prior2=prior;prior3=prior
  }else if(case==3){
    model_formula_y1 <- Formula::as.Formula(y ~ BS.1+BS.2+BS.3 | x)
    model_formula_y2 <- Formula::as.Formula(y ~ x | x)
    model_formula_y3 <- Formula::as.Formula(y ~ x | x)
    prior1=prior_bsn0;prior2=prior;prior3=prior
  }
  
  nrep_list <- list()
  for (i in 1:nrep) {
    nrep_list[[as.character(i)]] <- i
  }
  
  Est_UNL=future_lapply(X=nrep_list,FUN=fit_and_cal_unl_con_best,data_con=data_con,R=R,
                        burn_in=burn_in,prior1=prior1,prior2=prior2,prior3=prior3,xlength=xlength,H=H,
                        model_formula_y1=model_formula_y1,model_formula_y2=model_formula_y2,
                        model_formula_y3=model_formula_y3,case=case,future.seed = TRUE)
  
  return(Est_UNL)
}

tranform_list_simu_con<-function(simulation_object,R=5000,xlength=30){
  nrep=length(names(simulation_object))
  Est_UNL=array(0,dim=c(nrep,R,xlength))
  counter=1
  for (i in names(simulation_object)){
    Est_UNL[counter,,]=simulation_object[[i]]
    counter=counter+1
  }
  return(Est_UNL)
}



