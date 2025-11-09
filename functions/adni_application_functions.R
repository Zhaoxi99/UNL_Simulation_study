rm(list=ls()) 
library(LSBP)    # Load the LSBP package
library(ggplot2) # Graphical library
library(coda)    # For MCMC analysis
library(splines) # For computing the natural B-splines basis
library(pracma)
library(DIRECT)


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
fit_and_cal_inf_con_bs<-function(data_con1,data_con2,data_con3,response,R,burn_in,prior,H,
                                 model_formula_y1,model_formula_y2,model_formula_y3,nknots=0){
  #data_con1=adni_1_fitting;data_con2=adni_2_fitting;data_con3=adni_3_fitting;response="HCI"
  #prior=prior_bsn0
  #model_formula_y1=model_formula_HCI_bsn0;model_formula_y2=model_formula_HCI_bsn0;model_formula_y3=model_formula_HCI_bsn0
  data_fitting1=na.omit(data.frame(y=data_con1[[response]],age=data_con1$age,gender=data_con1$gender))
  data_fitting2=na.omit(data.frame(y=data_con2[[response]],age=data_con2$age,gender=data_con2$gender))
  data_fitting3=na.omit(data.frame(y=data_con3[[response]],age=data_con3$age,gender=data_con3$gender))
  
  #nknots=0
  if(nknots==0){
    age_bs1 <- bs(data_fitting1$age, intercept = FALSE, degree = 3, knots = c())
    age_bs2 <- bs(data_fitting2$age, intercept = FALSE, degree = 3, knots = c())
    age_bs3 <- bs(data_fitting3$age, intercept = FALSE, degree = 3, knots = c()) 
  }else if(nknots==1){
    age_bs1 <- bs(data_fitting1$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting1$age,probs = c(0.5)))
    age_bs2 <- bs(data_fitting2$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting2$age,probs = c(0.5)))
    age_bs3 <- bs(data_fitting3$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting3$age,probs = c(0.5))) 
  }else if(nknots==2){
    age_bs1 <- bs(data_fitting1$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting1$age,probs = c(0.33,0.66)))
    age_bs2 <- bs(data_fitting2$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting2$age,probs = c(0.33,0.66)))
    age_bs3 <- bs(data_fitting3$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting3$age,probs = c(0.33,0.66))) 
  }else if(nknots==3){
    age_bs1 <- bs(data_fitting1$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting1$age,probs = c(0.25,0.5,0.75)))
    age_bs2 <- bs(data_fitting2$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting2$age,probs = c(0.25,0.5,0.75)))
    age_bs3 <- bs(data_fitting3$age, intercept = FALSE, degree = 3, knots = quantile(data_fitting3$age,probs = c(0.25,0.5,0.75))) 
  }
  
  data_fitting1=data.frame(data_fitting1,BS=age_bs1)
  data_fitting2=data.frame(data_fitting2,BS=age_bs2)
  data_fitting3=data.frame(data_fitting3,BS=age_bs3)
  
  #y_merge=c(data_con$y1[,l],data_con$y2[,l],data_con$y3[,l])
  data_scaled1=data.frame(scale(data_fitting1[,-3]),gender=data_fitting1$gender)
  data_scaled2=data.frame(scale(data_fitting2[,-3]),gender=data_fitting2$gender)
  data_scaled3=data.frame(scale(data_fitting3[,-3]),gender=data_fitting3$gender)
  
  # Gibbs algorithm
  set.seed(23) # The seed is setted so that the Gibbs sampler is reproducible.
  fit_Gibbs1   <- LSBP_Gibbs(model_formula_y1, data=data_scaled1, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs2   <- LSBP_Gibbs(model_formula_y2, data=data_scaled2, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs3   <- LSBP_Gibbs(model_formula_y3, data=data_scaled3, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  
  ##the design matrix for computing information criteria
  X1_1_inf=as.matrix(cbind(1,data_scaled1[,3:(dim(data_scaled1)[2]-1)],ifelse(data_scaled1$gender==2,1,0)))
  X2_1_inf=cbind(1,data_scaled1$age,ifelse(data_scaled1$gender==2,1,0))
  
  X1_2_inf=as.matrix(cbind(1,data_scaled2[,3:(dim(data_scaled2)[2]-1)],ifelse(data_scaled2$gender==2,1,0)))
  X2_2_inf=cbind(1,data_scaled2$age,ifelse(data_scaled2$gender==2,1,0))
  
  X1_3_inf=as.matrix(cbind(1,data_scaled3[,3:(dim(data_scaled3)[2]-1)],ifelse(data_scaled3$gender==2,1,0)))
  X2_3_inf=cbind(1,data_scaled3$age,ifelse(data_scaled3$gender==2,1,0))
  
  sd_y1=sd(data_fitting1$y)
  sd_y2=sd(data_fitting2$y)
  sd_y3=sd(data_fitting3$y)
  
  dens=cal_density(fit_Gibbs1=fit_Gibbs1,fit_Gibbs2=fit_Gibbs2,fit_Gibbs3=fit_Gibbs3,
                   R=R,X1_1=X1_1_inf,X2_1=X2_1_inf,X1_2=X1_2_inf,X2_2=X2_2_inf,
                   X1_3=X1_3_inf,X2_3=X2_3_inf,
                   y_grid1=data_scaled1$y,y_grid2=data_scaled2$y,y_grid3=data_scaled3$y,
                   sd_y1=sd_y1,sd_y2=sd_y2,sd_y3=sd_y3)
  waic_1=waic_density(dens$pred_Gibbs1);waic_2=waic_density(dens$pred_Gibbs2)
  waic_3=waic_density(dens$pred_Gibbs3)
  lpml_1=lpml_density(dens$pred_Gibbs1);lpml_2=lpml_density(dens$pred_Gibbs2)
  lpml_3=lpml_density(dens$pred_Gibbs3)
  return(list(fit_Gibbs1=fit_Gibbs1,fit_Gibbs2=fit_Gibbs2,fit_Gibbs3=fit_Gibbs3,
              waic_1=waic_1,waic_2=waic_2,waic_3=waic_3,
              lpml_1=lpml_1,lpml_2=lpml_2,lpml_3=lpml_3))
}

fit_and_cal_inf_con<-function(data_con1,data_con2,data_con3,response,R,burn_in,prior,H,
                              model_formula_y1,model_formula_y2,model_formula_y3){
  #data_con1=adni_1_fitting;data_con2=adni_2_fitting;data_con3=adni_3_fitting;response="HCI"
  #prior=prior
  #model_formula_y1=model_formula_HCI;model_formula_y2=model_formula_HCI;model_formula_y3=model_formula_HCI
  data_fitting1=na.omit(data.frame(y=data_con1[[response]],age=data_con1$age,gender=data_con1$gender))
  data_fitting2=na.omit(data.frame(y=data_con2[[response]],age=data_con2$age,gender=data_con2$gender))
  data_fitting3=na.omit(data.frame(y=data_con3[[response]],age=data_con3$age,gender=data_con3$gender))
  
  
  #y_merge=c(data_con$y1[,l],data_con$y2[,l],data_con$y3[,l])
  data_scaled1=data.frame(scale(data_fitting1[,-3]),gender=data_fitting1$gender)
  data_scaled2=data.frame(scale(data_fitting2[,-3]),gender=data_fitting2$gender)
  data_scaled3=data.frame(scale(data_fitting3[,-3]),gender=data_fitting3$gender)
  
  # Gibbs algorithm
  set.seed(23) # The seed is setted so that the Gibbs sampler is reproducible.
  fit_Gibbs1   <- LSBP_Gibbs(model_formula_y1, data=data_scaled1, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs2   <- LSBP_Gibbs(model_formula_y2, data=data_scaled2, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  set.seed(23)
  fit_Gibbs3   <- LSBP_Gibbs(model_formula_y3, data=data_scaled3, H=H, prior=prior, 
                             control=control_Gibbs(R=R,burn_in=burn_in,method_init="random"),
                             verbose=FALSE)
  
  ##the design matrix for computing information criteria
  X1_1_inf=cbind(1,data_scaled1$age,ifelse(data_scaled1$gender==2,1,0))
  X2_1_inf=X1_1_inf
  
  X1_2_inf=cbind(1,data_scaled2$age,ifelse(data_scaled2$gender==2,1,0))
  X2_2_inf=X1_2_inf
  
  X1_3_inf=cbind(1,data_scaled3$age,ifelse(data_scaled3$gender==2,1,0))
  X2_3_inf=X1_3_inf
  
  sd_y1=sd(data_fitting1$y)
  sd_y2=sd(data_fitting2$y)
  sd_y3=sd(data_fitting3$y)
  
  dens=cal_density(fit_Gibbs1=fit_Gibbs1,fit_Gibbs2=fit_Gibbs2,fit_Gibbs3=fit_Gibbs3,
                   R=R,X1_1=X1_1_inf,X2_1=X2_1_inf,X1_2=X1_2_inf,X2_2=X2_2_inf,
                   X1_3=X1_3_inf,X2_3=X2_3_inf,
                   y_grid1=data_scaled1$y,y_grid2=data_scaled2$y,y_grid3=data_scaled3$y,
                   sd_y1=sd_y1,sd_y2=sd_y2,sd_y3=sd_y3)
  waic_1=waic_density(dens$pred_Gibbs1);waic_2=waic_density(dens$pred_Gibbs2)
  waic_3=waic_density(dens$pred_Gibbs3)
  lpml_1=lpml_density(dens$pred_Gibbs1);lpml_2=lpml_density(dens$pred_Gibbs2)
  lpml_3=lpml_density(dens$pred_Gibbs3)
  return(list(fit_Gibbs1=fit_Gibbs1,fit_Gibbs2=fit_Gibbs2,fit_Gibbs3=fit_Gibbs3,
              waic_1=waic_1,waic_2=waic_2,waic_3=waic_3,
              lpml_1=lpml_1,lpml_2=lpml_2,lpml_3=lpml_3))
}

effect_size_conditional_density_gender<-function(fit_Gibbs1,fit_Gibbs2,fit_Gibbs3,
                                                 R,X1_1,X2_1,X1_2,X2_2,X1_3,X2_3,
                                                 y_grid,y_grid1,y_grid2,y_grid3,
                                                 sd_y1,sd_y2,sd_y3){
  
  # Posterior density - Gibbs sampling
  #print("b")
  pred_Gibbs1 <- matrix(0,length(y_grid),dim(X1_1)[1])
  pred_Gibbs2 <- matrix(0,length(y_grid),dim(X1_2)[1])
  pred_Gibbs3 <- matrix(0,length(y_grid),dim(X1_3)[1])
  
  
  Est_UNL <- matrix(0,R,dim(X1_3)[1])
  #print("a")
  
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    #r=2
    for(i in 1:length(y_grid)){  # Cycle over the GAD grid
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
      inte_unl=simpson_inte(y_grid, pmax(pred_Gibbs1[,j],pred_Gibbs2[,j],pred_Gibbs3[,j]))
      
      Est_UNL[r,j]=inte_unl
    }
    
    print(r)
    
  }
  return(list(Est_UNL_male=Est_UNL[,1],Est_UNL_female=Est_UNL[,2]))
}


cdf_lsbp_origin<-function(y,fit_Gibbs,r,X1,X2){
  #fit_Gibbs=fit_Gibbs1_ptau;X1=X1_1;X2=X2_1
  #X1_orig=X1_1_orig;X2_orig=X2_1_orig
  #sd_y=sd_ptau1;mean_y=mean_ptau1;sd_x=sd_age1;mean_x=mean_age1
  #y=ptau_grid;y=ptau_grid1
  #r=1002
  H=fit_Gibbs$H
  y=y
  x_grid_length=dim(X2)[1]
  # cdfs=array(0,dim=c(R,x_grid_length,length(y)))
  cdf2=matrix(0,nrow=x_grid_length,ncol=length(y))
  weight_mixing=matrix(0,H,x_grid_length)
  eta=fit_Gibbs$param$beta_mixing[r,,]%*%t(X2)
  prob_equal=exp(eta)/(1+exp(eta))
  prob_next=1-prob_equal
  for(counter in 1:H){
    if(counter==1){
      prob_survive=rep(1,dim(eta)[2])
    }
    
    if(counter!=H){
      #counter=2
      weight_mixing[counter,]=prob_survive*prob_equal[counter,]
      prob_survive=prob_survive-weight_mixing[counter,]
    }else{
      weight_mixing[counter,]=1-apply(weight_mixing,MARGIN = 2,sum)
      if(any(weight_mixing[counter,]<0)){
        negative_index=weight_mixing[counter,]<0
        weight_mixing[counter,negative_index]=0
      }
    }
    
  }
  
  means_lsbp2=fit_Gibbs$param$beta_kernel[r,,]%*%t(X1)
  for(counter2 in 1:x_grid_length){
    #counter2=30
    aux2 <- norMix(mu = means_lsbp2[,counter2], sigma = (sqrt(1/fit_Gibbs$param$tau[r,])), 
                   w = weight_mixing[,counter2])
    #qnorMix(p=0.99, obj = aux1)
    cdf2[counter2,] <- pnorMix(q=y, obj = aux2)
  }
  
  #print(r)
  return(cdf2)
}


effect_size_conditional_density<-function(fit_Gibbs1,fit_Gibbs2,fit_Gibbs3,
                                          R,slct_index,
                                          X1_1,X2_1,X1_2,X2_2,X1_3,X2_3,
                                          y_grid,y_grid1,y_grid2,y_grid3,
                                          sd_y1,sd_y2,sd_y3){
  
  # Posterior density - Gibbs sampling
  #print("b")
  pred_Gibbs1 <- matrix(0,length(y_grid),dim(X1_1)[1])
  pred_Gibbs2 <- matrix(0,length(y_grid),dim(X1_2)[1])
  pred_Gibbs3 <- matrix(0,length(y_grid),dim(X1_3)[1])
  
  pred_Gibbs1_uncon <- rep(0,length(y_grid))
  pred_Gibbs2_uncon <- rep(0,length(y_grid))
  pred_Gibbs3_uncon <- rep(0,length(y_grid))
  
  pred_Gibbs1_slct <- array(0,dim=c(R,length(y_grid),length(slct_index)))
  pred_Gibbs2_slct <- array(0,dim=c(R,length(y_grid),length(slct_index)))
  pred_Gibbs3_slct <- array(0,dim=c(R,length(y_grid),length(slct_index)))
  
  Est_UNL <- matrix(0,R,dim(X1_3)[1])
  #Est_OVL_uncon <- rep(0,R)
  #Est_OVL_2class_a_uncon <- rep(0,R)
  #Est_OVL_2class_b_uncon <- rep(0,R)
  #print("a")
  
  for(r in 1:R){      # Cycle over the iterations of the MCMC chain
    #r=2
    for(i in 1:length(y_grid)){  # Cycle over the GAD grid
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
      pred_Gibbs1_slct[r,i,]=pred_Gibbs1[i,slct_index]
      pred_Gibbs2_slct[r,i,]=pred_Gibbs2[i,slct_index]
      pred_Gibbs3_slct[r,i,]=pred_Gibbs3[i,slct_index]
      
    }
    
    cdf1_st_end=cdf_lsbp_origin(y=c(y_grid1[1],y_grid1[length(y_grid1)]),fit_Gibbs=fit_Gibbs1,r=r,X1=X1_1,X2=X2_1)
    cdf2_st_end=cdf_lsbp_origin(y=c(y_grid2[1],y_grid2[length(y_grid2)]),fit_Gibbs=fit_Gibbs2,r=r,X1=X1_2,X2=X2_2)
    cdf3_st_end=cdf_lsbp_origin(y=c(y_grid3[1],y_grid3[length(y_grid3)]),fit_Gibbs=fit_Gibbs3,r=r,X1=X1_3,X2=X2_3)
    for(j in 1:dim(X1_1)[1]){
      #j=3
      inte_unl=simpson_inte(y_grid, pmax(pred_Gibbs1[,j],pred_Gibbs2[,j],pred_Gibbs3[,j]))
      Est_UNL[r,j]=inte_unl
    }
    #print(r)
    
  }
  return(list(pred_Gibbs1_slct=pred_Gibbs1_slct,pred_Gibbs2_slct=pred_Gibbs2_slct,pred_Gibbs3_slct=pred_Gibbs3_slct,
              Est_UNL=Est_UNL))
}