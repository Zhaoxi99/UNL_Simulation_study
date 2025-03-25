Tutorial for Simulation study in the paper “The Underlap Coefficient as
Measure of a Biomarker’s Discriminatory Ability in a Three-Class Disease
Setting”
================
Zhaoxi Zhang, Vanda Inácio, and Miguel de Carvalho

In this tutorial we describe the steps for obtaining the results of the
simulation study in the unconditional case in Section 5.1 of the paper
"The Underlap Coefficient as Measure of a Biomarker’s Discriminatory
Ability in a Three-Class Disease Setting".

As a preliminary step, we load on a clean environment all the required
libraries.

As a preliminary step, we load on a clean environment all the required
libraries.

``` r
library(LSBP)    # Load the LSBP package
library(ggplot2) # Graphical library
library(coda)    # For MCMC analysis
library(splines) # For computing the natural B-splines basis
library(pracma)
library(DIRECT)
library(sn)
library(nor1mix)
library(future)
library(future.apply)
library(rjags)
```

We need to load some predefined functions the simulation study in the
unconditional case.

``` r
source("../functions/unconditional_simulation_functions.R")
```

Because we conduct the estimation of underlap repeatedly for 100
datasets in each scenario, so the process could be accelerated by
parallelizing the computation. Here I use the future package to
demonstrate one way of parallelizing the code, the utilization of the
future package is incorporated in the simulation_unl_uncon function.

``` r
plan(multisession, workers = detectCores()-2)
nrep=100;nsave=5000;nburn=2000
nlist1=list(n1=100,n2=100,n3=100)
nlist2=list(n1=200,n2=200,n3=200)
nlist3=list(n1=500,n2=500,n3=500)
nlist4=list(n1=1000,n2=1000,n3=1000)
nlist5=list(n1=100,n2=300,n3=500)
##case1 simulation
simu_unconc1_s1=simulation_unl_uncon(case=1,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc1_s2=simulation_unl_uncon(case=1,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc1_s3=simulation_unl_uncon(case=1,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc1_s4=simulation_unl_uncon(case=1,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc1_s5=simulation_unl_uncon(case=1,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)

##case2 simulation
simu_unconc2_s1=simulation_unl_uncon(case=2,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc2_s2=simulation_unl_uncon(case=2,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc2_s3=simulation_unl_uncon(case=2,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc2_s4=simulation_unl_uncon(case=2,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc2_s5=simulation_unl_uncon(case=2,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case3 simulation
simu_unconc3_s1=simulation_unl_uncon(case=3,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc3_s2=simulation_unl_uncon(case=3,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc3_s3=simulation_unl_uncon(case=3,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc3_s4=simulation_unl_uncon(case=3,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc3_s5=simulation_unl_uncon(case=3,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case4 simulation
simu_unconc4_s1=simulation_unl_uncon(case=4,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc4_s2=simulation_unl_uncon(case=4,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc4_s3=simulation_unl_uncon(case=4,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc4_s4=simulation_unl_uncon(case=4,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc4_s5=simulation_unl_uncon(case=4,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case5 simulation
simu_unconc5_s1=simulation_unl_uncon(case=5,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc5_s2=simulation_unl_uncon(case=5,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc5_s3=simulation_unl_uncon(case=5,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc5_s4=simulation_unl_uncon(case=5,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc5_s5=simulation_unl_uncon(case=5,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case6 simulation
simu_unconc6_s1=simulation_unl_uncon(case=6,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc6_s2=simulation_unl_uncon(case=6,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc6_s3=simulation_unl_uncon(case=6,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc6_s4=simulation_unl_uncon(case=6,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc6_s5=simulation_unl_uncon(case=6,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case7 simulation
simu_unconc7_s1=simulation_unl_uncon(case=7,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc7_s2=simulation_unl_uncon(case=7,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc7_s3=simulation_unl_uncon(case=7,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc7_s4=simulation_unl_uncon(case=7,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc7_s5=simulation_unl_uncon(case=7,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case8 simulation
simu_unconc8_s1=simulation_unl_uncon(case=8,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc8_s2=simulation_unl_uncon(case=8,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc8_s3=simulation_unl_uncon(case=8,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc8_s4=simulation_unl_uncon(case=8,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc8_s5=simulation_unl_uncon(case=8,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
##case9 simulation
simu_unconc9_s1=simulation_unl_uncon(case=9,nlist=nlist1,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc9_s2=simulation_unl_uncon(case=9,nlist=nlist2,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc9_s3=simulation_unl_uncon(case=9,nlist=nlist3,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc9_s4=simulation_unl_uncon(case=9,nlist=nlist4,nrep=nrep,nsave=nsave,nburn=nburn)
simu_unconc9_s5=simulation_unl_uncon(case=9,nlist=nlist5,nrep=nrep,nsave=nsave,nburn=nburn)
```

We need to use the function tranform_list_simu to transform the results
from future_lapply to a more interpretable manner.

``` r
simu_unconc1_s1=tranform_list_simu(simu_unconc1_s1,nsave=nsave,nrep = nrep)
simu_unconc1_s2=tranform_list_simu(simu_unconc1_s2,nsave=nsave,nrep = nrep)
simu_unconc1_s3=tranform_list_simu(simu_unconc1_s3,nsave=nsave,nrep = nrep)
simu_unconc1_s4=tranform_list_simu(simu_unconc1_s4,nsave=nsave,nrep = nrep)
simu_unconc1_s5=tranform_list_simu(simu_unconc1_s5,nsave=nsave,nrep = nrep)

simu_unconc2_s1=tranform_list_simu(simu_unconc2_s1,nsave=nsave,nrep = nrep)
simu_unconc2_s2=tranform_list_simu(simu_unconc2_s2,nsave=nsave,nrep = nrep)
simu_unconc2_s3=tranform_list_simu(simu_unconc2_s3,nsave=nsave,nrep = nrep)
simu_unconc2_s4=tranform_list_simu(simu_unconc2_s4,nsave=nsave,nrep = nrep)
simu_unconc2_s5=tranform_list_simu(simu_unconc2_s5,nsave=nsave,nrep = nrep)

simu_unconc3_s1=tranform_list_simu(simu_unconc3_s1,nsave=nsave,nrep = nrep)
simu_unconc3_s2=tranform_list_simu(simu_unconc3_s2,nsave=nsave,nrep = nrep)
simu_unconc3_s3=tranform_list_simu(simu_unconc3_s3,nsave=nsave,nrep = nrep)
simu_unconc3_s4=tranform_list_simu(simu_unconc3_s4,nsave=nsave,nrep = nrep)
simu_unconc3_s5=tranform_list_simu(simu_unconc3_s5,nsave=nsave,nrep = nrep)

simu_unconc4_s1=tranform_list_simu(simu_unconc4_s1,nsave=nsave,nrep = nrep)
simu_unconc4_s2=tranform_list_simu(simu_unconc4_s2,nsave=nsave,nrep = nrep)
simu_unconc4_s3=tranform_list_simu(simu_unconc4_s3,nsave=nsave,nrep = nrep)
simu_unconc4_s4=tranform_list_simu(simu_unconc4_s4,nsave=nsave,nrep = nrep)
simu_unconc4_s5=tranform_list_simu(simu_unconc4_s5,nsave=nsave,nrep = nrep)

simu_unconc5_s1=tranform_list_simu(simu_unconc5_s1,nsave=nsave,nrep = nrep)
simu_unconc5_s2=tranform_list_simu(simu_unconc5_s2,nsave=nsave,nrep = nrep)
simu_unconc5_s3=tranform_list_simu(simu_unconc5_s3,nsave=nsave,nrep = nrep)
simu_unconc5_s4=tranform_list_simu(simu_unconc5_s4,nsave=nsave,nrep = nrep)
simu_unconc5_s5=tranform_list_simu(simu_unconc5_s5,nsave=nsave,nrep = nrep)

simu_unconc6_s1=tranform_list_simu(simu_unconc6_s1,nsave=nsave,nrep = nrep)
simu_unconc6_s2=tranform_list_simu(simu_unconc6_s2,nsave=nsave,nrep = nrep)
simu_unconc6_s3=tranform_list_simu(simu_unconc6_s3,nsave=nsave,nrep = nrep)
simu_unconc6_s4=tranform_list_simu(simu_unconc6_s4,nsave=nsave,nrep = nrep)
simu_unconc6_s5=tranform_list_simu(simu_unconc6_s5,nsave=nsave,nrep = nrep)

simu_unconc7_s1=tranform_list_simu(simu_unconc7_s1,nsave=nsave,nrep = nrep)
simu_unconc7_s2=tranform_list_simu(simu_unconc7_s2,nsave=nsave,nrep = nrep)
simu_unconc7_s3=tranform_list_simu(simu_unconc7_s3,nsave=nsave,nrep = nrep)
simu_unconc7_s4=tranform_list_simu(simu_unconc7_s4,nsave=nsave,nrep = nrep)
simu_unconc7_s5=tranform_list_simu(simu_unconc7_s5,nsave=nsave,nrep = nrep)

simu_unconc8_s1=tranform_list_simu(simu_unconc8_s1,nsave=nsave,nrep = nrep)
simu_unconc8_s2=tranform_list_simu(simu_unconc8_s2,nsave=nsave,nrep = nrep)
simu_unconc8_s3=tranform_list_simu(simu_unconc8_s3,nsave=nsave,nrep = nrep)
simu_unconc8_s4=tranform_list_simu(simu_unconc8_s4,nsave=nsave,nrep = nrep)
simu_unconc8_s5=tranform_list_simu(simu_unconc8_s5,nsave=nsave,nrep = nrep)

simu_unconc9_s1=tranform_list_simu(simu_unconc9_s1,nsave=nsave,nrep = nrep)
simu_unconc9_s2=tranform_list_simu(simu_unconc9_s2,nsave=nsave,nrep = nrep)
simu_unconc9_s3=tranform_list_simu(simu_unconc9_s3,nsave=nsave,nrep = nrep)
simu_unconc9_s4=tranform_list_simu(simu_unconc9_s4,nsave=nsave,nrep = nrep)
simu_unconc9_s5=tranform_list_simu(simu_unconc9_s5,nsave=nsave,nrep = nrep)
```

We can then calculate the posterior medians, the 95% credible intervals,
the width of the 95% credible intervals, and the coverage of the 95%
credible intervals of UNL estimates in each scenario.

``` r
for(case_i in 1:9){
  median_name=paste("median_uncon","c",case_i,sep="")
  simu_name1=paste("simu_uncon","c",case_i,"_s1",sep="")
  simu_name2=paste("simu_uncon","c",case_i,"_s2",sep="")
  simu_name3=paste("simu_uncon","c",case_i,"_s3",sep="")
  simu_name4=paste("simu_uncon","c",case_i,"_s4",sep="")
  simu_name5=paste("simu_uncon","c",case_i,"_s5",sep="")
  assign(median_name,data.frame(UNL=c(apply(get(simu_name1)[["Est_UNL1"]],FUN=median,MARGIN = 2),
                                    apply(get(simu_name2)[["Est_UNL1"]],FUN=median,MARGIN = 2),
                                    apply(get(simu_name3)[["Est_UNL1"]],FUN=median,MARGIN = 2),
                                    apply(get(simu_name4)[["Est_UNL1"]],FUN=median,MARGIN = 2),
                                    apply(get(simu_name5)[["Est_UNL1"]],FUN=median,MARGIN = 2)),
                              size=as.factor(c(rep("s1",nrep),rep("s2",nrep),rep("s3",nrep),
                                               rep("s4",nrep),rep("s5",nrep)))))
  
}


for(case_i in 1:9){
  width_name=paste("width_uncon","c",case_i,sep="")
  simu_name1=paste("simu_uncon","c",case_i,"_s1",sep="")
  simu_name2=paste("simu_uncon","c",case_i,"_s2",sep="")
  simu_name3=paste("simu_uncon","c",case_i,"_s3",sep="")
  simu_name4=paste("simu_uncon","c",case_i,"_s4",sep="")
  simu_name5=paste("simu_uncon","c",case_i,"_s5",sep="")
  assign(width_name,data.frame(width=c(apply(get(simu_name1)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.975)-
                                         apply(get(simu_name1)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.025),
                                       apply(get(simu_name2)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.975)-
                                         apply(get(simu_name2)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.025),
                                       apply(get(simu_name3)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.975)-
                                         apply(get(simu_name3)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.025),
                                       apply(get(simu_name4)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.975)-
                                         apply(get(simu_name4)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.025),
                                       apply(get(simu_name5)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.975)-
                                         apply(get(simu_name5)[["Est_UNL1"]],FUN=quantile,MARGIN = 2,probs=0.025)),
                               size=as.factor(c(rep("s1",nrep),rep("s2",nrep),rep("s3",nrep),
                                                rep("s4",nrep),rep("s5",nrep)))))
  

for(case_i in 1:9){
  for(size_i in 1:5){
    CI_name=paste("CI95_uncon","c",case_i,"_s",size_i,sep="")
    simu_name=paste("simu_uncon","c",case_i,"_s",size_i,sep="")
    inte_name="Est_UNL1"
    assign(CI_name,data.frame(lower=apply(get(simu_name)[[inte_name]],FUN=quantile,MARGIN = 2,probs=0.025),
                              higher=apply(get(simu_name)[[inte_name]],FUN=quantile,MARGIN = 2,probs=0.975)))
    
  }
}

coverage=matrix(0,nrow=9,ncol = 5)
for(case_i in 1:9){
  for(size_i in 1:5){
    CI_name=paste("CI95_uncon","c",case_i,"_s",size_i,sep="")
    coverage[case_i,size_i]=sum(true_unl_uncon(case=case_i)>=get(CI_name)[["lower"]]&
                                  true_unl_uncon(case=case_i)<=get(CI_name)[["higher"]])
  }
}
  
}
```

## Posterior median boxplots

We can make boxplots of the posterior median of the coefficient of
underlap across 100 simulated datasets for different parameter
configurations and sample sizes as following.

``` r
theme_set(theme_bw())
median_box_unconc1=ggplot(median_unconc1, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=1), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario I - UNL = ",2.792,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc2=ggplot(median_unconc2, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=2), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario I - UNL = ",1.919,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc3=ggplot(median_unconc3, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=3), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario I - UNL = ",1.139,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc4=ggplot(median_unconc4, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=4), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario II - UNL = ",2.527,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc5=ggplot(median_unconc5, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=5), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario II - UNL = ",1.855,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc6=ggplot(median_unconc6, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=6), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario II - UNL = ",1.191,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc7=ggplot(median_unconc7, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=7), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario III - UNL = ",2.508,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc8=ggplot(median_unconc8, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=8), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario III - UNL = ",1.933,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

median_box_unconc9=ggplot(median_unconc9, aes(x=size, y=UNL)) + 
  geom_boxplot()+geom_hline(yintercept = true_unl_uncon(case=9), linetype = "solid", color = "red")+
  ggtitle(paste("Scenario III - UNL = ",1.143,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(1,3,by=0.2),limits = c(1,3))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )
```

``` r
cowplot::plot_grid(median_box_unconc1,median_box_unconc2,median_box_unconc3,
                   median_box_unconc4,median_box_unconc5,median_box_unconc6,
                   median_box_unconc7,median_box_unconc8,median_box_unconc9,ncol = 3)
```

![](UNL_simu_tutorial_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->

## 95% CI width boxplots

We can make the boxplot of the width of the 95% credible intervals of
the coefficient of underlap across 100 simulated datasets for varying
parameter configurations and sample sizes as following.

``` r
##width box plots
theme_set(theme_bw())
width_box_unconc1=ggplot() + geom_boxplot(data=width_unconc1, aes(x=size, y=width))+
  ggtitle(paste("Scenario I - UNL = ",2.792,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc2=ggplot() + geom_boxplot(data=width_unconc2, aes(x=size, y=width))+
  ggtitle(paste("Scenario I - UNL = ",1.919,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc3=ggplot() + geom_boxplot(data=width_unconc3, aes(x=size, y=width))+
  ggtitle(paste("Scenario I - UNL = ",1.139,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc4=ggplot(width_unconc4, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario II - UNL = ",2.527,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc5=ggplot(width_unconc5, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario II - UNL = ",1.855,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc6=ggplot(width_unconc6, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario II - UNL = ",1.191,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc7=ggplot(width_unconc7, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario III - UNL = ",2.508,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc8=ggplot(width_unconc8, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario III - UNL = ",1.933,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )

width_box_unconc9=ggplot(width_unconc9, aes(x=size, y=width)) + 
  geom_boxplot()+
  ggtitle(paste("Scenario III - UNL = ",1.143,sep=""))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_continuous(breaks =  seq(0,0.5,by=0.1),limits = c(0,0.4))+
  scale_x_discrete(labels=c("(100,100,100)","(200,200,200)","(500,500,500)","(1000,1000,1000)",
                            "(100,300,500)"))+xlab("Sample Size")+ylab("Interval width")+
  theme(
    plot.title = element_text(size = 10,face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 3)
  )
```

``` r
cowplot::plot_grid(width_box_unconc1,width_box_unconc2,width_box_unconc3,
                   width_box_unconc4,width_box_unconc5,width_box_unconc6,
                   width_box_unconc7,width_box_unconc8,width_box_unconc9,ncol = 3)
```

![](UNL_simu_tutorial_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->
