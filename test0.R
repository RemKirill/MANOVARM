
#матрица R

source("/Users/bulochka/RProject/funcsAR_8.R")
library(MASS)
library('pracma')
# 
# моделирование данных для дисперсионного анализа повторных наблюдений

Model.ARM<-function(alpha,beta,sigma2,sigma2.1,n.a)
{
  I_<-length(alpha);I_
  T_<-length(beta);T_
  gamma<-matrix(rep(0,T_*I_),ncol=T_);gamma
  n<-sum(n.a)
  Gamma_<-gamma
  Gamma_<-sapply(seq(T_),function(t)Gamma_[,t]+beta[t])
  Gamma_<-t(sapply(seq(I_),function(i)Gamma_[i,]+alpha[i]))
  Gamma_

  A<-Reduce('rbind',
         lapply(seq(I_),function(i)matrix(rep(Gamma_[i,],n.a[i]),ncol=T_,byrow=TRUE)))
  A1<-apply(A,2,function(x)x+rnorm(n,0,sqrt(sigma2)))
  group<-do.call('c',apply(cbind(n.a,c(1,2)),1,function(x)rep(x[2],x[1])))
  data.frame(group=group,t(apply(A1,1,function(x)x+rnorm(1,0,sqrt(sigma2.1)))))
}


#-----------------
# пример работы программы на полных данных
alpha<-c(1,2);beta<-c(0,0,0); n.a<-c(5,7);
sigma2<-1;sigma2.1<-4
r<-length(alpha);n<-sum(n.a);T_<-length(beta)

#data<-Model.ARM(alpha,beta,sigma2,sigma2.1,n.a);

View(data)

data <- read.csv("/Users/bulochka/RProject/nir_data.csv", header = TRUE, as.is = FALSE, sep = ';')[-1]
data <- rbind(subset(data, job == 0), subset(data, job == 1), subset(data, job == 2))[c(1,2,3,4,5)]
head(data)
data$job <- data$job + 1
data <- data[!is.na(data[,2]), ]
head(data)
str(data)
names(data)[names(data) == 'job'] <- 'group'
Group<-TRUE
eps<-1/10^5

matrices <- matrices.list(as.matrix(data), data)
# 
ep <- Estimation.Parameters(matrices,Group)

B<-Bias(matrices,Group)
ii<-c(1,1)

# MC<-Mick.Computation(c(t=1,tau=2),B,ii=ii,Group);
# 
# Cov.F(ii,B,Group)
CV<-Covariance.matrix(B,Group)[[2]]
# 
#View(round(CV,2))
p.values.AU<-AU.2(data,eps,Group)$P.values
#eigen(CV)[[1]]

#-------------------------------------
#Стандартная программа
AR.One<-function(dat.AR)
{
  m<-ncol(dat.AR)-1;m
  dat.AR<-na.omit(dat.AR)
  dat.AR.T<-data.frame( stack(dat.AR[,-1]),
                        sub=as.factor(rep(seq(nrow(dat.AR)),m)), 
                        gr=as.factor(rep(dat.AR$group,m)))
  dat.AR.T$values<-as.numeric(dat.AR.T$values)
  
  formula<-values~gr*ind+Error(sub/ind)
  
  aov.out <- aov(formula, data=dat.AR.T)
  aov.S<-summary(aov.out)
  
  pp<-c(pp1<-aov.S[[1]][[1]][1,5],aov.S[[2]][[1]][seq(2),5])
  
  wm<-with(dat.AR.T, tapply(values, list(gr,ind), mean))
  wsd<-with(dat.AR.T, tapply(values, list(gr,ind), sd))
  wn<-with(dat.AR.T, tapply(values, list(gr,ind), length))
  
  wma<-with(dat.AR.T, tapply(values, gr, mean))
  wsda<-with(dat.AR.T, tapply(values, gr, sd))
  wna<-with(dat.AR.T, tapply(values, gr, length))
  
  Means<-list(wm=wm,wsd=wsd,wn=wn)
  MeansA<-list(wmA=wma,wsdA=wsda,wnA=wna)
  names(pp)<-c("p.gr","p.time","p.gr.time")
  list(pp=pp,Means=Means,MeansA=MeansA,dat.AR.T=dat.AR.T)
}
#----------

#-----------Сравнение со стандартной процедурой в R

cbind(AR.One(data)[[1]], 
      AU.2(data,eps,Group)$P.values[-1])
#---------------------------- 

data <- read.csv("/Users/bulochka/RProject/nir_data.csv", header = TRUE, as.is = FALSE, sep = ';')[-1]
data <- rbind(subset(data, Sex == 0), subset(data, Sex == 1))
data$Sex <- data$Sex + 1
names(data)[1] <- 'group'
df1 <- data[,c(1,2,3,4,5)]
df1 <- df1[!is.na(df1[,2]), ]
df2 <- data[,c(1,6,7,8,9)]
df2 <- df2[!is.na(df2[,2]), ]
df3 <- data[,c(1,10,11,12,13)]
df3 <- df3[!is.na(df3[,2]), ]
AU.2(df1,eps,Group)$P.values[-1]
AU.2(df2,eps,Group)$P.values[-1]
AU.2(df3,eps,Group)$P.values[-1]
data <- data[!is.na(data[,2]), ]

opt1 <- optim(c(0.5,0.5,0.3), y1)
opt1

opt2 <- optim(c(0.5,0.5,0.3), y2)
opt2

opt3 <- optim(c(0.5,0.5,0.3), y3)
opt3

y1 <- function(x){
  Group<-TRUE
  eps<-1/10^5
  n_1 <- x[1]*data$Lad.po.1 + x[2]*data$LVd.po.1 + x[3]*data$LVs.po.1
  n_5 <- x[1]*data$Lad.po5 + x[2]*data$LVd.po5 + x[3]*data$LVs.po5
  n_9 <- x[1]*data$Lad.po9 + x[2]*data$LVd.po9 + x[3]*data$LVs.po9
  n_13 <- x[1]*data$Lad.po13 + x[2]*data$LVd.po13 + x[3]*data$LVs.po13
  new_ind <- data.frame(cbind(data$group, n_1, n_5, n_9, n_13))
  names(new_ind)[1] <- 'group'
  new_ind <- new_ind[!is.na(new_ind[,2]), ]
  AU.2(new_ind,eps,Group)$P.values[-1][1]
}

y2 <- function(x){
  Group<-TRUE
  eps<-1/10^5
  n_1 <- x[1]*data$Lad.po.1 + x[2]*data$LVd.po.1 + x[3]*data$LVs.po.1
  n_5 <- x[1]*data$Lad.po5 + x[2]*data$LVd.po5 + x[3]*data$LVs.po5
  n_9 <- x[1]*data$Lad.po9 + x[2]*data$LVd.po9 + x[3]*data$LVs.po9
  n_13 <- x[1]*data$Lad.po13 + x[2]*data$LVd.po13 + x[3]*data$LVs.po13
  new_ind <- data.frame(cbind(data$group, n_1, n_5, n_9, n_13))
  names(new_ind)[1] <- 'group'
  new_ind <- new_ind[!is.na(new_ind[,2]), ]
  AU.2(new_ind,eps,Group)$P.values[-1][2]
}

y3 <- function(x){
  Group<-TRUE
  eps<-1/10^5
  n_1 <- x[1]*data$Lad.po.1 + x[2]*data$LVd.po.1 + x[3]*data$LVs.po.1
  n_5 <- x[1]*data$Lad.po5 + x[2]*data$LVd.po5 + x[3]*data$LVs.po5
  n_9 <- x[1]*data$Lad.po9 + x[2]*data$LVd.po9 + x[3]*data$LVs.po9
  n_13 <- x[1]*data$Lad.po13 + x[2]*data$LVd.po13 + x[3]*data$LVs.po13
  new_ind <- data.frame(cbind(data$group, n_1, n_5, n_9, n_13))
  names(new_ind)[1] <- 'group'
  new_ind <- new_ind[!is.na(new_ind[,2]), ]
  AU.2(new_ind,eps,Group)$P.values[-1][3]
}

