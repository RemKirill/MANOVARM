library(MASS)
library(pracma)
# 
# моделирование данных для дисперсионного анализа повторных наблюдений

Model.ARM<-function(alpha,beta,sigma2,sigma2.1,n.a)
{
  I_<-length(alpha);I_
  T_<-length(beta);T_
  gamma<-matrix(rep(0,T_*I_),ncol=T_);gamma
  n<-sum(n.a)
  Gamma_<-gamma
  # Gamma_<-sapply(seq(T_),function(t)Gamma_[,t]+beta[t])
  # Gamma_<-t(sapply(seq(I_),function(i)Gamma_[i,]+alpha[i]))
  # Gamma_
  
  A<-Reduce('rbind',
            lapply(seq(I_),function(i)matrix(rep(Gamma_[i,],n.a[i]),ncol=T_,byrow=TRUE)))
  A1<-apply(A,2,function(x)x+rnorm(n,0,sqrt(sigma2)))
  group<-do.call('c',apply(cbind(n.a,c(1,2)),1,function(x)rep(x[2],x[1])))
  data.frame(group=group,t(apply(A1,1,function(x)x+rnorm(1,0,sqrt(sigma2.1)))))
}

Mod.missing<-function(data){
  data[c(1),3]<-NA
  data[c(2),4]<-NA
  data[c(3),4]<-NA
  data[nrow(data),3]<-NA
  data[nrow(data),4]<-NA
  data[nrow(data)-1,4]<-NA
  data[nrow(data)-3,3]<-NA
  data[nrow(data)-4,3]<-NA
  data
}

#-------------------------------------
#преобразование данных для матричной формы дисперсионного анализа
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
#################### пример
dat.AR<-Model.ARM(alpha=c(1,-1),beta=c(0,1,-1),
                sigma2=1,sigma2.1=0.5,n.a=c(10,12))

AR.One(dat.AR)
