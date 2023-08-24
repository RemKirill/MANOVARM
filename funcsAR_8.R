library("magic")
library("corpcor")

Local.List<-function(X){
  J<-t(apply(X,1,function(x)ifelse( is.na(x),0,1)));J
  
  Nj<-apply(X,1,function(y)length(na.omit(y)));Nj
  Mt<-apply(X,2,function(y)length(na.omit(y)));Mt
  Counts<-list(Nj=Nj,Mt=Mt)
  
  Vi<-matrix(apply(X,1,function(x)mean(x,na.rm=TRUE)),ncol=1) 
  Ui<-matrix(apply(X,2,function(x)mean(x,na.rm=TRUE)),ncol=1) 
  UV<-list(U=Ui,V=Vi)
  
  Lambda_nu<-diag(1/Nj);Lambda_nu
  Lambda_iT<-diag(1/Mt);Lambda_iT
  
  Ri<- Lambda_nu %*% J;Ri
  Pi<- Ri %*% Lambda_iT  %*% t(J); Pi
  
  pi_<-Nj/sum(Nj);pi_
  Pi_<-matrix(rep(pi_,length(pi_)),
              byrow=TRUE,ncol=length(pi_))
  I_<-diag(rep(1,length(pi_)))
  
  Qi<-ginv( I_- Pi+Pi_ );Qi
  
  det(I_- Pi+Pi_)
  
  Ci<-Pi_+Qi;Ci
  Di<-Qi %*% Ri;Di
  
  Hi<-(I_-Ci) %*% Vi + Di %*% Ui ;Hi
  
  Zi<- Vi-Hi
  
  list(Q=Qi,C=Ci,D=Di,H=Hi,Counts=Counts,UV=UV,J=J,
       Lambda=list(Lambda_nu=Lambda_nu,Lambda_iT=Lambda_iT))
}

Group.List<-function(matrices)
{
  # групповая поправка
  
  MM<-Reduce('rbind',lapply(matrices,function(x)
    apply(x,2,function(y)length(na.omit(y)))));MM
  
  M<-t(apply(MM,1,function(x)x/sum(x)));M
  N<-t(apply(t(MM),1,function(x)x/sum(x)));N
  
  P<-M%*%N;P
  
  P.inf<-as.vector(rowSums(MM)/sum(rowSums(MM)));
  P.inf<-matrix(rep(P.inf,nrow(M)),ncol=nrow(M),byrow=TRUE)
  
  
  Q<-solve(eye(nrow(M))-P+P.inf);Q
  A<-P.inf+Q-eye(nrow(M));A
  B_<-Q %*% M; B_
  
  
  XT<-apply(Reduce('rbind',matrices),2,function(x)mean(x,na.rm=TRUE));XT
  
  XI<-unlist(lapply(matrices,function(x)mean(as.vector(x),na.rm=TRUE)));XI
  
  XIT<-Reduce('rbind',lapply(matrices,function(x)apply(x,2,function(y)mean(y,na.rm=TRUE))));XIT
  
  XX<-sum(XI*rowSums(MM))/sum(sum(MM))  
  
  Means<-list(XT=XT,XI=XI,XIT=XIT,XX=XX)
  
  L<-matrix(XT-XX,ncol=1);L
  K<-matrix(XI-XX,ncol=1);K
  G<- -A %*% K + B_ %*% L; G
  list(MM=MM,G=G,AB=list(A=A,B=B_),Means=Means)
}

Bias<-function(matrices,Group)
{
  First<-lapply(matrices,function(mx)Local.List(mx))
  Second<-Group.List(matrices)
  names(Second)
  
  H<-lapply(First,function(x) x$H);H

  V<-lapply(First,function(x)x$UV$V);V;
  
  G<-rep(0,length(matrices))
  if(Group==TRUE) G<-Second$G
 
  
  Z<-list(length(matrices))
  for(i in seq(length(matrices))) Z[[i]]<-V[[i]]-H[[i]]-G[i]
  
  matrices_<-matrices
  for(i in seq(length(matrices)))
  {
    mx<-matrices[[i]]
    mx_<-apply(cbind(mx,Z[[i]]),1,function(x) x[-length(x)]-x[length(x)])
     
    mx_<-t(mx_);mx_
    matrices_[[i]]<-mx_
  }
  
  list(matrices_=matrices_,First=First,Second=Second,Z=Z)
}

LF.Comp<-function(B)
{
  Lf<-lapply(B$matrices_,function(X)
  na.omit(cbind(stack(data.frame(X)),
                j=rep(seq(nrow(X)),ncol(X)))))
Lf
Y<-Reduce('rbind',Lf)[,"values"]
c<-unlist(lapply(lapply(B$First,function(x)x$Counts$Nj),function(y)sum(y)))
cc<-NULL;for(i in seq(length(c)))cc<-c(cc,rep(i,c[i]));cc

Lf<-Reduce('rbind',Lf)

dff<-data.frame(rq=names(table(Lf$ind)))
t<-as.numeric(sapply(Lf$ind,function(x)row.names(dff)[dff[,1]==x]))

Lf<-cbind(Lf,i=cc,t=t)
          
          
list(Lf=Lf,Y=Y)
}

delta.w0<-function(B) 
{
  CC<-lapply(B$First,function(x)x$C);CC
  DD<-lapply(B$First,function(x)x$D)
  
  Nj<-lapply(B$First,function(x)x$Counts$Nj)
  EVVTi<-lapply(Nj,function(x)diag(1/x))
  
  Mt<-lapply(B$First,function(x)x$Counts$Mt)
  EUUTi<-lapply(Mt,function(x)diag(1/x));EUUTi
  
  JJ<-lapply(B$First,function(x)x$J)
  
  
  
  EVUTi<-list();length(EVUTi)<-length(Nj)
  for (i in seq(length(Nj)))
    EVUTi[[i]]<-sapply(Mt[[i]],function(m_it) 
      sapply(Nj[[i]],function(n_ij)1/n_ij/m_it)) * JJ[[i]]
  EVUTi
  
  
  
  Delta<-list();length(Delta)<-length(Nj)
  for (i in seq(length(Nj)))
  {
    
    x1 <-    CC[[i]] %*% EVVTi[[i]] %*% t(CC[[i]])
    x2 <-    CC[[i]] %*% EVUTi[[i]] %*% t(DD[[i]])
    x3<-     t(x2)
    x4 <-    DD[[i]] %*% EUUTi[[i]] %*% t(DD[[i]])
    
    Delta[[i]]<-x1-x2-x3+x4
  }
  
  Delta
}

Mick.Computation<-function(tt,B,ii,Group)
{
  
  Nj<-lapply(B$First,function(x)x$Count$Nj);Nj             
  
  nu<-unlist(lapply(Nj,function(x)length(x)));nu
  C<-B$First[[ii[1]]]$C;D<-B$First[[ii[1]]]$D;C;D
  Mt<-lapply(B$First,function(x)x$Counts$Mt)
  MM<-Reduce('rbind',Mt);MM;colSums(MM)
  JJ<-lapply(B$First,function(x)x$J);JJ

  
  if(ii[1]==ii[2])
    W.A<-sapply(seq(nu[ii[1]]),function(l)
      sapply(seq(nu[ii[2]]),function(j) 
        C[j,l]/Nj[[ii[1]]][l]- D[j,tt[2]]/MM[ii[1],tt[2]] )) else
                                    W.A<-zeros(nu[[ii[1]]],nu[[ii[2]]])
  
  if(Group==TRUE )
  W.B<- B$Second$AB$A[ii[1],ii[2]]/rowSums(MM)[ii[2]]-
       B$Second$AB$B[ii[1],tt[2]]/colSums(MM)[tt[2]] else W.B<-0
  
  W<-W.A+W.B;W
  Cond1<-JJ[[ii[1]]][,tt[1]]==1;Cond2<-JJ[[ii[2]]][,tt[2]]==1
  
  W_<-t(subset(t(subset(W,Cond1)),Cond2));dim(W_);round(W_,2)
  W_
  
}

Covariance.matrix <- function(B,Group)
{
  MM<-Reduce('rbind',lapply(B$First,function(x)x$Counts$Mt));MM
  
  Cov<-Reduce('rbind',lapply(seq(nrow(MM)), function(i1)
    Reduce('cbind',lapply(seq(nrow(MM)),function(i2)Cov.F(c(i1,i2),B,Group))
           )
    )
    )
  
  Cov
  
  
   Lf<-LF.Comp(B)$Lf
   
    list(Y=LF.Comp(B)$Y,Lambda=Cov)
}

Cov.F<-function(ii,B,Group)
{
  i1<-ii[1];i2<-ii[2]
  
  MM<-Reduce('rbind',lapply(B$First,function(x)x$Counts$Mt));MM
  T_<-ncol(B$matrices_[[1]])
  Lambda.1 <- NULL
  
  for(t in seq(T_))
  {
    Lambda.0 <- NULL
    for(tau in seq(T_))
  {
    D.tilda<- delta.tilda(c(t,tau),B,ii=c(i1,i2),Group);
    Mi<-   Mick.Computation(c(t,tau),B,ii=c(i1,i2),Group);Mi
    
    Mi.t<-   Mick.Computation(c(tau,t),B,ii=c(i2,i1),Group);
    Lam<-  eye(nrow(D.tilda),ncol(D.tilda))*ifelse(i1==i2,1,0)*ifelse(t==tau,1,0)-
      Mi-t(Mi.t)+D.tilda;
    Lambda.0 <- cbind(Lambda.0,Lam)
    }
    Lambda.1<-rbind(Lambda.1,Lambda.0)
  }
  Lambda.1
  
}

delta.f<-function(tt,B,ii,Group) 
{
  Nj<-lapply(B$First,function(x)x$Count$Nj);Nj
  nu<-unlist(lapply(Nj,function(x)length(x)));nu
  C<-B$First[[ii[1]]]$C;D<-B$First[[ii[1]]]$D;
  Mt<-lapply(B$First,function(x)x$Counts$Mt)
  MM<-Reduce('rbind',Mt);MM
  JJ<-lapply(B$First,function(x)x$J);JJ
  D.A<-DDD__(ii[1],ii[2],B,Group)
  Cond1<-JJ[[ii[1]]][,tt[1]]==1;Cond2<-JJ[[ii[2]]][,tt[2]]==1
  
 
if(sum(Cond2)==1) resu<-t(Trun(t(D.A),Cond2,Cond1)) else resu<-Trun(D.A,Cond1,Cond2)
  resu
  
}


Trun<-function(D.A,Cond1,Cond2)subset(D.A[,Cond2],Cond1)

delta.tilda<-function(tt,B,ii,Group) 
{
  f_<-delta.f(tt,B,ii,Group) ;f_
  JJ<-lapply(B$First,function(x)x$J);JJ
  Delta<-delta.w0(B);Delta
  Cond2<-JJ[[ii[1]]][,tt[2]]==1;
  Cond1<-JJ[[ii[1]]][,tt[1]]==1
  if(ii[1]==ii[2]) {
    if(sum(Cond2)==1) resu<-t(Trun(Delta [[ii[1]]],Cond2,Cond1)) else resu<-Trun(Delta [[ii[1]]],Cond1,Cond2)
   

    F_<-resu+f_}  else F_<-f_
  F_
}

Matrix.Plan<-function(B)
{
  L<-list(Y=lapply(B$matrices_,function(mx)
    as.vector(na.omit(stack(data.frame(mx))[,1]))),
    m=lapply(B$First,function(x)x$Counts$Mt),
    n=lapply(B$First,function(x)x$Counts$Nj))
  
  matrix.plan (L)
}

Estimation.Parameters <- function(matrices,Group)
{
  B<-Bias(matrices,Group)
  B$First
  
  plan<-Matrix.Plan(B);
  
  CV<-Covariance.matrix(B,Group);Lambda<-CV$Lambda;
  Y<-matrix(CV$Y,ncol=1);Y
  
 
  
  Lambda_reverse <- pseudoinverse(Lambda)
  
  
  b <- t(plan) %*% Lambda_reverse %*% plan
  c <- qr(b)
  
  b_reverse <- qr.coef(c, diag(dim(plan)[2]))
  
 
  teta <- ginv(b) %*% t(plan) %*% Lambda_reverse %*% Y
  
  
  teta
     l <- list("Y" = Y, "teta" = teta, "plan" = plan, 
              "Lambda_r" = Lambda_reverse,B=B,Lambda=Lambda)
}

Estimation.Parameters.Beta <- function(E.P)
{
  MT<-lapply(E.P$B$First,function(x)x$Count$Mt);MT
  T_<-length(MT[[1]]);T_
  I<-length(MT);I
  
  Lambda_reverse <- E.P$Lambda_r
  
  
  plan <- E.P$plan;plan
  
  H_beta <- H.beta.computing(plan, T_, I)
  b <- t(H_beta) %*% Lambda_reverse %*% H_beta
  c <- qr(b)
  solve(b)
  b_reverse <- qr.coef(c, diag(dim(H_beta)[2]))
  
  beta <- b_reverse %*% t(H_beta) %*% Lambda_reverse %*% E.P$Y
  
  beta
}

Estimation.Parameters.Gamma <- function(E.P)
{
  
  I<-length(E.P$B$matrices_) 
  
  T_<-ncol(E.P$B$matrices[[1]])
 
  
  Lambda_reverse <- E.P$Lambda_r
  
  plan <- E.P$plan;plan
  
  H_gamma <- H.Gamma.Computing(plan, I, T_)
  b <- t(H_gamma) %*% Lambda_reverse %*% H_gamma
  c <- qr(b)
  
  b_reverse <- qr.coef(c, diag(dim(H_gamma)[2]))
  
  gamma <- b_reverse %*% t(H_gamma) %*% Lambda_reverse %*% E.P$Y
  
  gamma
}

Beta.Computing <- function(E.P)
{
  teta <- E.P$teta
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  MT<-lapply(E.P$B$First,function(x)x$Count$Mt)
  
  beta <- numeric(I)
  beta <- sapply(1:(T_ - 1), function(t) beta[t] <- teta[1 + (t - 1) * I])
  
  q <- Q.Computing(MT)
  
  beta <- c(beta, sum(q * beta[1:(T_ - 1)]))
  beta
}

Q.Computing <- function(MT)
{
  T_ <- length(MT[[1]])
  I <- length(MT)
  q <- numeric(T_ - 1)
  
  s <- sum(sapply(1:I, function(i) MT[[i]][T_]));s
  
  for(t in 1:(T_ - 1))
  {
    q[t] <- - sum(sapply(1:I, function(i) MT[[i]][t])) / s
  }
  
  q
}

Q_.Computing <- function(MT)
{
  T_ <- length(MT[[1]])
  I <- length(MT)
  q <- matrix(nrow = I - 1, ncol = T_ - 1)
  
  s <- sapply(1:I, function(i) MT[[i]][T_])
  
  for(i in 1:(I - 1))
    for(t in 1:(T_ - 1))
      q[i, t] <- - MT[[i]][t] / s[i]
  
  q
}

Gamma.Computing <- function(E.P)
{
  teta <- E.P$teta
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  MT<-lapply(E.P$B$First,function(x)x$Count$Mt)
  
  gamma <- matrix(nrow = I, ncol = T_)
  
  for(t in 1:(T_ - 1))
    gamma[1:(I - 1), t] <- teta[(2 + (t - 1) * I):(I * t)]
  
  p <- matrix(0, nrow = I - 1, ncol = T_ - 1)
  r <- matrix(0, nrow = I - 1, ncol = T_ - 1)
  
  for(i in 1:(I - 1))
    for(t in 1:(T_ - 1))
    {
      p[i,t] <- -MT[[i]][t] / MT[[I]][t]
      r[i,t] <- MT[[i]][t] / MT[[I]][T_]
    }
  
  gamma[I, 1:(T_ - 1)] <- sapply(1:(T_ - 1), function(t) sum(p[, t] *
                                                               gamma[1:(I - 1), t]))
  
  q <- Q_.Computing(MT)
  
  gamma[1:(I - 1), T_] <- sapply(1:(I - 1), function(i) q[i, ] %*%
                                   gamma[i, 1:(T_ - 1)])
  
  gamma[I, T_] <- sum(sapply(1:(I - 1), function(i) sum(r[i,] *
                                                          gamma[i, 1:(T_ - 1)])))
  row.names(gamma)<-paste("gr",seq(I),sep="_")
  gamma
}

Q.2e.Gamma.Computation <- function(E.P, beta)
{
  
  Y <- E.P$Y
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  
  H_beta <- H.beta.computing(E.P$plan, T_, I)
  Lambda_r <- E.P$Lambda_r
  
  t(Y - H_beta %*% beta) %*% Lambda_r %*% (Y - H_beta %*% beta)
}

Q.2e.Beta.Computation <- function(E.P, gamma)
{
  Y <- E.P$Y
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  
  H_gamma <- H.Gamma.Computing(E.P$plan, I, T_)
  
  Lambda_r <- E.P$Lambda_r
  
  t(Y - H_gamma %*% gamma) %*% Lambda_r %*% (Y - H_gamma %*% gamma)
}

Q.2e.Computation <- function(E.P)
{
  
  Y <- E.P$Y
  H <- E.P$plan
  
  Lambda_r <- E.P$Lambda_r
  
  t(Y - H %*% E.P$teta) %*% Lambda_r %*% (Y - H %*% E.P$teta)
}

F.Gamma <- function(E.P)
{
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  m__<-do.call(sum,lapply(E.P$B$First,function(x)x$Count$Mt))
  n<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
    
  Q_2e <- Q.2e.Computation(E.P)
  df2 <- (m__ - n - I * (T_ - 1))
  df1 <- (I * T_ - T_ - I + 1)
  Q_g<-(Q.2e.Gamma.Computation(E.P, Estimation.Parameters.Beta(E.P)) - Q_2e)
  stat <- (Q.2e.Gamma.Computation(E.P, Estimation.Parameters.Beta(E.P)) - Q_2e) *
    df2 / Q_2e / df1
  
  l <- list("stat" = stat, "df1" = df1, "df2" = df2, Q_2e=Q_2e, Q_g=Q_g)
}

F.Beta <- function(E.P)
{
  
  I<-length(E.P$B$matrices_) 
  T_<-ncol(E.P$B$matrices[[1]])
  m__<-do.call('sum',lapply(E.P$B$First,function(x)x$Count$Mt) )
  
  n<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  
  Q_2e <- Q.2e.Computation(E.P)
  
  df2 <- m__ - n - I * (T_ - 1)
  df1 <- T_ - 1
  
  stat <- (Q.2e.Beta.Computation(E.P, 
                                 Estimation.Parameters.Gamma(E.P)) -
             Q_2e) *
    df2 / Q_2e / df1
  
  l <- list("stat" = stat, "df1" = df1, "df2" = df2)
}

Z.Computating  <- function(E.P)
{
  Reduce('rbind',E.P$B$Z)
}


Matrix.Plan.Mu <- function(E.P)
{
  Matrix.Plan.Main.Effect(E.P)[, -1]
}

Matrix.Plan.Main.Effect <- function(E.P)
 
{

  I <- length(E.P$B$matrices_) 
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
 
  nu<-unlist(lapply(E.P$B$matrices_,function(x)nrow(x)));nu
  n<-sum(nu);n
  
  matrix_plan <- matrix(nrow = n, ncol = I, data = 0)
  matrix_plan[, 1] <- 1
  
  beg <- 1
  end <- nu[1]
  for(i in 1:(I - 1))
  {
    matrix_plan[beg:end, i + 1] <- 1
    beg <- end + 1
    end <- beg + nu[i+1] -1
  }
  
  p <- numeric(I - 1)
  p <- sapply(1:(I - 1), function(i) -nu[i]/sum(nu) )
              
  
  end <- n
  for(i in 1:(I - 1))
    matrix_plan[beg:end, i + 1] <- p[i]
  matrix_plan
}

Matrix.Plan.Alpha <- function(E.P)
{
  Matrix.Plan.Main.Effect(E.P)[, 1]
}

Matrix.Plan.mu <- function(E.P)
{
  Matrix.Plan.Main.Effect(E.P)[, -1]
}

Sigma.Sq.Estimation <- function(E.P)
{
  
  m__<- do.call(sum,lapply(E.P$B$First,function(x)x$Count$Mt))
  I <- length(E.P$B$matrices_)
  T_ <- ncol(E.P$B$matrices_[[1]])
  n<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  
  Q_2e <- Q.2e.Computation(E.P)
  df1 <- (m__ - n - I * (T_ - 1))
  
  Q_2e / df1
}

Sigma1.First.Calculating <- function(E.P)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  T_<-length(m[[1]]);T_
  I<-length(m);I
  n_<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  Nj<-lapply(E.P$B$First,function(x)x$Count$Nj)
  
 
  Z<-Reduce('rbind',E.P$B$Z);Z
  H <- Matrix.Plan.Main.Effect(E.P)
  
  Omega <- diag(dim(H)[1])
  R_ <- pseudoinverse(Omega)
  teta <- Param.Estimating(H, R_, Z)
  
  m__ <- sum(sapply(1:length(m), function(i) m[[i]]));m__
  df_ <- m__ - n_ - I
  Q_2e <- t(Z - H %*% teta) %*% Omega %*% (Z - H %*% teta)
  df_<-n_-I
  Q_2e / df_
}

Param.Estimating <- function(H, R_, Z)
{
  param <- pseudoinverse( t(H) %*% R_ %*% H) %*% t(H) %*% R_ %*% Z
  
  param
}

Sigma1.Sq.Estimation <- function(E.P,eps,Group)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  T_<-length(m[[1]]);T_
  I<-length(m);I
  n<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)));n
  
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  
  H <- Matrix.Plan.Main.Effect(E.P)
  
  sigma_sq <-  Sigma.Sq.Estimation(E.P)
  sigma1_sq <- Sigma1.First.Calculating(E.P)
  
  sigmas1_sq <- sigma1_sq
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R<-Omega/as.numeric(sigma1_sq)
  R_ <- solve(R)
  
  Theta <- Param.Estimating(H, R_, Z);
  
  
  for(i in 1:100)
  {
    Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
    R<-Omega/as.numeric(sigma1_sq)
    R_ <- solve(R)
    
    teta <- Param.Estimating(H, R_, Z);teta
    Theta<-rbind(Theta,teta)
    
    sigma1_sq <- Estimating.Step(Z, H, teta, R_, n, I);sigma1_sq
    if(abs(sigmas1_sq[length(sigmas1_sq)]-sigma1_sq)<0.00001)break else
    sigmas1_sq <- c(sigmas1_sq,sigma1_sq);sigmas1_sq
  }
  
  sigma1_sq
  
}

Estimating.Step <- function(Z, H, teta, R_, n, I)
{
  sigma1_sq <- t(Z - H %*% teta) %*% R_ %*% (Z - H %*% teta) / (n - I)
  sigma1_sq
}

Q.1e_.Calculating <- function(E.P)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  T_<-length(m[[1]]);T_
  I<-length(m);I
  n_<-do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  Nj<-lapply(E.P$B$First,function(x)x$Count$Nj)
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  
  x__ <- E.P$B$Second$Means$XX
  sum(do.call('c',Nj)*(Z-x__)^2)
  
}

Estimation.Parameters.Main <- function(E.P,eps,Group)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  m__ <- sum(sapply(1:length(m), function(i) m[[i]]))
  T_ <- length(m[[1]]);T_
  I <- length(m);I
  n <- do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  Nj<-lapply(E.P$B$First,function(x)x$Count$Nj)
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)

  H <- Matrix.Plan.Main.Effect(E.P)
  
  sigma_sq <- Sigma.Sq.Estimation(E.P)
  sigma1_sq <- Sigma1.Sq.Estimation(E.P,eps,Group)
  
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R_ <- pseudoinverse(Omega/as.numeric(sigma1_sq))
  
  teta <- Param.Estimating(H, R_, Z)
  
  teta
}

Estimation.Alpha <- function(E.P,Theta)
{
  teta <- Theta
  alpha <- teta[-1]
  
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  
  mI <- lapply(m, sum)
  T_ <- length(m[[1]]);T_
  I <- length(m);I
 
  alpha <- c(alpha, -sum(sapply(1:(I - 1), 
             function(i) mI[[i]] * alpha[i])) / mI[[I]])
  alpha
}

Estimation.Mu <- function(E.P,Theta)
{
  teta <- Theta
  alpha <- teta[1]
}

Alpha.Statistics <- function(E.P,eps,Theta,Group)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  m__ <- sum(sapply(1:length(m), function(i) m[[i]]))
  T_ <- length(m[[1]]);T_
  I <- length(m);I
  n <- do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  Nj<-lapply(E.P$B$First,function(x)x$Count$Nj)
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  CV
  
 
 
  H <- Matrix.Plan.Main.Effect(E.P)
  H_ <- Matrix.Plan.Alpha(E.P)
  
  
  sigma_sq <- Sigma.Sq.Estimation(E.P)
  sigma1_sq <- Sigma1.Sq.Estimation(E.P,eps,Group)
  Sigma<-c(sigma_sq,sigma1_sq)
  
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R_ <- pseudoinverse(Omega/as.numeric(sigma1_sq))
  R_[seq(5),seq(5)]
  teta <- Theta
  mu <- Param.Estimating(H_, R_, Z)
  
  df1 <- I - 1
  df2 <- n - I
  
  stat <- (t(Z - H_ %*% mu) %*% R_ %*% (Z - H_ %*% mu) - 
             t(Z - H %*% teta) %*% R_ %*% (Z - H %*% teta)) * df2 / 
    (t(Z - H %*% teta) %*% R_ %*% (Z - H %*% teta)) / df1
  
  list("stat" = stat, "df1" = df1, "df2" = df2,Sigma=Sigma,
       R_=R_,
       HT=H %*% teta,
       HM=H_ %*% mu)
}

P.Value.Computating <- function(stat, n = n, nm = nm)
{
  1 - pf(q = stat$stat/(n*sqrt(log(n+1, base = 2))), df1 = stat$df1, df2 = stat$df2)
}

Main.Statistics <- function(E.P,Theta,Sigma,Group)
{
  m<-lapply(E.P$B$First,function(x)x$Count$Mt);m
  m__ <- sum(sapply(1:length(m), function(i) m[[i]]))
  T_ <- length(m[[1]]);T_
  I <- length(m);I
  n <- do.call(sum,lapply(E.P$B$matrices,function(x)nrow(x)))
  Nj<-lapply(E.P$B$First,function(x)x$Count$Nj)
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  
  H <- Matrix.Plan.Main.Effect(E.P)
  H_ <- Matrix.Plan.Alpha(E.P)
  H__ <- Matrix.Plan.Mu(E.P)

  sigma_sq <- Sigma[1]
  sigma1_sq <- Sigma[2]
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R_ <- pseudoinverse(Omega/as.numeric(sigma1_sq))
  
  teta<-Theta
  alpha <- Param.Estimating(H__, R_, Z)
  
  df1 <- 1
  df2 <- n - I
  
  stat <- (t(Z - H__ %*% alpha) %*% R_ %*% (Z - H__ %*% alpha) - t(Z - H %*% teta) %*% R_ %*%
             (Z - H %*% teta)) * df2 / (t(Z - H %*% teta) %*% R_ %*% (Z - H %*% teta)) / df1

  list("stat" = stat, "df1" = df1, "df2" = df2)
}

Rest.Shapiro<-function(data,ep,nm)
{
  y_<-NULL
  for(i in seq(length(ep$B$Z)))
    y_<-c(y_,ep$B$Z[[i]]-ep$B$Second$Means$XI[i])
  
  XIT<-ep$B$Second$Means$XIT
  
  n_<-ceiling(length(y_)/nm)
  
  LL<-lapply(seq(30),function(c){
    list(y__<-y_[sample(length(y_))[seq(n_)]],p=shapiro.test(y__)$p.value) 
  })
  
  ll0<-LL[which.max(unlist(lapply(LL,function(x)x[[2]])))]
   
  
  hist(ll0[[1]][[1]],freq=FALSE,main="histogram of rest 1",xlab="",
       sub = paste("p" , format(ll0[[1]][[2]],3,2),sep="="))
  curve(dnorm(x,mean(y_),sd(y_)),add=TRUE,col=2)
  
    Y_<-ep$Y-ep$plan %*% ep$teta
    
    y_<-Y_
    LL<-lapply(seq(30),function(c){
      list(y__<-y_[sample(length(y_))[seq(n_)]],p=shapiro.test(y__)$p.value) 
    })
    
    ll0<-LL[which.max(unlist(lapply(LL,function(x)x[[2]])))]
    
    
    hist(ll0[[1]][[1]],freq=FALSE,main="histogram of rest 2",xlab="",
         sub = paste("p" , format(ll0[[1]][[2]],3,2),sep="="))
    curve(dnorm(x,mean(y_),sd(y_)),add=TRUE,col=2)
 
 }

AU.2<-function(data,eps,Group, n = 1, nm = 2)
{
  matrices <- matrices.list(as.matrix(data), data)
  length(matrices)
  #оценивание параметров
  ep <- Estimation.Parameters(matrices,Group)
  ep$B$Second$Means$XIT
  names(ep)
  B<-ep$B
  
#эффект времени
  beta <- Beta.Computing(ep);beta
  #статистика критерия со степенями свободы
  F_b <- F.Beta(ep)
  #p-value
  p_value_b <- P.Value.Computating(F_b, n = n, nm = nm)
  p_value_b
  #эффект взаимодействия группы и времени
  gamma <- Gamma.Computing(ep)
  #статистика критерия со степенями свободы
  F_g <- F.Gamma(ep)
  #p-value
  p_value_g <- P.Value.Computating(F_g, n = n, nm = nm)
  p_value_g
  #эффект группы
  Theta<-Estimation.Parameters.Main(ep,eps,Group)
  alpha <- Estimation.Alpha(ep,Theta)
  #статистика критерия со степенями свободы
  F_alpha <- Alpha.Statistics(ep,eps,Theta,Group)
  Sigma<-F_alpha$Sigma
  #p-value
  p_value_alpha <- P.Value.Computating(F_alpha, n = n, nm = nm)
  
  #генеральное среднее
  mu <- Estimation.Mu(ep,Theta)
  #статистика критерия со степенями свободы
  F_mu <- Main.Statistics(ep,Theta,Sigma,Group)
  #p-value
  p_value_mu <- P.Value.Computating(F_mu, n = n, nm = nm)
  
  MU<-list(mu=mu,Alpha=alpha,Beta=beta,gamma=gamma)
  P.values<-c(
    p_value_mu=p_value_mu,
    p_value_alpha=p_value_alpha,
    p_value_b=p_value_b,
    p_value_g=p_value_g)
  list(ep=ep,MU=MU,P.values=P.values,data=data)
}

matrix.plan <- function(L)
{
  I <- length(L$m)
  T_ <- length(L$m[[1]])
  Y_length <- sum(sapply(1:I, function(i) length(L$Y[[i]])))
  plan <- matrix(data = 0, ncol = (T_ - 1) * I)
  
  q <- numeric(T_ - 1)
  
  s <- sum(sapply(1:I, function(i) L$m[[i]][T_]))
  
  for(t in 1:(T_ - 1))
  {
    q[t] <- - sum(sapply(1:I, function(i) L$m[[i]][t])) / s
  }
  
  f <- function(i)
  {
    Y <- L$Y[[i]]
    m <- L$m[[i]]
    plan_i <- matrix(data = 0, nrow = length(Y), ncol = (T_ - 1) * I)
    
    q_i <- -m / m[length(m)]
    q_i <- q_i[-T_]
    
    beg <- 0
    end <- 0
    
    for(j in 1:(T_ - 1))
    {
      beg <- end + 1
      end <- end + m[j]
      ind <- 1 + (j - 1) * I
      if(I > 1)
      {
        plan_i[beg:end, c(ind, ind + i)] <-  1
      }
      plan_i[beg:end, ind] <- 1
    }
    
    for(k in 1:(T_ - 1))
    {
      if(I > 1) plan_i[(end + 1):length(Y), 1 + i + (k - 1) * I] <- q_i[k]
      plan_i[(end + 1):length(Y), 1 + (k - 1) * I] <- q[k]
    }
    
    plan_i
  }
  
  if(I > 1)
  {
    for(i in 1:(I - 1))
    {
      plan <- rbind(plan, f(i))
    }
    
    Y <- L$Y[[I]]
    m <- L$m[[I]]
    plan_I <- matrix(data = 0, nrow = length(Y), ncol = (T_ - 1) * I)
    p <- matrix(0, nrow = I - 1, ncol = T_ - 1)
    r <- matrix(0, nrow = I - 1, ncol = T_ - 1)
    
    for(i in 1:(I - 1))
      for(t in 1:(T_ - 1))
      {
        p[i,t] <- -L$m[[i]][t] / m[t]
        r[i,t] <- L$m[[i]][t] / m[T_]
      }
    
    beg <- 0
    end <- 0
    
    for(j in 1:(T_ - 1))
    {
      beg <- end + 1
      end <- end + m[j]
      ind <- 1 + (j - 1) * I
      plan_I[beg:end,ind] <-  1
      for(k in 1:(I - 1))
      {
        plan_I[beg:end, 2 + (j - 1) * I + (k - 1)] <-  p[k,j]
      }
    }
    
    for(k in 1:(T_ - 1))
    {
      plan_I[(end + 1):length(Y), 1 + (k - 1) * I] <- q[k]
    }
    
    for(j in 1:(T_ - 1))
      for(k in 1:(I - 1)) plan_I[(end + 1):length(Y), 2 + (j - 1) * I + (k - 1)] <- r[k,j]
    
    plan <- rbind(plan, plan_I)
    
  }
  
  else
  {
    plan <- rbind(plan, f(I))
  }
  
  plan <- plan[-1,]
  plan
}

matrices.list <- function(mx, df)
{
  fact <- as.factor(df$group)
  I <- length(levels(fact))
  if(I > 1)
    l <- lapply(1:I, function(i) subset(mx, df$group == i, -1)) else list(mx[ , -1])
}

H.Gamma.Computing <- function(plan, I, T_)
{
  H_gamma <- 0
  
  for(t in 1:(T_ - 1))
    H_gamma <- cbind(H_gamma, plan[, (2 + (t - 1) * I):(I * t)])
  
  H_gamma <- subset(H_gamma,select=-H_gamma)
  
  H_gamma
}

H.beta.computing <- function(plan, T_, I)
{
  H_beta <- sapply(1:(T_ - 1), function(t) plan[, 1 + (t - 1) * I])
}

f0.f<-function(iota,i,j,t,Second){
  A<-Second$AB$A;MM<-Second$MM
  ifelse(iota==i,1,0)*(Second$AB$B[iota,t]/colSums(MM)[t]-A[iota,i]/rowSums(MM)[i])
}

F1_<-function(i,First,Second){
  MM<-Second$MM
  
  BB<-Second$AB$B %*% diag(1/colSums(MM))
  C1<-First[[i]]$Lambda$Lambda_nu %*% First[[i]]$J %*% t(BB)
  
  C_<-matrix(rep(Second$AB$A[,i]/ rowSums(MM) [i],nrow(C1)),
             byrow=TRUE,ncol=nrow(BB))
  C1-C_ 
}

F2_<-function(i,Second)
{
  MM<-Second$MM
  BB<-Second$AB$B %*% diag(1/colSums(MM));BB
  AA<-Second$AB$A %*% diag(1/rowSums(MM));AA
  BB-AA[,i]
}

DDD__<-function(i,k,B,Group){
  Nu<-unlist(lapply(B$First,function(x)length(x$Counts$Nj)));Nu
  SS<-zeros(Nu[i],Nu[k]);SS
  if(Group==TRUE) SS<-SS+F0_(B$Second)[i,k] 
  SS}

Covariance.matrix.Omega.1<-function(B,Group)
{
  I<-length(lapply(B$First,function(x)x$Counts$Nj));I
  Reduce('rbind',lapply(seq(I),function(i1)
    Reduce('cbind',lapply(seq(I),function(i2)
      WT(i1,i2,B,Group)
      ))))
    
}

WT<-function(i1,i2,B,Group){
  if(i1==i2){
    res<-delta.w0(B)[[i1]]+delta.f(c(1,1),B,c(i1,i2),Group) } else 
    res<-delta.f(t=c(1,1),B,c(i1,i2),Group);
  res
}

F0_<-function(Second)
{
  MM<-Second$MM;MM
  B_<-Second$AB$B;BB<-B_ %*% diag(1/(colSums(MM)));BB
  
  A<-Second$AB$A;AA<-A %*% diag(1/(rowSums(MM)))
  B_%*%t(BB)-A%*%t(AA)
}

Plot<-function(res,Name,data,Inv)
{
  if(Inv==TRUE)res$MU$gamma<-res$MU$gamma[order(seq(nrow(res$MU$gamma)),decreasing = TRUE),]
  plot(seq(ncol(res$MU$gamma)),c(min(res$MU$gamma),rep(0,ncol(res$MU$gamma)-2),max(res$MU$gamma)),
       type="n",ylab=paste("gamma",Name,sep="."),lty=2,xlab="Time")
  for(i in 1:nrow(res$MU$gamma))lines(res$MU$gamma[i,],type="b",
                                      pch=20,col=i)
  Names<-names(table(data$group))
  legend("right",Names,col=seq(3),pch=20,cex=0.6)
  title(sub=paste("p",round(res$P.values[4],digits=4),sep="="))
}

Plot.A<-function(res,Name)
{
  plot(res$MU$Alpha,type="b",ylab=paste("alpha",Name,sep="."),xaxt = "n")
  axis(1, 1:2, c(1,3))
  title(sub=paste("p",round(res$P.values[2],digits=4),sep="="))
}

Plot.B<-function(res,Name,Inv)
{
  y<-res$MU$Beta
  z<-seq(length(y))
  if(Inv==TRUE) z<-z[order(z,decreasing=TRUE)]
  plot(y[z],type="b",ylab=paste("beta",Name,sep="."))
  
  title(sub=paste("p",round(res$P.values[3],digits=4),sep="="))
}

file_<-function(data.,group,Name)
  {
  c<-sort(na.omit(unique(group)));c
  group_<-sapply(group,function(x)ifelse(!is.na(x),which(c==x),NA))
  
  data_<-data.frame(group=group_,
  data.[,which(sapply(colnames(data.),
                     function(names)
                    substr(names,1,nchar(names)-2))==Name)])
  data_
}

Plot.Total<-function(data,group,Name)
{
  dat.0<-file_(data,group,Name)
  aqq<-sapply(names(table(dat.0$group)),function(nn)
    colMeans(dat.0[dat.0$group==nn,-1],na.rm=TRUE))
  
  c<-c(rep(min(aqq,na.rm=TRUE),nrow(aqq)-1),max(aqq,na.rm=TRUE))
  
  plot(seq(nrow(aqq)),c,type="n",
       ylab=Name,xlab="time",axes=FALSE)
  for(i in seq(ncol(aqq)) )
      lines(aqq[,i],type="b",col=i)
 axis(1,seq(nrow(aqq)),seq(nrow(aqq)))
 axis(2)
  legend('right',LETTERS[seq(ncol(aqq))],pch=1,col=seq(ncol(aqq)),cex=0.7)
}

Plot.Total_<-function(dat.0)
{
  aqq<-sapply(names(table(dat.0$group)),function(nn)
    colMeans(dat.0[dat.0$group==nn,-1],na.rm=TRUE))
  colnames(aqq)<-paste("m",seq(3))
  sqq<-sapply(names(table(dat.0$group)),function(nn)
    
    apply(dat.0[dat.0$group==nn,-1],2,function(x)sd(x,na.rm=TRUE))
    )
  colnames(sqq)<-paste("sd",seq(3))
  
  nqq<-sapply(names(table(dat.0$group)),function(nn)
    
    apply(dat.0[dat.0$group==nn,-1],2,function(x)length(na.omit(x)))
  )
  colnames(nqq)<-paste("n",seq(3))
  c<-c(rep(min(aqq),nrow(aqq)-1),max(aqq))
  eqq<-sqq/sqrt(nqq)
  colnames(eqq)<-paste("err",seq(3))
  
  plot(seq(nrow(aqq)),c,type="n",
       ylab=Name,xlab="time",axes=FALSE)
  for(i in seq(ncol(aqq)) )
    lines(aqq[,i],type="b",col=i)
  axis(1,seq(nrow(aqq)),seq(nrow(aqq)))
  axis(2)
  legend('right',as.character(seq(ncol(aqq))),pch=1,col=seq(ncol(aqq)),cex=0.7)
  cbind(aqq,eqq,nqq)
}

Plot.Total_2<-function(dat.0,Name,Name.gr)
{
  aqq<-sapply(names(table(dat.0$group)),function(nn)
    colMeans(dat.0[dat.0$group==nn,-1],na.rm=TRUE))
  colnames(aqq)<-seq(ncol(aqq))-1
  c<-c(rep(min(aqq),nrow(aqq)-1),max(aqq));c
  
  plot(seq(nrow(aqq)),c,type="n",
       ylab=Name,xlab="time",axes=FALSE)
  for(i in seq(ncol(aqq)) )
    lines(aqq[,i],type="b",col=i)
  axis(1,seq(nrow(aqq)),seq(nrow(aqq)))
  axis(2)
  leg<-paste(Name.gr,names(table(dat.0$group)))
  legend('right',leg,pch=1,col=seq(ncol(aqq)),cex=0.7)
}


Plot.Total_3<-function(data,group,Name,Name.gr,tit.m,tit.s)
{
  dat.0<-file_(data,group,Name)
  
  aqq<-sapply(names(table(dat.0$group)),function(nn)
    colMeans(dat.0[dat.0$group==nn,-1],na.rm=TRUE))
  
  c<-c(rep(min(aqq),nrow(aqq)-1),max(aqq))
  plot(seq(nrow(aqq)),c,type="n",
       ylab=Name,xlab="time",axes=FALSE)
  for(i in seq(ncol(aqq)) )
    lines(aqq[,i],type="b",col=i)
  axis(1,seq(nrow(aqq)),colnames(data))
  axis(2)
  leg<-paste(Name.gr,names(table(dat.0$group)))
  legend('bottomright',leg,pch=1,col=seq(ncol(aqq)),cex=0.5)
  title(sub=tit.s,main=tit.m)
}

ARM<-function(Base,group,Name)

{
  data<-file_(Base,group,Name)
  data<-data[,c("group",as.character(na.omit(sapply(colnames(data),
                                                    function(x)if(substr(x,nchar(x),
                                                                         nchar(x)) %in% as.character((seq(10)-1))) x else NA))))]
  data<-data[complete.cases(data[,seq(2)]),]
  if(ncol(data)>2) {
    pr<-apply(data[,-1],2, 
              function(x)length(x[complete.cases(x)])/nrow(data) )
    if( min(pr)>0.2 & length(table(data[,2]))>10)
    {
      res..<-AU.2(data,eps=0.000001,Group=TRUE)
      res..$MU
      res<-round(res..$P.values,4)} else res<-rep(NA,4)} else 
        res<-rep(NA,4)
  res
}

ARM.Covar<-function(Base,group,Name,Name.Covar)
{
  data<-file_(Base,group,Name)
  View(data)
  data<-cbind(data[,c("group",
                as.character(na.omit(sapply(colnames(data),
            function(x)if(substr(x,nchar(x),
          nchar(x)) %in% as.character((seq(10)-1))) x else NA)))
          )],cov=Base[,Name.Covar])
  
  data<-data[complete.cases(data[,seq(2)]),]
  
  data<-data.frame(group=data$group,
                   apply(subset(data,select=-c(group,cov)),2,
       function(x){x-predict(lm(x~data$cov),
                             subset(data,select=-c(group,cov)))+
                     mean(x,na.rm=TRUE)}))
  
  
  if(ncol(data)>2) {
    pr<-apply(data[,-1],2, 
              function(x)length(x[complete.cases(x)])/nrow(data) )
    if( min(pr)>0.2 & length(table(data[,2]))>10)
    {
      data<-data[complete.cases(data[,seq(2)]),]
      res..<-AU.2(data,eps=0.000001,Group=TRUE)
      res<-round(res..$P.values,4)} else res<-rep(NA,4)}else 
        res<-rep(NA,4)
  res
}

ARM.Covar.Full<-function(Base,group,Name,Name.Covar)
{
  data<-file_(Base,group,Name)
  data.0<-data.frame(data[,c("group",
                      as.character(na.omit(sapply(colnames(data),
                                                  function(x)if(substr(x,nchar(x),
                                                                       nchar(x)) %in% as.character((seq(10)-1))) x else NA)))
  )])
  if(length(Name.Covar)==1) data.1<-data.frame(cov=Base[,Name.Covar])else
  data.1<-Base[,Name.Covar]

  data.0<-subset(data.0,select=-group)
  
  data<-data.frame(group=data$group,
       apply(data.0,2,
      function(x){
        dat_<-data.frame(x=x,data.1);
        x-predict(lm(x~.,dat_),dat_)+mean(x,na.rm=TRUE)}))
                                  
                                       
  
  data<-data[complete.cases(data[,seq(2)]),]
  if(ncol(data)>2) {
    pr<-apply(data[,-1],2, 
              function(x)length(x[complete.cases(x)])/nrow(data) )
    if( min(pr)>0.2 & length(table(data[,2]))>10)
    {
      data<-data[complete.cases(data[,seq(2)]),]
      res..<-AU.2(data,eps=0.000001,Group=TRUE)
      res<-res..} else res<-rep(NA,4)}else 
        res<-rep(NA,4)
  res
}

Mick.Computation.0<-function(t,B,ii,Group)
{
  
  Nj<-lapply(B$First,function(x)x$Count$Nj);Nj             
  
  nu<-unlist(lapply(Nj,function(x)length(x)));nu
  C<-B$First[[ii[1]]]$C;D<-B$First[[ii[1]]]$D;C;D
  Mt<-lapply(B$First,function(x)x$Counts$Mt)
  MM<-Reduce('rbind',Mt);MM;colSums(MM)
  JJ<-lapply(B$First,function(x)x$J);JJ
  
  
  if(ii[1]==ii[2])
    W.A<-sapply(seq(nu[ii[1]]),function(l)
      sapply(seq(nu[ii[2]]),function(j) 
        C[j,l]/Nj[[ii[1]]][l]- D[j,t]/MM[ii[1],t] )  ) else
          W.A<-zeros(nu[[ii[1]]],nu[[ii[2]]])
  W.A
  if(Group==TRUE & ii[1]==ii[2])
    W.B<- B$Second$AB$A[ii[1],ii[2]]/rowSums(MM)[ii[2]]-
    B$Second$AB$B[ii[1],t]/colSums(MM)[t] else W.B<-0
  
  W<-W.A+W.B;W
  W
}

Mick.T<-function(tt,B,ii,Group)
{
  JJ<-lapply(B$First,function(x)x$J);JJ
  W<-(Mick.Computation.0(t=tt[1],B,ii=ii,Group)+
  t(Mick.Computation.0(t=tt[2],B,ii=ii[c(2,1)],Group)))*(-1);
  Cond1<-JJ[[ii[1]]][,tt[1]]==1;Cond2<-JJ[[ii[2]]][,tt[2]]==1
  W_<-t(subset(t(subset(W,Cond1)),Cond2));dim(W_);round(W_,2)
  W_
}

CovEE<-function(i,j,t,k,l,tau,B,Group)
{
  JJ<-lapply(B$First,function(x)x$J);JJ
  Ja<-JJ[[i]][,t]*JJ[[k]][,tau]
  
  Delta<-delta.w0(B);Delta
  
  d<-ifelse(i==k,Delta[[i]][j,l],0)+F0_(B$Second)[i,k]+ifelse(j==l & i==k & t==tau,1,0);d

  d<-d-(Mick.Local.0(i,j,k,l,tau,B,Group)+Mick.Local.0(k,l,i,j,t,B,Group))#*Ja[j]*Ja[l]
  d

}

CovET<-function(tt,ii,B,Group){
  Nj<-lapply(B$First,function(x)x$Count$Nj);Nj 
  nu<-unlist(lapply(Nj,function(x)length(x)));nu
  JJ<-lapply(B$First,function(x)x$J);JJ
  i<-ii[1];k<-ii[2];
  if(i==k)
    {W<-sapply(seq(nu[i]),function(j)
      sapply(seq(nu[k]),function(l) 
       CovEE(i,j,tt[1],k,l,tt[2],B,Group)  ))} else W<-zeros(nu[1],nu[2])
  W
  
}

CovETC<-function(JJ,W,tt)
{
  Cond1<-JJ[[ii[1]]][,tt[1]]==1;Cond2<-JJ[[ii[2]]][,tt[2]]==1
  W_<-t(subset(t(subset(W,Cond1)),Cond2));dim(W_);round(W_,2)
  W_
}


Mick.Local.0<-function(i,j,k,l,tau,B,Group)
{
  
  Nj<-lapply(B$First,function(x)x$Count$Nj);Nj             
  JJ<-lapply(B$First,function(x)x$J);JJ
  nu<-unlist(lapply(Nj,function(x)length(x)));nu
  C<-B$First[[i]]$C;
  D<-B$First[[i]]$D;D
 
  Mt<-lapply(B$First,function(x)x$Counts$Mt)
  MM<-Reduce('rbind',Mt);MM;colSums(MM)
  
  #*JJ[[i]][j,t]
  if(i==k ) w.a<-(C[j,l]/Nj[[i]][l]- D[j,t]/MM[i,t]) else w.a<-0
  w.a
  
  if(Group==TRUE & i==k)
    w.b<- B$Second$AB$A[i,k]/rowSums(MM)[k]-
    B$Second$AB$B[i,t]/colSums(MM)[t] else w.b<-0
  
  w<-w.a+w.b

  w
  
}

Cov.F_<-function(ii,B,Group)
{
  MM<-Reduce('rbind',lapply(B$First,function(x)x$Counts$Mt));MM
  T_<-ncol(MM)
  C<-do.call('rbind',sapply(seq(T_),function(t)
  list(do.call('cbind',
          sapply(seq(T_),function(tau)
      {
            JJ<-lapply(B$First,function(x)x$J)
            W<-CovET(c(t,tau),ii,B,Group)
            list(  CovETC( JJ,W, c(t,tau))   )}
      )
      )
  ))
  )
  
  
  
C
}

Wilks <- function(p, v_H, v_E, E.P) #p-- число признаков; v_H-- степени свободы H; v_E-- степени свободы E
{
  Lambda_r <- E.P$Lambda_r
  lambda_ <- (p*v_H-2)/4
  if (p*p + v_H*v_H > 5){
    d <- sqrt((p*p*v_H*v_H-4)/(p*p + v_H*v_H - 5))
  }
  else{
    d <- 1
  }
  f <- v_E - (p - v_H + 1)/2
  (1-Lamda_r^(1/d))/Lamda_r^(1/d)*(f*d-2*lambda_)/(p*v_H)
}

R_0 <- function(Z, H, Teta) #H=матрица плана, для R_0 и учеченная для R_1
{
  eps <- 0.001 #для критерия останова
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  
  H <- Matrix.Plan.Main.Effect(E.P)
  
  sigma_sq <-  Sigma.Sq.Estimation(E.P)
  sigma1_sq <- Sigma1.First.Calculating(E.P)
  
  sigmas1_sq <- sigma1_sq
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R<-Omega/as.numeric(sigma1_sq)
  R_ <- solve(R)
  
  Theta <- Param.Estimating(H, R_, Z);
  while(abs(sigmas1_sq[length(sigmas1_sq)]-sigma1_sq) > eps){ #расстояние Фробениуса
    Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
    R<-Omega/as.numeric(sigma1_sq)
    R_ <- solve(R)
    
    teta <- Param.Estimating(H, R_, Z);teta
    Theta<-rbind(Theta,teta)
    
    sigma1_sq <- Estimating.Step(Z, H, teta, R_, n, I);sigma1_sq
    sigmas1_sq <- c(sigmas1_sq,sigma1_sq);sigmas1_sq
  }
  t(Y-H%*%Teta)%*%sigmas1_sq%*%(Y-H%*%Teta)
}

R_1 <- function(Z, H, Teta) #H=матрица плана, для R_0 и учеченная для R_1
{
  eps <- 0.001 #для критерия останова
  Z<-Reduce('rbind',E.P$B$Z);Z
  CV<-Covariance.matrix.Omega.1(E.P$B,Group)
  
  H <- Matrix.Plan.Main.Effect(E.P)
  
  sigma_sq <-  Sigma.Sq.Estimation(E.P)
  sigma1_sq <- Sigma1.First.Calculating(E.P)
  
  sigmas1_sq <- sigma1_sq
  Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
  R<-Omega/as.numeric(sigma1_sq)
  R_ <- solve(R)
  
  Theta <- Param.Estimating(H, R_, Z);
  while(abs(sigmas1_sq[length(sigmas1_sq)]-sigma1_sq) > eps){ #расстояние Фробениуса
    Omega <-  as.numeric(sigma1_sq)*eye(n)+as.numeric(sigma_sq)*CV
    R<-Omega/as.numeric(sigma1_sq)
    R_ <- solve(R)
    
    teta <- Param.Estimating(H, R_, Z);teta
    Theta<-rbind(Theta,teta)
    
    sigma1_sq <- Estimating.Step(Z, H, teta, R_, n, I);sigma1_sq
    sigmas1_sq <- c(sigmas1_sq,sigma1_sq);sigmas1_sq
  }
  t(Y-H%*%pseudo_inverse(t(H)%*%sigmas1_sq*H))%*%sigmas1_sq%*%(Y-H%*%pseudo_inverse(t(H)%*%sigmas1_sq*H))
}