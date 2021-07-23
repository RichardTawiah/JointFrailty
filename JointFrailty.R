options(max.print=1000000)
library(dplyr)
#setwd("C:/Users/r.tawiah/OneDrive/comorbidity_data_Angus/Biostatistics_paper/Review/Clean_code")
data<-read.table("data_GitHub.txt",header=TRUE)
dat=as.matrix(data)

joint.frailty<-function(dat,patient,theta01,theta02,rho0,itmax){
  
  p<-ncol(dat)-5
  n<-nrow(dat)
  m<-length(unique(patient))
  
  Z<-matrix(0,ncol=m,nrow=n)
  for(j in 1:m){
    Z[,j]<-ifelse(patient==j,1,0)
  }

  XZ<-data.frame(cbind(dat,Z))  # all recurrences
  XZ.R<-as.matrix(XZ)
  XZ.R.r<-XZ.R[sort.list(XZ.R[,3]),]  # reorder recurrent gap times 
  time1<-as.vector(XZ.R.r[,3]) # gap time and its censoring indicator
  indi1<-as.vector(XZ.R.r[,4])
  
  nrpp<-table(dat[,2])
  XZ.D<-data.frame(cbind(dat,Z))   
  XZ.D<-XZ.D %>% 
    group_by(id) %>% 
    filter( ((delta2)==0 & row_number()==n())| 
              ((delta2)==1 & row_number()==n())) 

  XZ.D<-as.matrix(XZ.D)
  XZ.D.r<-XZ.D[sort.list(XZ.D[,3]),] 
  time2<-as.vector(XZ.D.r[,3])            
  indi2<-as.vector(XZ.D.r[,5])
  
  Z.R<-as.matrix(XZ.R.r[,-(1:(5+p))])
  Z.D<-as.matrix(XZ.D.r[,-(1:(5+p))])
  X.R<-as.matrix(XZ.R.r[,6:(5+p)])
  X.D<-as.matrix(XZ.D.r[,6:(5+p)])

  M1<-diag(0,n)
  for( j in 1:n)
    for(i in 1:j) M1[j,i]<-1
  
  M2<-diag(0,m)
  for( j in 1:m)
    for(i in 1:j) M2[j,i]<-1

#initial values
  beta0<-as.vector(rep(0,p))
  gamma0<-as.vector(rep(0,p))
  u<-as.vector(rep(0,m))
  v<-as.vector(rep(0,m))

  par0.lat<-as.vector(c(beta0,gamma0,u,v))
  
  eta.R<-as.vector(X.R%*%beta0+Z.R%*%u)
  eta.D<-as.vector(X.D%*%gamma0+Z.D%*%v)
  
  flag.var<-0
  eps.reg<-0.0001
  eps.var<-0.005
# change eps.var from 0.005 to 0.01 for decreasing the number of iterations
  
# Lai 2008 approach
K1<-rbind(cbind(diag(m),diag(0,m)),cbind(diag(0,m),diag(0,m)))
K2<-rbind(cbind(diag(0,m),diag(m)),cbind(diag(m),diag(0,m)))
K3<-rbind(cbind(diag(0,m),diag(0,m)),cbind(diag(0,m),diag(m)))
  
  for(iteration in 1:itmax){

    cat("iteration=",iteration,"\n")
    
    flag.reg<-0
    
    UG<-diag(0,(p+p+m+m))
    UG[(p+p+1):(p+p+m),(p+p+1):(p+p+m)]<-(theta02)*diag(m)
    UG[(p+p+m+1):(p+p+m+m),(p+p+1):(p+p+m)]<--(rho0*(sqrt(theta01*theta02)))*diag(m)
    UG[(p+p+1):(p+p+m),(p+p+m+1):(p+p+m+m)]<--(rho0*(sqrt(theta01*theta02)))*diag(m)
    UG[(p+p+m+1):(p+p+m+m),(p+p+m+1):(p+p+m+m)]<-(theta01)*diag(m)
    UG<-(1/(theta01*theta02*(1-rho0^2)))*UG
    
    for(it.m.step in 1:itmax){
      
      eta.R<-as.vector(X.R%*%beta0+Z.R%*%u)
      eta.D<-as.vector(X.D%*%gamma0+Z.D%*%v)
  
    survbase.R<-exp(-M1%*%(indi1/(t(M1)%*%(exp(eta.R)))))
    survbase.D<-exp(-M2%*%(indi2/(t(M2)%*%(exp(eta.D)))))

    bb<- as.vector(rep(0,m))
    for( j in 1:m) bb[j]<-sum(exp(eta.D[which(time2>=time2[j])]))
    survDnew<- exp(-M2%*%(indi2/(bb)))

    dd<- as.vector(rep(0,n))
    for( j in 1:n) dd[j]<-sum(exp(eta.R[which(time1>=time1[j])]))
    survRnew<- exp(-M1%*%(indi1/(dd)))
      
    w.R<-as.vector(exp(eta.R))

    A.R<-indi1/(dd)
    B.R<-cumsum(A.R)
    amw<-(t(M1*w.R)*A.R)

    A.Rnew<-indi1/(dd)
    B.R<-cumsum(A.Rnew)
    amw<-(t(M1*w.R)*A.Rnew)      
      
    f1.eta.R<-as.vector(indi1-t(amw)%*%rep(1,n))
    f2.eta.R<-diag(w.R*B.R)-t(amw)%*%amw  
    w.D<-as.vector(exp(eta.D))
    A.D<-indi2/(bb)      

    B.D<-cumsum(A.D)
    amw<-(t(M2*w.D)*A.D)

    A.Dnew<-indi2/(bb)
    B.D<- cumsum(A.Dnew)
    amw<-(t(M2*w.D)*A.Dnew)

    f1.eta.D<-as.vector(indi2-t(amw)%*%rep(1,m))
    f2.eta.D<-diag(w.D*B.D)-t(amw)%*%amw   
      
      dl.dbeta<-t(X.R)%*%f1.eta.R
      dl.dgamma<-t(X.D)%*%f1.eta.D
      dl.du<-t(Z.R)%*%f1.eta.R-(u*theta02-v*rho0*sqrt(theta01*theta02))/
        (theta01*theta02*(1-rho0^(2)))
      dl.dv<-t(Z.D)%*%f1.eta.D-(v*theta01-u*rho0*sqrt(theta01*theta02))/
        (theta01*theta02*(1-rho0^(2)))
      
      O.np<-matrix(rep(0,n*p), nrow=n,byrow=TRUE)
      O.mp<-matrix(rep(0,m*p), nrow=m,byrow=TRUE)
      O.nm<-matrix(rep(0,n*m), nrow=n,byrow=TRUE)
      O.mm<-matrix(rep(0,m*m), nrow=m,byrow=TRUE)
      
      XX.R<-cbind(X.R,O.np,Z.R,O.nm)
      XX.D<-cbind(O.mp,X.D,O.mm,Z.D)
      XX<-rbind(XX.R,XX.D)
      O.f1<-matrix(rep(0,n*m), nrow=n,byrow=TRUE)
      O.f2<-matrix(rep(0,m*n), nrow=m,byrow=TRUE)
      one<-cbind(f2.eta.R,O.f1)
      two<-cbind(O.f2,f2.eta.D)

      f2.eta.R.D<-rbind(one,two)
      HH<-t(XX)%*%f2.eta.R.D%*%XX+UG
      H1<-solve(HH)
      
      Svec.lat<-as.vector(c(dl.dbeta,dl.dgamma,dl.du,dl.dv))
      
      par.lat<-par0.lat+H1%*%Svec.lat     
    
      if(max(abs(c((par.lat-par0.lat))))<eps.reg){flag.reg<-1;break}
      
      par0.lat<-par.lat
      beta0<-par.lat[1:p]
      gamma0<-par.lat[(p+1):(p+p)]
      u<-par.lat[(p+p+1):(p+p+m)]
      v<-par.lat[(p+p+m+1):(p+p+m+m)]
   
      q<-par.lat[(p+p+1):(p+p+m+m)]
      
      cat("beta0=",beta0,"gamma0=",gamma0,'\n')
      
    } 
    if(flag.reg==0)stop("not reach convergence")
    flag.reg<-0


  tau<-(q%*%t(q)+H1[(p+p+1):(p+p+m+m),(p+p+1):(p+p+m+m)])  
  
  B1<-sum(diag(K1%*%tau))
  B2<-sum(diag(K2%*%tau))/2
  B3<-sum(diag(K3%*%tau))
  
  theta1<-B1/m; 
  theta2<-B3/m; 
  rho<-B2/sqrt(B1*B3)

 cat("theta1=",theta1,"theta2=",theta2,"rho=",rho,'\n')

if(max(abs(c((theta1-theta01),(theta2-theta02),(rho0-rho))))<eps.var){flag<-1;break}
 
      theta01<-theta1
      theta02<-theta2
      rho0<-rho
  }
  A<-theta01*diag(m); B<-rho0*sqrt(theta01*theta02)*diag(m);
  C<-rho0*sqrt(theta01*theta02)*diag(m); D<-theta02*diag(m)
  A1<-cbind(A,B); A2<-cbind(C,D)

  Sigma<-matrix(c(A,B,C,D),nrow=m+m,byrow=TRUE)


dif.is.th1<-(1/(2*theta01^(2)*theta02*(1-rho0^2)))*(rbind(cbind(-2*theta02*diag(m),
                                                                rho0*sqrt(theta01*theta02)*diag(m)),
                                                                cbind(rho0*sqrt(theta01*theta02)*diag(m),
                                                                diag(0,m))))

dif.is.th2<-(1/(2*theta01*theta02^(2)*(1-rho0^2)))*(rbind(cbind(diag(0,m),
                                                                rho0*sqrt(theta01*theta02)*diag(m)),
                                                          cbind(rho0*sqrt(theta01*theta02)*diag(m),
                                                                -2*theta01*diag(m))))

dif.is.rho<-(1/(theta01*theta02*(1-rho0^2)^2))*(rbind(cbind(2*rho0*theta02*diag(m),
                                                            -(1+rho0^(2))*sqrt(theta01*theta02)*diag(m)),
                                                      cbind(-(1+rho0^(2))*sqrt(theta01*theta02)*diag(m),
                                                            2*rho0*theta01*diag(m))))

 J1<-H1[(p+p+1):(p+p+m+m),(p+p+1):(p+p+m+m)]%*%dif.is.th1
 J2<-Sigma%*%dif.is.th1
 J3<-H1[(p+p+1):(p+p+m+m),(p+p+1):(p+p+m+m)]%*%dif.is.th2
 J4<-Sigma%*%dif.is.th2
 J5<-H1[(p+p+1):(p+p+m+m),(p+p+1):(p+p+m+m)]%*%dif.is.rho
 J6<-Sigma%*%dif.is.rho

a11<-sum(diag((J1-J2)%*%(J1-J2)))
a12<-sum(diag(J1%*%J3+J2%*%J4-2*J1%*%J4))
a13<-sum(diag(J1%*%J5+J2%*%J6-2*J1%*%J6))
a22<-sum(diag((J3-J4)%*%(J3-J4)))
a23<-sum(diag(J3%*%J5+J4%*%J6-2*J3%*%J6))
a33<-sum(diag((J5-J6)%*%(J5-J6)))

varmat<-2*solve(matrix(c(a11,a12,a13,a12,a22,a23,a13,a23,a33),ncol=3))
    
    se.var.par<-sqrt(diag(varmat))

stdvar<-cbind(c(theta1=theta1,theta2=theta2,rho=rho),
              sqrt(diag(varmat)),
              2*(1-pnorm(abs(c(theta1=theta1,theta2=theta2,rho=rho)/sqrt(diag(varmat))))))
dimnames(stdvar)<-list(c("theta1","theta2","rho"),c("estimate","s.e.","p-value"))
stdvar<-round(stdvar,3)

# output results
    
    beta<-par.lat[1:p]
    se.beta<-sqrt(diag(H1)[1:p])
 
   gamma<-par.lat[(p+1):(p+p)]
    se.gamma<-sqrt(diag(H1)[(p+1):(p+p)])

    list(eta.R=eta.R,eta.D=eta.D)  

#Beta
ebeta<-cbind(beta,se.beta,2*(1-pnorm(abs(beta/se.beta))))
    names(beta)<-c("Age","Male")
    dimnames(ebeta)<-list(names(beta),c("Estimate","SE","p-value"))
    ebeta<-round(ebeta,3)
    options(digits=3)   
    
    cat('\n',"----------------------------------",'\n')
    
#Gamma
egamma<-cbind(gamma,se.gamma,2*(1-pnorm(abs(gamma/se.gamma))))
    names(gamma)<-c("Age","Male")
    dimnames(egamma)<-list(names(gamma),c("Estimate","SE","p-value"))
    egamma<-round(egamma,3)
    options(digits=3)   
    
    cat('\n',"----------------------------------",'\n')

if(flag==0) stop("not reach the converge")
    
    cat('\n',"Bivariate frailty model for recurrent events and death:\n\n")
    
    cat('\n',"Recurrent events",'\n')
    
   print(ebeta)
    
   cat('\n')
    
    cat('\n',"Survival (time to death)",'\n')
    print(egamma)
    
   cat('\n')
   
   cat('\n',"Frailty variance parameters:",'\n')
    
    print(stdvar[(1:3),])
    
    cat('\n',"----------------------------------",'\n')
    
    
}
joint.frailty(dat,patient=dat[,2],theta01=0.2,theta02=0.2,rho0=0.2,itmax=300)




