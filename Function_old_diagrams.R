#Function to plot previous Christen diagrams

#mat is the pam, with species in columns and sites in rows. 
#view is the per-sites view (1) or the per-species view(2)
#limits is 0 to 1 (1) or the range of values of phi (2)

rdp <- function(mat, view, limits) {
  # First, removes empty columns or rows
  mat=as.matrix(mat)
  j=which(colSums(mat)>0)
  mat2=mat[,j]
  i=which(rowSums(mat2)>0)
  mat3=mat2[i,]
  
  if(view==1){
    # View per sites
    n=dim(mat3)[[1]]
    s=dim(mat3)[[2]]
    alfas<<-rowSums(mat3)
    omegas<<-colSums(mat3)
    alfast<<-alfas/s
    omegast<<-omegas/n
    fist<<-mat3%*%omegast
    fistprom<<-fist/alfas

    if(limits==1) {xl=1;yl=1} else {xl=max(fistprom)*1.1;yl=max(alfast)*1.1}
    
    betty<-1/mean(alfast)
    rho<-alfast*(fistprom-1/betty)  
    xM<-seq(1.01/betty,1,length=100)
    yM<-max(rho)/(xM-1/betty)
    xm<-seq(0,.99/betty,length=50)
    ym<-min(rho)/(xM-1/betty)
    
    plot(fistprom,alfast,xlim=c(0,xl),ylim=c(0,yl),
         col = "gray45",
         xlab="Mean normalized dispersion field",ylab="Normalized richness")
    
    abline(v=1/betty, lty = 2)
    abline(h = min(alfast))
    abline(h = max(alfast))
    lines(xM,yM)
    lines(xM,ym)
  } else {
    #View per species
    mat3t=t(mat3)
    nt=dim(mat3t)[[1]]
    st=dim(mat3t)[[2]]
    alfast<-rowSums(mat3t)
    omegast<-colSums(mat3t)
    alfastt<-alfast/st
    omegastt<-omegast/nt
    fistt<-mat3t%*%omegastt
    fistpromt<-fistt/alfast
    
    # Type of limits: 1 -> limits form 0 to 1, not 1, limits near the maximum of the range of x and y.
    if(limits==1) {xl=1;yl=1} else {xl=max(fistpromt)*1.1;yl=max(alfastt)*1.1}
    
    bettyt<-1/mean(alfastt)
    rhot<-alfastt*(fistpromt-1/bettyt)
    xM<-seq(1.01/bettyt,1,length=100)
    yM<-max(rhot)/(xM-1/bettyt)
    xm<-seq(0,.99/bettyt,length=50)
    ym<-min(rhot)/(xm-1/bettyt)
    
    plot(fistpromt,alfastt,xlim=c(0,xl),ylim=c(0,yl),
         col = "gray45",
         xlab="Normalized richness",ylab="Mean normalized dispersion field")
    
    abline(v=1/bettyt)
    lines(xM,yM)
    lines(xM,ym)
  }
}
