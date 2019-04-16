dstat <-
function(y,qs=c(1/3,2/3),gamma=1,f=NULL,fscore=NULL,fr=1,alpha=0.05){

  # Check input
  stopifnot((length(fr)==1)&(fr>=0))
  stopifnot(is.vector(gamma)&(length(gamma)==1)&(gamma>=1))
  stopifnot(is.vector(y))
  stopifnot(is.vector(qs)&all((qs>0)&(qs<1)))
  if (!is.null(f)) {
    stopifnot(is.factor(f))
    stopifnot(length(y)==length(f))
  }
  if (!is.null(fscore)) stopifnot(!is.null(f))
  if (!is.null(fscore)) stopifnot(length(unique(f))==length(fscore))
  if (!is.null(fscore)) stopifnot(all(round(fscore)==fscore))

  # Housekeeping
  if (length(qs)>1) qs<-sort(qs)
  if ((!is.null(f))&(!is.null(fscore))){
    scores<-fscore
    names(scores)<-levels(f)
  }
  else if ((!is.null(f))&(is.null(fscore))) {
    scores<-rep(1,length(unique(f)))
    names(scores)<-levels(f)
  }
  else scores<-NULL

  # Define needed functions

  gconv<-function(g1,g2){
    #convolution of two generating functions
    g<-stats::convolve(g1,rev(g2),type="o")
    names(g)<-0:(length(g)-1)
    g
  }

  gmult<-function(g,k){
    #mutliply an integer random variable by k
    stopifnot((length(k)==1)&(k>=1)&(k==round(k)))
    g0<-g[1]
    g<-g[2:length(g)]
    m<-length(g)
    gk<-matrix(0,k-1,m)
    gk<-as.vector(rbind(gk,g))
    gk<-c(g0,gk)
    names(gk)<-0:(length(gk)-1)
    gk
  }

  gcond<-function(g,x){
    #condition on being >=x
    stopifnot(length(x)==1)
    v<-0:(length(g)-1)
    g<-g*(v>=x)
    g<-g/sum(g)
    names(g)<-v
    g
  }

  gcum<-function(g,k=NULL){
    #Cumulative probability.  Pr(Y>=k) if k is not null.
    v<-0:(length(g)-1)
    if (!is.null(k)){
      stopifnot((length(k)==1)&(k>=0)&(k<=(length(g)-1))&(k==round(k)))
      sum(g[v>=k])
    }
    else{
      cum<-rev(cumsum(rev(g)))
      names(cum)<-v
      cum
    }
  }

  gbinomMC<-function(ns,pr,cond=0,mult=1){
    #Generating function for a binomial of sample size ns
    #conditional on being >= cond, then multiplied by mult
    stopifnot((length(mult)==1)&(mult>=0)&(mult==round(mult)))
    stopifnot(length(cond)==1)
    stopifnot((length(pr)==1)&(pr>=0)&(pr<=1))
    stopifnot((length(ns)==1)&(ns>=0)&(ns==round(ns)))
    if ((mult>0)&(ns>0)){
      g<-stats::dbinom(0:ns,ns,pr)
      g<-gcond(g,cond)
      g<-gmult(g,mult)
    }
    else g<-1
    names(g)<-0:(length(g)-1)
    g
  }

  gnull<-function(tab,pr,A,cond=0,mult=1){
    #tab is 2xL with L binomials
    stopifnot(is.matrix(tab)&(2==dim(tab)[1]))
    L<-dim(tab)[2]
    stopifnot((length(A)==L)&all((A==TRUE)|(A==FALSE)))
    stopifnot(!all(A=FALSE))
    stopifnot((length(pr)==1)|(length(pr)==L))
    stopifnot((length(cond)==1)|(length(cond)==L))
    stopifnot((length(mult)==1)|(length(mult)==L))
    if(length(pr)==1) pr<-rep(pr,L)
    if(length(cond)==1) cond<-rep(cond,L)
    if(length(mult)==1) mult<-rep(mult,L)
    first<-TRUE

    for (k in 1:L){
      if (A[k]){
        g0<-gbinomMC(tab[1,k]+tab[2,k],pr[k],cond=cond[k],mult=mult[k])
        if (first) {
          first<-FALSE
          g<-g0
        }
        else g<-gconv(g,g0)
      }
    }
    g
  }

  aquant<-function(y,qs){
    # Computes an absolute quantile factor and associated
    # integer scores
    ay<-abs(y)
    qay<-stats::quantile(ay,qs)
    resf<-cut(ay,unique(c(0,qay,max(ay))),include.lowest=TRUE,ordered_result=TRUE)
    res<-(as.integer(resf)-1)*(ay>0)
    list(factor=resf,score=res)
  }

  tabularversion<-function(tab,gamma=1,fr=0,mult=1,scores=NULL){
    stopifnot(is.matrix(tab)&(2==dim(tab)[1]))
    stopifnot((length(gamma)==1)&(gamma>=1))
    stopifnot((length(fr)==1)&(fr>=0))
    if (is.vector(tab)) tab<-matrix(tab,2,1)
    L<-dim(tab)[2]
    stopifnot((length(mult)==1)|(length(mult)==L))
    if(length(mult)==1) mult<-rep(mult,L)
    ns<-tab[1,]+tab[2,]
    pr<-gamma/(1+gamma)
    cond<-ceiling(ns*pr*fr)
    A<-tab[1,]>=cond
    stat<-sum(mult*tab[1,]*A)
    if (all(!A)) {
      g<-1
      names(g)<-0
    }
    else g<-gnull(tab,pr,A,cond=cond,mult=mult)
    pct<-tab[1,]/(tab[1,]+tab[2,])
    tab<-rbind(tab,pct)
    pval<-round(gcum(g,stat),10)
    probability<-rep(round(pr,3),length(cond))
    tab2<-cbind(t(tab),A,mult,cond,probability)
    tab2[,3]<-round(tab2[,3],3)
    tab2[,2]<-tab2[,1]+tab2[,2]
    colnames(tab2)[1]<-"B"
    colnames(tab2)[2]<-"I"
    colnames(tab2)[3]<-"B/I"
    colnames(tab2)[5]<-"w"
    colnames(tab2)[6]<-"b"
    colnames(tab2)[7]<-"prob"
    tab2<-tab2[,c(1,2,6,4,5,3,7)]
    tab2<-as.data.frame(tab2)
    if (pval<=alpha){
      txt<-paste("H0 is rejected at level ",alpha," in the presence of a bias of at most Gamma=",gamma)
    }
    else {
      txt<-paste("H0 cannot be rejected at level ",alpha," in the presence of a bias of Gamma=",gamma)
    }
    list(T=stat,pval=pval,scores=scores,table=tab2,summary=txt)
  }

  # Begin program

  qu<-aquant(y,qs)
  if (is.null(f)) {
    group<-qu$factor
    mult<-sort(unique(qu$score))
  }
  else if (is.null(fscore)) {
    group<-f:qu$factor
    mult<-rep(sort(unique(qu$score)),length(levels(f)))
  }
  else{
    group<-f:qu$factor
    mult<-as.vector(outer(sort(unique(qu$score)),fscore,"*"))
  }
  tb<-table(y>0,group)[2:1,]
  if (is.vector(tb)) {
    tb<-matrix(tb,2,1)
    colnames(tb)<-levels(group)[1]
  }
  tabularversion(tb,gamma=gamma,fr=fr,mult=mult,scores=scores)
}
