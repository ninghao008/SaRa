# Instructions  (Pseudo-code)
# Step 1: run all the functions defined below;
# step 2: Obtain y (usually Log R Ratio for SNP genotyping data) from your data set
# step 3: Let n=length(y); choice bandwith h.
# step 4: For simple SaRa, CP =  SARAp(Y=y, h=h)
#         then CP$index gives a candidate set of pontential change point positions
#         use pV = CP$pV; CP$index[order(pV)] to obtain a sequence of pontentail change points ranked by local P-values
#         and use sort(pV) to see the local P-values for the sequence
# step 5: FDR approach.
# step 5-1: get an empirical distribution of local P-values
#          FINV    = fInverse(n=n, h= h)     # this function is based on the normaility assumption of the noise
# step 5-2: changeP = SARAFDR(Y=y, h=h, FINV = FINV, fdr=c(0.05,0.1,0.15))
#           give a location set of change points selected by FDRSaRa

##############################################################################
# 1) Local diagnostic function D_h
##############################################################################
llce <- function(y,h) {
    yy = c(rep(0,h-1),y,rep(0,h))            # add some zeros on the head and tail
    n = length(y)
    z = rep(0,n)
    for (i in 1:h) {z = z + yy[i:(n+i-1)]-yy[(h+i):(n+h+i-1)]}
    return(z/h)                                                                 #
}                                                                               #


#############################################################################
#2) Local Max
#localMax(x,span=5)
#return a set of indices which are local max within an area of length 2*span+1
#############################################################################
localMax = function(y,span=5){
    if (length(y)<span*2+1)  return(NULL)
    n  = length(y)
    index = NULL
    for (i in (span+1):(n-span) ) {
         if ( y[i]==max(y[(i-span):(i+span)]) ) index=c(index,i)
    }
    return(index)
}

#############################################################################
#3)estimate sd
#############################################################################
estimateSigma = function(Y,h=10){                   #constant case
  n     = length(Y)                                 # can we make it faster?
  YBar  = rep(0,n)
  for (i in 1:n) {
     a       = min(n,i+h)
     b       = max(1,i-h)
     YBar[i] = mean(Y[b:a])
  }
  return(sd(Y-YBar))
}
#4)Distance from one vector to another
#compare two vectors: calculate 'distance' defined by
#compareVec = max_{x\in vec1}(min_{y\in vec2}(abs(x-y)))
#compare two vectors: calculate average 'distance' defined by
#compareVec2 = mean_{x\in vec1}(min_{y\in vec2}(abs(x-y)))
#compare two vectors: calculate median 'distance' defined by
#compareVec3 = median_{x\in vec1}(min_{y\in vec2}(abs(x-y)))

compareVec   = function(vec1,vec2){
	l1     = length(vec1)
	l2     = length(vec2)
	if ( (l2 == 0) | (l1 == 0) ) return (10000)
	dist   = min(abs(vec2-vec1[1]))
      if (l1>1) {
	    for (i in 2:l1){	dist = max(dist,min(abs(vec2-vec1[i]))) }
	}
	return(dist)
}

compareVec2  = function(vec1,vec2){
	l1     = length(vec1)
	l2     = length(vec2)
	if ( (l2 == 0) | (l1 == 0) ) return (10000)
	dist   = min(abs(vec2-vec1[1]))
      if (l1>1) {
	for (i in 2:l1){ 	dist = sum(dist,min(abs(vec2-vec1[i]))) }
	}
	return(dist/l1)
}

compareVec3  = function(vec1,vec2){
	l1     = length(vec1)
	l2     = length(vec2)
	if ( (l2 == 0) | (l1 == 0) ) return (10000)
	dist   = NULL
	for (i in 1:l1){
		dist = c(dist,min(abs(vec2-vec1[i])))
	}
      mddist = median(dist)
	return(mddist)
}
#############################################################################
#5) distribution of local min of p-value                                   #
#############################################################################

fInverse = function(n=10000, h= 20, hh=20, precise=10000, simT=100){        #need to be faster
   set.seed(2011)
   empirical = NULL
   for (i in 1:simT){
     Y   = rnorm(n)
     LDF  =  abs(llce(Y,h))
     sigma =1
     index = localMax(y=LDF,span=hh)
     pV   = 2*(1-pnorm(LDF[index]/(sqrt(2/h)*sigma)))
     empirical = c(empirical,pV)
     if (length(empirical)>100000) break
   }
   return(quantile(empirical,probs=c(0:precise)/precise))
}

############################################################################
#6) simple SARA  p-value returned      #    # return index and local min of p-value
############################################################################
SARAp = function(Y, h, hh, sigma=NULL){
   n        = length(Y)
   LDF      =  abs(llce(Y,h))
   if (is.null(sigma)) sigma = estimateSigma(Y, h=max(3,2*floor(log(n))))
   pV       = 2*(1-pnorm(LDF/(sqrt(2/h)*sigma)))
   LocalMin = localMax(LDF,span=hh)                                      # aviod the case that several P-values are zero
   #LocalMin = localMax(-pV,span=hh)                                     # the index set of all local  Min
   LocalMinValue = pV[LocalMin]                                       # all local min value

	 return(pV=list(index=LocalMin,pV=LocalMinValue))
}



############################################################################
#7) FDR control                        #
############################################################################
SARAFDR = function(Y, h=10, hh=10, sigma=NULL, FINV = NULL, precise=10000, fdr=c(0.05,0.1,0.15)){
   object = SARAp(Y=Y, h=h, hh=hh, sigma=sigma)
   index  = object$index
   pV     = object$pV
   if (is.null(FINV)) FINV   = fInverse(n=length(Y), h= h, hh=hh, precise=precise, simT=100)  # could not afford to calculate too many times.
   pVcorr = pV
   for (i in 1:length(pV)){
         #pVcorr[i] = (max(length(which(FINV<pV[i])),1)-1)/precise
         pVcorr[i] = (length(which(FINV<pV[i])))/(precise+1)
   }
   plot(sort(pVcorr))
   pVcorrS=sort(pVcorr)
   K      = length(pVcorr)
   thresh =rep(0,length(fdr))
   changeP = list()
   for (i in 1:length(fdr)){
      #thresh[i]= pVcorrS[max(which(pVcorrS-c(1:K)*fdr[i]/K <=0))]
      aaa = which(pVcorrS-c(1:K)*fdr[i]/K <=0)
      if (length(aaa)==0) {thresh[i]=0} else {thresh[i]= pVcorrS[max(aaa)]}
      if (length(which(pVcorr<=thresh[i]))>0) {changeP[[i]] = index[pVcorr<=thresh[i]]}
        else   changeP[[i]] = 0
   }
   return(changeP)
}


#############################################################################
#simulation I
#############################################################################
#model generation

#################################################################################
set.seed(2011)                                                                  #
n = 30000                                                                       #
J = 50                                                                          #
sigma  = 1                                                                      #
fdr=c(0.05,0.1,0.15)                                                            #
                                                                                #
a= unique(c(0,floor(runif(50)*n/5)*5,n))       # random change-point + 0 and n  #
while (length(a)<J+2) {a=unique(c(a, floor(runif(1)*n/5)*5))}                   #
a = sort(a)                                                                     #
                                                                                #
#> a                                                                            #
# [1]     0   650   855  1260  1600  2835  3845  4105  4225  5510  6160  6315   #
#[13]  7055  7150  7685  7860  9425 10110 10390 11070 11085 11935 12050 12690   #
#[25] 12810 14655 14970 15725 17095 17215 17610 17650 17780 17990 19070 20930   #
#[37] 21060 21170 21365 22515 22800 23380 24090 24295 24450 25645 26295 26895   #
#[49] 27160 28365 29630 30000                                                   #
                                                                                #
a1=a[-(J+2)]                                                                    #
a2=a[-1]                                                                        #
aDiff=a2-a1                                                                     #
Y0 = NULL                                                                       #
for (i in 0:J) { Y0 = c(Y0,rep(i%%2, aDiff[i+1]))}                              #
#################################################################################

#Delta = 1.5
#Y0 = Y0*Delta
#h = 20
#hh = 20
#precise = 1000
#FINV   = fInverse(n=n, h= h, hh=hh, precise=precise, simT=100)

simulI = function(Y0, sigma = 1, Delta = 1.5, h = 20, hh = 20, precise = 1000) {
set.seed(85715)
Y1      = Y0*Delta
n       = length(Y0)
FINV    = fInverse(n=n, h= h, hh=hh, precise=precise, simT=100)
CCount1 = NULL
CCount2 = NULL
CCount3 = NULL
TP1 = NULL
TP2 = NULL
TP3 = NULL
FDR1 = NULL
FDR2 = NULL
FDR3 = NULL
for (ii in 1:100){

   Y = Y1+rnorm(n,sd=sigma)
   changeP = SARAFDR(Y, h=h, hh=hh, sigma=NULL, FINV = FINV, precise=precise, fdr=c(0.05,0.1,0.15))
   f1=f2=f3=0
   for (i in 2:{J+1}){
      if (compareVec(a[i],changeP[[1]])<11) f1=f1+1
      if (compareVec(a[i],changeP[[2]])<11) f2=f2+1
      if (compareVec(a[i],changeP[[3]])<11) f3=f3+1                   # 11 can be changed
   }

   CCount1 = c(CCount1,length(changeP[[1]]))
   CCount2 = c(CCount2,length(changeP[[2]]))
   CCount3 = c(CCount3,length(changeP[[3]]))
   TP1     = c(TP1,f1)
   TP2     = c(TP2,f2)
   TP3     = c(TP3,f3)

   fdc     = 0
   #print(changeP[[1]])
   for ( i in 1:length(changeP[[1]])){
      if (compareVec(changeP[[1]][i],a)>10) fdc=fdc+1                 #10 can be changed
   }
   FDR1 = c(FDR1,fdc/length(changeP[[1]]))

   fdc     = 0
   for ( i in 1:length(changeP[[2]])){
      if (compareVec(changeP[[2]][i],a)>10) fdc=fdc+1
   }
   FDR2 = c(FDR2,fdc/length(changeP[[2]]))

   fdc     = 0
   for ( i in 1:length(changeP[[3]])){
      if (compareVec(changeP[[3]][i],a)>10) fdc=fdc+1
   }
   FDR3 = c(FDR3,fdc/length(changeP[[3]]))
}

return (list(CCount1=mean(CCount1), CCount2=mean(CCount2), CCount3=mean(CCount3),
 TP1 = mean(TP1), TP2 = mean(TP2), TP3 = mean(TP3),
 FDR1 = mean(FDR1), FDR2 = mean(FDR2), FDR3 = mean(FDR3)))
}

 
A1 = simulI(Y0=Y0, sigma = 1,Delta = 1.5, h = 10, hh = 10, precise = 1000)
A2 = simulI(Y0=Y0, sigma = 1,Delta = 1.5, h = 20, hh = 20, precise = 1000)
A3 = simulI(Y0=Y0, sigma = 1,Delta = 1.5, h = 30, hh = 30, precise = 1000)

 
B1 = simulI(Y0=Y0, sigma = 1,Delta = 3, h = 10, hh = 10, precise = 1000)
B2 = simulI(Y0=Y0, sigma = 1,Delta = 3, h = 20, hh = 20, precise = 1000)
B3 = simulI(Y0=Y0, sigma = 1,Delta = 3, h = 30, hh = 30, precise = 1000)

save(A1,A2,A3,B1,B2,B3,file='simulationI.Rdata')
load('simulationI.Rdata')
###############################################################################
#table
kkk=c(1,4,7,2,5,8,3,6,9)
Table =NULL
Row   = NULL
for (i in kkk) {Row=c(Row,A1[[i]])}
Table = rbind(Table, Row); Row= NULL
for (i in kkk) {Row=c(Row,A2[[i]])}
Table = rbind(Table, Row); Row= NULL
for (i in kkk) {Row=c(Row,A3[[i]])}
Table = rbind(Table, Row); Row= NULL
for (i in kkk) {Row=c(Row,B1[[i]])}
Table = rbind(Table, Row); Row= NULL
for (i in kkk) {Row=c(Row,B2[[i]])}
Table = rbind(Table, Row); Row= NULL
for (i in kkk) {Row=c(Row,B3[[i]])}
Table = rbind(Table, Row);
library(xtable)
xtable(Table,digit=3)


 
