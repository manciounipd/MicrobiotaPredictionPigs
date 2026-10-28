

f90_to_BLGR=function(geno,map) {
  cat("loading ...\n")
  load(geno)
  load(map)
  cat("converting ...\n")
  id=geno[,1]
  gx=list(geno[,3])
  g=function(x) strsplit(x, "",fixed = TRUE)
  gx = do.call("g",gx)
  gm=as.matrix(as.data.frame(gx))
  gm[1:10,1:10]
  dim(gm);dim(map)
  gm=t(gm)
  gm=apply(gm,2,as.integer)
  colnames(gm) = map[,1]
  rownames(gm) = id
  gm[gm==5] = NA
  return(gm)  
}


qcG = function(G,maf_t,msnp_t,mind_t) {
    
    Fmaf=function(x) {
      H=na.omit(x)
      A1=(length(H[H==0])*2+length(H[H==1]))/(length(H)*2)
      A2=(length(H[H==2])*2+length(H[H==1]))/(length(H)*2)
      return(min(c(A1,A2)))
      } 

    
    maf=apply(G,1,Fmaf)
    
    Gx=G[,maf > maf_t]
    cat("lost by maf",ncol(Gx)-ncol(G),"\n")

    msnp=apply(G,2,function(x)sum(!is.na(x)))/nrow(G)
    Gxx=Gx[,msnp > msnp_t]
    cat("lost by allele msnp",(abs(ncol(Gxx)-ncol(Gx))),"\n")    
    
    mind=apply(G,1,function(x) sum(!is.na(x)))/ncol(G)
    Gxxx=Gxx[mind > mind_t,]
   cat("individaul lost by  msnpid", print(nrow(Gxx)-nrow(Gx)),"\n")   
   return(Gxxx)

} 

impute2 = function(DT) {
    DT=data.table(G)
    means = (colMeans((DT), na.rm = T)) #boh magari sbalgio, ,a va piu veloce
    for (j in 1:ncol(DT))
      {  set(DT, which(is.na(DT[[j]])), j, means[j])}
      return(DT)    
}


Ikernel = function(G) {
  colSdApply <- function(x, ...)apply(X=x, MARGIN=2, FUN=sd, ...)
  cat("imputing missing value..\n")
  Gi=impute2(G) # impute missing values
  Gsx=as.matrix(Gi) # 
  csd=colSdApply(Gsx) # calculate sd
  Gsx = Gsx[,csd>0.01] # remove low sd
  Gsx = as.matrix(scale(Gsx)) #
  row.names(Gsx)=row.names(G)
  cat("make inner kernle..\n")
  GS=crossprod(as.matrix(t(Gsx)))
  return(GS[]/ncol(G))

  }


minfo = function(file,period) {
        
        load(file)
        # count NAs per animal
        Abun=AbundanceTransf[,-c(1:2)]
        missingid=apply(Abun,1,function(x) length(x[is.na(x)]))
        limite=.10*ncol(Abun)
        summary(missingid)
        row.names(Abun) = AbundanceTransf$animGLDS
        u=paste0("S",period)
        keep=grep(paste(u, collapse = "|"), names(Abun), value = TRUE)
        Abun=Abun[,keep]
        Abun=Abun[missingid < limite,]
        Abun[1:10,1:10]
        Abun=scale(Abun)
        Abun[1:10,1:10]
        Km=crossprod(t(Abun))
        Km = Km/ncol(Abun)
        return(list(Km,Abun))        
}

