# correct data 
rm(list=ls())
# questo scipt fa la parte di matematica dove mi scompone le
# matrici ecc ..
# carirca le prime informazini
setwd("/home/enrico/Work/prediction_micro_v2/")

load("rdata/info_r1.R")

print("questo scirpt ci da tutte le info necessario per le analisi come la creazione delle matricole e anche la loro SVD")

# attach all inof
attach(pass_info_1)
names(pass_info_1)


ords=match(c(cb$animGLDS,pb$animGLDS),row.names(K))
Ktmp=K[ords,ords]

ords=match(c(cb$animGLDS,pb$animGLDS),row.names(Km))
KMtmp=Km[ords,ords]


gE=G[match(row.names(Ktmp),row.names(G)),]
gE=gE[,colSums(is.na(gE)) == 0]

SNP_cb=gE[match(cb$animGLDS,row.names(gE)),]
SNP_pb=gE[match(pb$animGLDS,row.names(gE)),]

SNP_pb=SNP_pb[,!apply(SNP_pb,2,sum)==0]
#SNP_pb=SNP_pb[,!apply(SNP_pb,2,sd)==0]
dim(SNP_pb)
SNP_cb=SNP_cb[,colnames(SNP_cb) %in% colnames(SNP_pb)]

dim(SNP_pb)[2]==dim(SNP_cb)[2]


mE=M[match(row.names(KMtmp),row.names(M)),]
mE=mE[,colSums(is.na(mE)) == 0]

OTU_cb=mE[match(cb$animGLDS,row.names(mE)),]
OTU_pb=mE[match(pb$animGLDS,row.names(mE)),]


cat("svd cb G\n")
SVD=svd(SNP_cb)

UGC=SVD$u
DGC=diag(SVD$d)
VGC=SVD$v

cat("svd Pb G\n")
SVD=svd(SNP_pb)

UGP=SVD$u
DGP=diag(SVD$d)
VGP=SVD$v

cat("svd cb M\n")
SVD=svd(OTU_cb)

UMC=SVD$u
DMC=diag(SVD$d)
VMC=SVD$v

cat("svd pb M\n")
SVD=svd(OTU_pb)

UMP=SVD$u
DMP=diag(SVD$d)
VMP=SVD$v

print("attach all SNPs e make SVD")
SNP=rbind(SNP_cb,SNP_pb)
all_SNP=svd(SNP)

print("attach all OTUs e make SVD")
OTU=rbind(OTU_cb,OTU_pb)
all_OTU=svd(OTU)


f1=formula("~CIG_Gender + CIG_Rep2")
fixC <- model.matrix(f1, data = cb)
row.names(fixC)=cb$animGLDS
head(fixC)
fixC=scale(fixC)
fixC[,1]=1

f1=formula("~CIG_Pen")
randC <- model.matrix(f1, data = cb)[, -1]
row.names(randC)=cb$animGLDS
randC=scale(randC)

cat("prepare formula and matrix cross breed \n")

f1=formula("~PIG_Gender + PIG_BirthCG")
fixP <- model.matrix(f1, data = pb)
row.names(fixP)=pb$animGLDS
fixP=scale(fixP)
fixP[,1]=1

f1=formula("~PIG_BirthLitterID")
randP <- model.matrix(f1, data = pb)[, -1]
row.names(randP)=pb$animGLDS
randP=scale(randP)


#SVD=svd(OTU_pb)
if(FALSE){
library(matlib)
ridge_factor <- 1e-8  # Un piccolo valore positivo
Q=(OTU_cb)
Qinv_ <- inv(t(Q) %*% Q)# + ridge_factor * diag(ncol(Q)))
Qinv=Qinv_ %*% t(Q)
dim(Qinv) # take a whileeeee!!!!
}
#SVD=svd(SNP_pb)



pass=list(cb=cb,pb=pb,
          randC=randC,randP=randP,
          fixC=fixC,fixP=fixP, 
 #         Qinv=Qinv,
          SNP=SNP,OTU=OTU,
          all_SNP=all_SNP,
          all_OTU=all_OTU,
          sire=sire,
          OTU_cb=OTU_cb,OTU_pb=OTU_pb,      
          SNP_cb=SNP_cb,SNP_pb=SNP_pb,
          UMP=UMP,DMP=DMP,VMP=VMP,
          UMC=UMC,DMC=DMC,VMC=VMC,
          UGP=UGP,DGP=DGP,VGP=VGP,
          UGC=UGC,DGC=DGC,VGC=VGC)

# save imagine
save(x=pass,file=paste0("rdata/start",Sys.Date(),".RData"))

cat("donee")