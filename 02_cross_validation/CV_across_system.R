cat("CROSS BREED POPULATION PREDICT PUREBREED")

rm(list=ls())
source('functions/function_cv.R')  

setwd("/home/enrico/Work/prediction_micro_v2/")

require("BGLR")
require("parallel")
require("tidyverse")

# vedere se con questi si possono
#source("script/swine_pred/swine_def/other_function.R")
#source("script/swine_pred/swine_def/exctract_results.R")

system("ls rdata/*")



load("rdata/start2025-11-18.RData")
attach(pass)
ls()

  

names(pass)

# Paremtro inziale 

m_verbose = FALSE
dfig = 2
dfigr = 2
nIter=100000
burnIn = 50000
thin = 50
catene=round((nIter-burnIn)/thin)



cat("number of chain",catene,"\n")

cat("to accelarate create matrix UD of SVD\n")

SVDG_c=UGC%*%DGC
SVDM_c=UMC%*%DMC

SVDM_p=UMP%*%DMP
SVDG_p=UGP%*%DGP

print(SVDG_p[1:10,1:10])
cat("prepare database for MTM\n")

FIX=plyr::rbind.fill(as.data.frame(fixC),as.data.frame(fixP))
FIX[is.na(FIX)]=0  #occhio alla media 
FIX=FIX[,-1]
gv=c(rep(1,ncol(fixC)-1),rep(2,ncol(fixP)-1))
length(gv)==ncol(FIX[])
   
all_SNPx=all_SNP$u%*%diag(all_SNP$d)  
all_OTUx=all_OTU$u%*%diag(all_OTU$d)

print(SVDG_c[1:10,1:10])

cat("Prepare all models \n")

# GPred Preidction 
# CB
Gp_CB = list( F=list(X =   fixC, df = dfig, model = "FIXED"),
               G = list(X = SVDG_c, df = dfig, model = "BRR") )
# PB
Gp_PB = list(   F=list(X =   fixP, df = dfig, model = "FIXED"),
               G = list(X = SVDG_p, df = dfig, model = "BRR") )
# CB+PB
G_W = list(G=list(X=all_SNPx,model="BRR"),
            F=list(X=as.matrix(FIX),model="FIXED",common=FALSE, idColumns=gv))

# oMPred Preidction 
#CB
oM_CB = list(M = list(X = SVDM_c, df = dfig, model = "BRR"))
# PB
oM_PB = list( M = list(X = SVDM_p, df = dfig, model = "BRR") )
# CB+PB
oM_W  = list(M=list(X=all_OTUx,model="BRR"))

# MPred Preidction 
#CB
M_CB = list(F=list(X =fixC, model = "FIXED"),
            M = list(X = SVDM_c, df = dfig, model = "BRR") )
#PB
M_PB = list(F=list(X =fixP, model = "FIXED"),
            M = list(X = SVDM_p, df = dfig, model = "BRR") )
#CB+PB
M_W  = list(F=list(X=as.matrix(FIX),model="FIXED",common=FALSE, idColumns=gv),
            M=list(X=all_OTUx, model="BRR"))


#GMPred Preidction 
#cb
GM_CB = list(G = list(X = SVDG_c, df = dfig, model = "BRR") ,
             F=list(X =fixC, df = dfig, model = "FIXED"),
             M = list(X = SVDM_c, df = dfig, model = "BRR") )
#PB
GM_PB = list(G = list(X = SVDG_p, df = dfig, model = "BRR"),
             F=list(X =fixP, model = "FIXED"),
             M = list(X = SVDM_p, df = dfig, model = "BRR") )

#CB+PB
GM_W  = list(G=list(X=all_SNPx,model="BRR"),
             F=list(X=as.matrix(FIX),model="FIXED",common=FALSE, idColumns=gv),
             M=list(X=all_OTUx,model="BRR"))

# put all models in a function

ETA_CB=list("Gpred"=Gp_CB,"Mpred"=M_CB,"oMpred"=oM_CB,"GMpred"=GM_CB)
ETA_PB=list("Gpred"=Gp_PB,"Mpred"=M_PB,"oMpred"=oM_PB,"GMpred"=GM_PB)
ETA_W =list("Gpred"=G_W,"Mpred"=M_W,"oMpred"=oM_W,"GMpred"=GM_W)
NESTEDw=list("Gpred"=G_W,"Mpred"=G_W,"oMpred"=G_W,"GMpred"=G_W)

pbcb=readxl::read_excel(("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/models1.xlsx"))
pbcb=pbcb[pbcb$CB%in%names(cb),]
pbcb=pbcb[pbcb$PB%in%names(pb),]

y_pb=pbcb$PB
y_cb=pbcb$CB

match(y_pb,names(pb))

SVD=svd(OTU_cb)
SVDMCB=list(U=UMC,D=DMC,V=VMC)

getwd()


dir.create("CV_multi")
setwd("CV_multi")

getwd()
# function per fare il bitraits models
#setwd("..")
result=list()
uu=7
GenomicPred =  list()
PhenoPred  =  list()

# CAPIRE CHE LLOOOP FARE
phenos=1:7
paste0(y_cb[phenos]," => ",y_pb[phenos])


for (uu in phenos) {

 # la tre è un alternativa alla 2.
   cat("prepare data and phenotype ", y_pb[uu] ,"\n")
# solo per essere sicuri fenotipo CB
   yC=ready_pheno(cb,y_cb[uu])
   yP=ready_pheno(pb,y_pb[uu])
   
   colnames(yC)="yC";colnames(yP)="yP"
#  Crea un fenotipico uno
   Y=as.matrix(plyr::rbind.fill(as.data.frame(yC),as.data.frame(yP)))
   Y[is.na(Y)]=NA
   row.names(Y) = c(row.names(yC),row.names(yP))
   # pensare se fare ognuno in unsa sua caetwlla
   n=length(ETA_W)
  
   cat("########### PREDICT WITH TRAINING DATA ########### \n")

   cat("Step 1) solve fist equations number 1\n")
   Partial = parallel::mclapply(ETA_CB,run, mc.cores = n) 
   
   cat("Step 2) extract M effect as speudo fenotype and refit the model ") 
   myCB=lapply(1:4,function(x) obt_predictor(x,Partial,SVDM_c,"M",yC))     
   Partial_eq2 = parallel::mclapply(myCB,run2, mc.cores = n)

   cat("Step 3) generate 3 predictor\n")
   gm_pb=parallel::mclapply(1:4,function(x) projection_svd(x,Partial_eq2,SNP_pb,VGC,"G",yP),mc.cores=n)
   my_pb=parallel::mclapply(1:4,function(x) projection_svd(x,Partial,OTU_pb,VMC,"M",yP),mc.cores=n)
   gy_pb=parallel::mclapply(1:4,function(x) projection_svd(x,Partial,SNP_pb,VGC,"G",yP),mc.cores=n)
   
   cat(" combine 3 predictor to obtain G/M deinfirirve\n")
   partial_g=lapply(1:4,function(x) gm_pb[[x]]+gy_pb[[x]])
   partial_my=lapply(1:4,function(x) my_pb[[x]]+gy_pb[[x]])


   cat("make alternative prediction in that case breeding values estimated directly in the model\n") 
   Partial_eq2V2 = parallel::mclapply(my_pb,run2PB, mc.cores = n)
   gm_pb2=lapply(1:4,function(x) obt_predictor(x,Partial_eq2V2,SVDG_p,"G",yP))
   
   for(i in 1:4) row.names(gm_pb2[[i]])=row.names(yP)
    
   partial_g2=lapply(1:4,function(x) gm_pb2[[x]]+gy_pb[[x]])
   # fai idea mia 
   residual = lapply(1:4,function(i) Partial_eq2V2[[i]][["y"]]-Partial_eq2V2[[i]][["yHat"]] ) 
   
   partial_predv2=lapply(1:4,function(x) residual[[x]]+partial_g2[[x]])
   
   #cor(merge(partial_g2[[ii]],partial_g[[ii]],by="row.names")[,-1])     
   
   if(F) {
         cat("##########    MAKE PREDICTION WITH WHOLE DATA #########\n")
         cat("               slowest part                              ")
         
         cat("Step 1\n")
         Whole = parallel::mclapply(ETA_W,MTMrun, mc.cores = n)
         cat("Step 2a)")
         Whole_my=lapply(1:4,function(x) obt_predictor(x,Whole,all_OTUx,"M",Y))
         Whole_eq2 = parallel::mclapply(Whole_my,MTMrun2, mc.cores = n)
         
         cat("Step3a) generate g  predictor\n")
         Whole_gy=lapply(1:4,function(x) obt_predictor(x,Whole,all_SNPx,"G",Y))
         Whole_gm=lapply(1:4,function(x) obt_predictor(x,Whole_eq2,all_SNPx,"G",Y))
         
         Whole_g = lapply(1:4,function(s) Whole_gy[[s]]+Whole_gm[[s]])
   }

   cat("##########   MAKE PREDICTION WITH PB DATA ############\n")
  
   Pb_only = parallel::mclapply(ETA_PB,runPB, mc.cores = n)
   Pb_only_my=lapply(1:4,function(x) obt_predictor(x,Pb_only,SVDM_p,"M",yP))   
   Pb_eq2 = parallel::mclapply(Pb_only_my,run2PB, mc.cores = n)
   
   Pb_gy=lapply(1:4,function(x) obt_predictor(x,Pb_only,SVDG_p,"G",yP))
   Pb_gm=lapply(1:4,function(x) obt_predictor(x,Pb_eq2,SVDG_p,"G",yP))

   Pb_g=lapply(1:4,function(x) Pb_gm[[x]]+Pb_gy[[x]])
   
   cat("saving imagine..\n")
   save.image(paste0("pred_",Sys.Date(),"_",y_pb[[uu]],".Rdata"))

   cat("#################   FINE   ##################\n")
   cat("inizio tirare su i dati (tiro su piu risultati poi al massimo scremo)\n")

   # next step bootstsp
   if(FALSE) {
      Gpredtmp=do.call("rbind",lapply(1:4,function(x) maing(x)))
      GenomicPred[[y_pb[uu]]]=Gpredtmp
      save(x=Gpredtmp,file=paste0("Gpred",y_pb[[uu]],".Rdata"))
   }   
   Ppredtmp=do.call("rbind",lapply(1:4,function(x) mainP(x)))

   save(x=Ppredtmp,file=paste0("Ppred",y_pb[[uu]],".Rdata"))


   
   PhenoPred[[y_pb[uu]]]=Ppredtmp

   getwd()
}


save.image(paste0("fin",Sys.Date(),".Rdata"))
print("FINEEEEEEEEEEEEEEE")
stop()

