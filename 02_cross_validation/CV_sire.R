rm(list=ls())


print("part predict mixed population SIRE")

require("BGLR")
require("parallel")
require("tidyverse")
source("functions/function_CV_sire.R") 

cat("loading all data \n")

setwd("/home/enrico/Work/prediction_micro_v2/")
load("rdata/start2025-11-18.RData")
attach(pass)



m_verbose = FALSE
dfig = 4
dfigr = 3
nIter=100000#/100
burnIn = 20000#/100
thin = 50#/10
catene=round((nIter-burnIn)/thin)


cat("number of chain",catene,"\n")

skip_ebv=TRUE # skippig ebv prediction.


# Crearte whole dataset
FIXED=plyr::rbind.fill(as.data.frame(fixC),as.data.frame(fixP))
FIXED[is.na(FIXED)]=0
gv=c(rep(1,ncol(fixC)),rep(2,ncol(FIXED)-ncol(fixC)))
FIXED=as.matrix(FIXED)


pbcb=readxl::read_excel(("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/models1.xlsx"))
pbcb=pbcb[pbcb$CB%in%names(cb),]
pbcb=pbcb[pbcb$PB%in%names(pb),]

y_pb=pbcb$PB
y_cb=pbcb$CB
#

names(pass)
all=plyr::rbind.fill(cb,pb)
head(all)
row.names(FIXED) = c(all$animGLDS)

run = function(Y,model) return(Multitrait(y=Y,resCov=list(type="DIAG"),ETA=model,nIter = nIter, burnIn = burnIn, thin = thin,verbose = FALSE))

all_SNPx=all_SNP$u%*%diag(all_SNP$d)
all_OTUx=all_OTU$u%*%diag(all_OTU$d)

allmodel=list(
   Whole_Gpred = list(G=list(X=all_SNPx,model="BRR"),F=list(X=FIXED,  model="FIXED",common=FALSE, idColumns=gv)),
   Whole_Mpred = list(M=list(X=all_OTUx,model="BRR"),F=list(X=FIXED,  model="FIXED",common=FALSE, idColumns=gv)),
   Whole_oMpred = list(M=list(X=all_OTUx,model="BRR"),F=list(X=FIXED,  model="FIXED",common=FALSE, idColumns=gv)),
   Whole_GMpred  = list(G=list(X=all_SNPx,model="BRR"),
                        M=list(X=all_OTUx,model="BRR"),F=list(X=FIXED,  model="FIXED",common=FALSE, idColumns=gv))
   ) 
   

run2=function(model) return(run(Y,model))
run2=function(Y) return(run(Y,allmodel[["Whole_Gpred"]]))

#####

all_information=list()

dir.create("CV_sire")
setwd("CV_sire")

getwd()

set.seed(1973)

n=length(allmodel)

#-- parementr of cv ---#
proportion=0.2
n_times=15
#
# create training and test set 
fx = function() sample(sire, size = length(sire) * proportion)
# select  sire
test_s <- replicate(n_times, fx(), simplify = FALSE)
# creta test and trainig
cv_cort=lapply(test_s,test_train)
skip_ebv=TRUE

for(uu in 7)  { #

  cat("prepare data for each phenotype (bit trait model) ...")

   print(y_pb[uu])
   l_pb_fen=y_pb[uu] 
   l_cb_fen=y_cb[uu] 

   print("Prpare the phenotype MTM\n")

   yC=as.matrix(cb[,l_cb_fen]); yC=scale(as.double(yC));yC=as.matrix(yC);row.names(yC) = cb$animGLDS
   yP=as.matrix(pb[,l_pb_fen]);yP=scale(as.double(yP));yP=as.matrix(yP);row.names(yP) = pb$animGLDS
   colnames(yC)="yC";colnames(yP)="yP"
   Y=plyr::rbind.fill(as.data.frame(yC),as.data.frame(yP))
   Y[is.na(Y)]=NA
   Y=as.matrix(Y)
   row.names(Y) = c(row.names(yC),row.names(yP))

   cat(" done\n")
   if(!skip_ebv) {
            cat("equazione whole step 1 ..")
            
            Whole=parallel::mclapply(allmodel, function(x) run(as.matrix(Y),x),mc.cores =n)  
            xc=names(Whole)[2:4]
            
            cat("and extract preidction\n")
            wgdir=lapply(c("Whole_Gpred","Whole_GMpred"),function(x) all_SNPx%*% Whole[[x]]$ETA$G$beta)
            wmy=lapply(xc,function(i) all_OTUx%*%Whole[[i]]$ETA$M$b)
            
            cat("don3\n")

            names(wmy) = xc
            names(wgdir)= c("Whole_Gpred","Whole_GMpred")
            
            cat("equazione whole step 2 mediated ..")
            Whole_eq2=parallel::mclapply(wmy, run2,mc.cores =n)  
            wgm=lapply(xc,function(i) all_SNPx%*%Whole_eq2[[i]]$ETA$G$b)
            
           # Ypass=as.data.frame(cbind( c(rep("cb",length(yC)), rep("pb",length(yP))), c(Y[,1],Y[,2])))
           # Ypass$V2=as.double(Ypass$V2)
   } else {
            cat("skkiping ebv LR CV \n")
   }
   
   # RUN THE CV 
   result=parallel::mclapply(cv_cort, mainCV, mc.cores = 10)
   
   print("finitoooo saving ....")
   save.image(paste0("pred_",Sys.Date(),"_",l_cb_fen,".Rdata"))
}  





#save.image(paste0("fin_",Sys.Date(),"_all.Rdata"))
stop("finito")






stop()

phen=
data=
load(paste0("CV_sire/results",data,y_pb[uu],".RData"))


correlation[3,1]
correlation[4,2]



cor(,,use="complete.obs")

# load data
load("sire/results.RData")
require("tidyverse")
tot=all_information[["PIG_OffTestBF_cms"]][["phen"]]
lr=all_information[["PIG_OffTestBF_cms"]][["phen"]]




head(lr)
tot %>% group_by(V1,V4) %>% summarize(mean=mean(as.double(r1)),sd=sd(r1)) #%>% View()
lr %>% group_by(V4) %>% summarize(cor=mean(as.double(lr1)),slope=mean(as.double(lr2)),bias=mean(as.double(lr3)))
lr_tot %>% group_by(V4) %>% summarize(cor=sd(as.double(lr1)),slope=sd(as.double(lr2)),bias=sd(as.double(lr3)))


ggplot(Ytmp[,c("V2","V1","pred")],aes(x=V2,y=pred,color=V1))+geom_point()
cor(Ytmp[,c("pred","V2")],use="complete.obs")[1,2]

length(portacasa)





#####

tmp=portacasa[[10]]
names(tmp)

[["Y"]][mask,1]



# fai lr 


result <- run_pipeline(Ytmp, allmodel_P, OTUxx, SNP, SNPtrainSVD, OTUi, OTUtrainsSVD, n)



# orginal version
 Partial=parallel::mclapply(allmodel_P, run2,mc.cores = n)  
    names(Partial) = names(allmodel_P)

   # estrai la predizione del microbiota
    xc=names(Partial)[2:n]
    pmy=list();for(i in xc) pmy[[i]] = OTUxx%*%Partial[[i]]$ETA$M$b
    
    run2=function(Y) return(run(Y,allmodel_P[["Partial_Gpred"]]))
    Partial_eq2=parallel::mclapply(pmy, run2,mc.cores =n)  
    
    pgm=list();for(i in names(Partial_eq2)) pgm[[i]] = SNP%*%(SNPtrainSVD$v%*%Partial_eq2[[i]]$ETA$G$b)
    pmy=list();for(i in names(Partial_eq2))  pmy[[i]] = OTUi%*%(OTUtrainsSVD$v%*%Partial[[i]]$ETA$M$b)  

    pdir=list()
    pdir[["Partial_Gpred"]] = SNP%*%(SNPtrainSVD$v%*%Partial[["Partial_Gpred"]]$ETA$G$b)
    pdir[["Partial_GMpred"]] =SNP%*%(SNPtrainSVD$v%*%Partial[["Partial_GMpred"]]$ETA$G$b) " and then if worng skip

}"

}