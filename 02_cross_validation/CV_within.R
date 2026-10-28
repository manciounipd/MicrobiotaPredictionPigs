rm(list=ls())


print("part predict cage !")
source('functions/function_within.R') 
require("BGLR")
require("parallel")
require("tidyverse")


setwd("/home/enrico/Work/prediction_micro_v2/")
data="2025-11-18" # invece (2025-04-14s) perche contegnono solo 3


load(paste0("rdata/start",data,".RData"))
attach(pass)


m_verbose = FALSE
dfig = 2
dfigr = 2
nIter=70000  #/100
burnIn = 30000 /100
thin = 50 #/10
catene=round((nIter-burnIn)/thin)

cat("number of chain",catene,"\n")

# cv parameters
combinazioni=15
poportion=0.3

#
skip_ebv=TRUE

table(fixP[,2])
table(fixC[,2])




print("to accelarate\n")
SVDG_c=UGC%*%DGC
SVDM_c=UMC%*%DMC

sex=as.matrix(fixC[,1:2])
cage=as.matrix(fixC[,3:ncol(fixC)])

table(sex)

#head(sex)
MDIA=matrix(rep(1,nrow(fixC)))

Gpred  <-  list(
                F=list(X =   fixC[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_c, df = dfig, model = "BRR",saveEffects=TRUE)
                )

oMpred <- list( M = list(X = SVDM_c, df = dfig, model = "BRR",saveEffects=TRUE))

Mpred <-  list(
                F = list(X = fixC[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),             
                M = list(X = SVDM_c, df = dfig, model = "BRR",saveEffects=TRUE))

GMpred <- list( F = list(X = fixC[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_c,  df = dfig, model = "BRR",saveEffects=TRUE), 
                M = list(X = SVDM_c, df = dfig, model = "BRR",saveEffects=TRUE))


Gpred2  <-  list(
                F=list(X =   fixC[,-1],df=dfig,model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDM_c, df=dfig,model = "BRR",saveEffects=TRUE)
                )



ETAsMicro=list("Gpred"=Gpred,"Mpred"=Mpred,"oMpred"=oMpred,"GMpred"=GMpred)

rungibbs_eq2=function(y,ETA) {
     return(BGLR(y = y , ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin,df = dfigr, verbose =FALSE))
}

whole_ebv = function(x) return(rungibbs(yC,ETAsMicro[[x]]))

whole_ebv2 = function(x) return(rungibbs_eq2(myall[[x]],Gpred2))

rungibbs=function(y,ETA)  return(BGLR(y = y , ETA = ETA, nIter = nIter*2, burnIn = burnIn, thin = thin, df = dfigr, verbose =FALSE))



# phenotype / next steps bit-traits
pbcb=readxl::read_excel(("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/models1.xlsx"))
pbcb=pbcb[pbcb$CB%in%names(cb),]
pbcb=pbcb[pbcb$PB%in%names(pb),]
y_pb=pbcb$PB
y_cb=pbcb$CB

# function per fare il bitraits models
SVDMCB=list(U=UMC,D=DMC,V=VMC)

set.seed(1224)


cat("fai le i gruppi test and training\n")

comb <- combn(unique(cb$CIG_Rep2), m = 2, simplify = FALSE)
test_l=list();train_l=list()


for(x in 1:combinazioni) {
tmpall=c()
for(i in unique(cb$CIG_Rep2)) {

        tmp1=sample(cb[cb$CIG_Rep2==i & cb$CIG_Gender=="Male",]$animGLDS)
        tmp1=sample(tmp1,round(poportion*length(tmp1)))
        tmp2=sample(cb[cb$CIG_Rep2==i & cb$CIG_Gender=="Female",]$animGLDS,)
        tmp2=sample(tmp2,round(poportion*length(tmp2)))

        tmpall=c(tmp1,tmp2,tmpall)
        
}

test_l[[x]]=!cb$animGLDS %in% tmpall 
train_l[[x]]=cb$animGLDS %in% tmpall 
}

collect=list()


dir.create("CV_withincage")
setwd("CV_withincage")

getwd()


collectall=list()




for(pehno in 1:7) {

    l_cb_fen=y_cb[pehno] 

    cat("###################################################",l_cb_fen,"\n")
    yC=as.matrix(cb[,l_cb_fen])
    yC=scale(as.double(yC))
    yC=as.matrix(yC)
    row.names(yC) = cb$animGLDS

    if(!skip_ebv) {
          cat("########### EBV WITH WHOLE DATA ########### \n")
          cat("equation 1\n")  
          ebvw=parallel::mclapply(names(ETAsMicro),whole_ebv,mc.cores=4)
          names(ebvw) = names(ETAsMicro)
          cat("solve my\n")
          myall=lapply(names(ETAsMicro)[2:4],function(i) SVDM_c%*%ebvw[[i]]$ETA$M$b) # pseufo phenotype
          names(myall) = names(ETAsMicro)[2:4]   

          cat("solve equation 2 \n ")   
          ebvw2=parallel::mclapply(names(myall),whole_ebv2,mc.cores=3)
          
          cat("crea i whole ebv Mediated and diretto\n")
          wgm=lapply(1:3,function(i) SVDG_c%*%ebvw2[[i]]$ETA$G$b)
          wgy=lapply(c(1,4),function(i) SVDG_c%*%ebvw[[i]]$ETA$G$b)
          
          for(i in 1:3) row.names(wgm[[i]])=row.names(yC)
          for(i in 1:2) row.names(wgy[[i]])=row.names(yC)

          names(wgm) = names(myall);names(wgy) =names(ETAsMicro)[c(1,4)]

          empty_pred=as.data.frame(rep(0,nrow(wgy[[1]])))
          row.names(empty_pred)=cb$animGLDS

          wgm[["Gpred"]] = empty_pred
          wgy[["Mpred"]] = empty_pred
          wgy[["oMpred"]] = empty_pred
          
          
          wg=lapply(names(wgy),function(i) wgm[[i]]+wgy[[i]])
          names(wg)=names(wgy)
   }


   cat(getwd()," ")
   cat(combinazioni,"\n")

   results=parallel::mclapply(1:combinazioni,predict_within,mc.cores=combinazioni)
  #do.call("cbind",results)
   save.image(paste0("pred_",Sys.Date(),"_",y_cb[[pehno]],".Rdata"))
}  



#save.image(paste0("fin_",Sys.Date(),"_all.Rdata"))
stop("finito")

















setwd("..")

getwd()
system("ls CV_withincage/*")
rm(list=ls())

phen="CIG_LiveAverageDailyGain_W2M_kgs"
data="2025-04-22" # qiesyta contitne solo solo microbioma nel tempo 3


cat("reading...")
load(paste0("CV_withincage/pred_",data,"_",phen,".Rdata"))
cat("\n")

db=list()
iter=1
for(ix in 1:15) {
   
tmp=results[[ix]]
test=tmp[[length(tmp)]]

pred=names(tmp)[-length(names(tmp))]

for(predictor in pred ) {
   tmp_tmp=tmp[[predictor]]
   nomi=names(tmp_tmp) 
   dbtmp=data.frame(do.call("rbind",lapply(nomi,function(x) cor(merge(yC,tmp_tmp[[x]],by="row.names")[,-1])[1,2])))
   dbtmp$model=nomi
   dbtmp$pred=predictor
   cat(predictor,"\n")
   db[[iter]]=dbtmp
   iter=iter+1
}
}

require("tidyverse")

db=as.data.frame(do.call("rbind",db))
names(db) = c("cor","model","pred")
db=as.data.frame(db)
row.names(db) = 1:nrow(db)
db$cor=as.double(db$cor)
head(db)

rb <- db %>% filter(!is.na(cor)) %>%
  group_by(pred, model) %>%
  summarize(
    mean = mean(cor),
    sd = sd(cor),
    q025 = quantile(cor, 0.025),
    q50  = quantile(cor, 0.5),
    q975 = quantile(cor, 0.975),
    .groups = "drop"
  )


rb %>%
  filter(model != "Gpred") %>%
  ggplot(aes(x = pred, y = mean, fill = model, label = round(mean, 3))) +
  geom_col(position = position_dodge2()) +
  #facet_grid(phen ~ .) +
  geom_text(position = position_dodge2(width = 0.9), vjust = -0.3) +
  theme_bw()

ls()

# bredingv va bene
names(tmp)

tmp_tmp=g_allv1
iter=1

uu=list()

for(pk in 1:4) {
tmp_tmp=results[[pk]][["gv1"]]
nomi=names(tmp_tmp)


for(i in 1:4) {
   mod=nomi[i]
   print(mod)
   m=merge(tmp_tmp[[mod]],wg[[mod]],by="row.names")[,-1]
   tmpdf=data.frame(do.call("rbind",getlr(m)))
   tmpdf$model=nomi[i]
   tmpdf$metric=c("acc","slope","bias")
   uu[[iter]]=tmpdf
   iter=iter+1
}
#
}

dfx=do.call("rbind",uu)

names(dfx) = c("val","model","metrics")
dfx$val=as.double(dfx$val)
dfx %>% group_by(model,metrics) %>% filter(!is.na(val)) %>% summarize(
    mean = mean(val),
    sd = sd(val),
    q025 = quantile(val, 0.025),
    q50  = quantile(val, 0.5),
    q975 = quantile(val, 0.975),
    .groups = "drop"
  )

   collectall[[l_cb_fen]]=list(partial=results,whole=list("y"=yC,"gy"=wgy,"gmed"=wgm))

}



save(x=collectall,file="results.Rdata")
stop()


load("within_cage/pred_2025-04-21.Rdata")

   # CAPIRE CHE RISULTATI PASSARE #
   
 ###  
   lrall=data.frame()
   for(i in 1:20) {
      df=as.data.frame(do.call("cbind",results[[i]]$lr))
      df$metrics=row.names(df)      
      lrall=rbind(lrall,df)
   }

row.names(lrall)=1:nrow(lrall)
allinf=data.frame()

for(u in unique(lrall$metrics))  {
  fx=function(x) quantile(as.double(x),probs=c(0.025,0.5,0.975))
  info=data.frame(cbind(metrics=u,t(apply(lrall[lrall$metrics==u,1:4],2,fx))))
  allinf=rbind(info,allinf)
}


allcor=data.frame()
nomi=names(results[[1]][[1]])[-ncol(results[[1]][[1]])]
for(i in 1:20) {
  tmp=results[[i]][[1]]
  allcor=rbind(allcor,cor(tmp)[ncol(tmp),-ncol(tmp)])
}

names(allcor)=nomi
apply(allcor,2,fx)

getwd()
save.image("results.RData")

