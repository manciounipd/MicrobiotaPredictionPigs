rm(list=ls())


print("part predict cage !")

require("BGLR")
require("parallel")
require("tidyverse")
require("data.table")


data="2025-11-20" # invece (2025-04-14s) perche contegnono solo 3

setwd("/home/enrico/Work/prediction_micro_v2/")
system("ls rdata/*")

load(paste0("rdata/start",data,".RData"))
attach(pass)

#dim(Km)

m_verbose = FALSE
dfig = 2
dfigr = 2
nIter=70000 #/100
burnIn = 30000#/100
thin = 50#/10
catene=round((nIter-burnIn)/thin)

cat("number of chain",catene,"\n")

SVDG_c=UGC%*%DGC
SVDM_c=UMC%*%DMC

sex=as.matrix(fixC[,1:2])
cage=as.matrix(fixC[,3:ncol(fixC)])
MDIA=matrix(rep(1,nrow(fixC)))

cat("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@   SET ALL MODELS  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@")

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
                G = list(X = SVDG_c, df=dfig,model = "BRR",saveEffects=TRUE)
                )



ETAsMicro=list("Gpred"=Gpred,"Mpred"=Mpred,"oMpred"=oMpred,"GMpred"=GMpred)



rungibbs=function(y,ETA) {
     return(BGLR(y = y , ETA = ETA, nIter = nIter, burnIn = burnIn, thin = thin, verbose =FALSE))
}



get_vce <- function(model,PHENO) {

  #cat(">>> Estimating VCE for model:", model_name, "\n")
  # https://github.com/gdlc/BGLR-R/blob/master/inst/md/heritability.md
  # Extract the ETA model and effect names
  
  effect_names <- names(model)

  # 1 Get all solution files
  files <- grep("_b\\.", list.files(), value = TRUE)

  # 2 Initialize containers
  effects <- list()

  # Step 1: Compute predictions for each effect
  for (eff in effect_names) {
    cat("Processing effect:", eff, "\n")
  
    X <- as.matrix(model[[eff]]$X)
    matching_file <- grep(eff, files, value = TRUE)
    cat("  Matched file:", matching_file, "\n")

    # Step 2: Read posterior samples
    if (length(grep("\\.dat$", matching_file)) > 0) {
      beta <- as.matrix(data.table::fread(matching_file))
      beta <- beta[(1 + burnIn / thin):nrow(beta), ]
    } else if (length(grep("\\.bin$", matching_file)) > 0) {
      beta <- readBinMat(matching_file)
    } else {
      stop(paste("No valid file found for effect", eff))
    }

    # Step 3: Compute predicted values per sample
    predicted_values <- lapply(1:nrow(beta), function(i) X %*% as.matrix(beta[i, ]))
    effects[[eff]] <- do.call("cbind", predicted_values)
  }

  # Step 4: Combine all effects to get total prediction yHat
  yhat <- matrix(0, nrow = nrow(effects[[1]]), ncol = ncol(effects[[1]]))
  for (eff in names(effects)) {
    yhat <- yhat + effects[[eff]]
  }

  # Step 5: Calculate residuals and add as a new effect
  residuals <- lapply(1:ncol(yhat), function(i) PHENO - yhat[, i])
  effects[["RES"]] <- do.call("cbind", residuals)
  

  # Step 6: Compute variance for each effect (across columns)
  vce_list <- lapply(names(effects), function(eff) apply(effects[[eff]], 2, function(x) var(x,na.rm=TRUE)))
  names(vce_list) <- names(effects)
  vce_df <- as.data.frame(vce_list)

  return(vce_df)
}





pbcb=readxl::read_excel(("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/models1.xlsx"))
pbcb=pbcb[pbcb$CB%in%names(cb),]
pbcb=pbcb[pbcb$PB%in%names(pb),]
y_pb=pbcb$PB
y_cb=pbcb$CB

#if(FALSE) {

dir.create("VCE_CB")
setwd("VCE_CB")

# phenotype / next steps bit-traits


getwd()
collect1=list()
collect2=list()

y_cb[c(1,2,5,7)]
pehno=5

getwd()
#load("../rdata/all_vceCB.RData") # cosi se ci sono problemi updaro i valori

y_cb[]

for(pehno in c(2,5,7)) {

    l_cb_fen=y_cb[pehno] 

    cat("################## ##########",l_cb_fen,"\n")
    yC=as.matrix(cb[,l_cb_fen])
    yC=scale(as.double(yC))
    yC=as.matrix(yC)
    row.names(yC) = cb$animGLDS
 

    vce_1=list()
    vce_2=list()

    cat(l_cb_fen,"with: ")

    for (model in names(ETAsMicro)) {
    
    MODEL1=paste0("FIRST_Model_",model,l_cb_fen)
    MODEL2=paste0("SECOND_Model_",model,l_cb_fen)
    
    dir.create(MODEL1)
    dir.create(MODEL2)


    setwd(MODEL1)
    getwd()
    cat("#######################",model)

            modello=ETAsMicro[[model]]
            ym=rungibbs(yC,modello)
            vce_1[[model]]=get_vce(modello,yC)        
            
     setwd("..")

        if(length(grep("M",model))>0) {

             cat("second equations \n")
             setwd(MODEL2)
              
                my=ETAsMicro[[model]][["M"]]$X%*%ym[["ETA"]][["M"]]$b
                system("rm *")
                fm2=rungibbs(my,Gpred2)
                vce_2[[model]]=get_vce(Gpred2,my) 
            
            setwd("..")
         } else {
            cat("\n")
         }

    }
     
     collect1[[l_cb_fen]]=vce_1
     collect2[[l_cb_fen]]=vce_2
     
}


names(collect1)


save.image("../rdata/all_vceCB.RData")

setwd("..")


print("START I PUREBREEED !!!!!!!!!!!!!!")
getwd()
dir.create("VCE_PB")
setwd("VCE_PB")
getwd()
collect1=list()
collect2=list()

SNP_pb=SNP_pb[,apply(SNP_cb,2,sd)!=0]
SVD_p=svd(SNP_pb)

UGP=SVD_p$u;VGP=SVD_p$v;DGP=diag(SVD_p$d)

SVDG_p=UGP%*%DGP
SVDM_p=UMP%*%DMP

head(fixP)

Gpred.pb  <-  list(
                F=list(X =   fixP[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_p, df = dfig, model = "BRR",saveEffects=TRUE)
                )

oMpred.pb <- list( M = list(X = SVDM_p, df = dfig, model = "BRR",saveEffects=TRUE))

Mpred.pb <-  list(
                F = list(X = fixP[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),             
                M = list(X = SVDM_p, df = dfig, model = "BRR",saveEffects=TRUE))

GMpred.pb  <- list( F = list(X = fixP[,-1], df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_p,  df = dfig, model = "BRR",saveEffects=TRUE), 
                M = list(X = SVDM_p, df = dfig, model = "BRR",saveEffects=TRUE))


Gpred2.pb  <-  list(
                F=list(X =   fixP[,-1],df=dfig,model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_p, df=dfig,model = "BRR",saveEffects=TRUE)
                )



ETAsMicro.PB=list("Gpred"=Gpred.pb,"Mpred"=Mpred.pb,"oMpred"=oMpred.pb,"GMpred"=GMpred.pb)

pehno=7
y_pb

for(pehno in c(2,5,7)) {

    l_cb_fen=y_pb[pehno] 

    cat("################## ##########",l_cb_fen,"\n")
    
    yP=as.matrix(pb[,l_cb_fen])
    yP=scale(as.double(yP))
    yP=as.matrix(yP)
    row.names(yP) = pb$animGLDS
 

    vce_1=list()
    vce_2=list()

    cat(l_cb_fen,"with => ")

    for (model in names(ETAsMicro.PB)) {
    
    MODEL1=paste0("FIRST_Model_",model,l_cb_fen)
    MODEL2=paste0("SECOND_Model_",model,l_cb_fen)
    
    dir.create(MODEL1)
    dir.create(MODEL2)


    setwd(MODEL1)
    cat(model)

            modello=ETAsMicro.PB[[model]]
            ym=rungibbs(yP,modello)
            vce_1[[model]]=get_vce(modello,yP)        
            
     setwd("..")

        if(length(grep("M",model))>0) {

             cat("second equations \n")
             setwd(MODEL2)
              
                cat("-------->   do second    <--------\n")
                my=ETAsMicro.PB[[model]][["M"]]$X%*%ym[["ETA"]][["M"]]$b
                system("rm *")
                fm2=rungibbs(my,Gpred2.pb)
                vce_2[[model]]=get_vce(Gpred2.pb,my) 
            
            setwd("..")
         } else {
            cat("\n")
         }

    }
     
     collect1[[l_cb_fen]]=vce_1
     collect2[[l_cb_fen]]=vce_2
     
}


save.image("../rdata/all_vcePB.RData")


print("FINITOOOOOOOOOOOOOOOOOO\n")
