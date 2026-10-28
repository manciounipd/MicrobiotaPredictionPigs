

##############################
#  FUNCTION TO EXTRACT THE VALUES
##############################

get_mediate_prj=function(Partial_eq2,SNPtest,SNPtrainSVD) {
    pgm <- list()
    for (i in names(Partial_eq2)) {
      if (!is.null(Partial_eq2[[i]])) {
        pgm[[i]] <- SNPtest %*% (SNPtrainSVD$v %*% Partial_eq2[[i]]$ETA$G$b)
      }
    }
  return(pgm)
  }

get_dir_prj=function(Partial,SNPtest,SNPtrainSVD) {
    pgm <- list()
    for (i in names(Partial)) {
      if (!is.null(Partial[[i]]$ETA$G$b)) {
        pgm[[i]] <- SNPtest %*% (SNPtrainSVD$v %*% Partial[[i]]$ETA$G$b)
      }
    }
  return(pgm)
  }


get_genmed <- function(Partial_eq2v2,SNPtestx) {
  pgm <- list()
  for (i in names(Partial_eq2v2)) {
    if (!is.null(Partial_eq2v2[[i]]$ETA$G)) {
      print(i)
      pgm[[i]] <- SNPtestx %*% Partial_eq2v2[[i]]$ETA$G$b
    }
  }
  return(pgm)
}

get_my <- function(Partial,OTUxx) {
  pgm <- list()
  for (i in names(Partial)) {
    if (!is.null(Partial[[i]]$ETA$M)) {
      print(i)
      pgm[[i]] <- OTUxx %*% Partial[[i]]$ETA$M$b
    }
  }
return(pgm)
}


get_my_prj=function(Partial,OTUtest,OTUtrainSVD) {
    pmy <- list()
    for (i in names(Partial)) {
      cat("model:",i,": ")
      if (!is.null(Partial[[i]]$ETA$M$b)) {
           cat("c'è\n")
          pmy[[i]] <-  OTUtest %*% (OTUtrainSVD$v %*% Partial[[i]]$ETA$M$b)
    }
  }
  return(pmy)
  }

###########################################
#    fineeeee
#########################################À
            
run_pipeline <- function(allmodel_P,Ytmp,GpredTR,OTUtest,OTUtrainSVD,SNPtestx,SNPtrainSVD,OTUxx,SNPtest) {

  library(parallel)  # For mclapply

  cat(".   Run first equation")
  Partial = lapply(allmodel_P, function(x) run(Ytmp,x)) # main equation
  
  cat(". my on trainig and")
  Pmy=get_my(Partial,OTUxx) # stima il my
  
  cat("... and my project on test\n")
  PredPmy=get_my_prj(Partial,OTUtest,OTUtrainSVD) # non giaralo
  
  cat(".  run my on trainig\n")
  Partial_eq2 =  lapply(Pmy,function(x) run(Ytmp,allmodel_P[["Gpred"]]) ) # scompioni my
   
  cat(".  run my on test \n")
  Partial_eq2v2 = lapply(PredPmy, function(x) run(x,GpredTR)) 
  # getresidual
  residula=lapply(c("Mpred","oMpred","GMpred"),function(x) PredPmy[[x]] - Partial_eq2v2[[x]][["ETAHat"]])
  names(residula) = c("Mpred","oMpred","GMpred")
  # extract the results
  gmv1=get_mediate_prj(Partial_eq2,SNPtest,SNPtrainSVD)
  gmv2=get_genmed(Partial_eq2v2,SNPtestx)

  gydir=get_dir_prj(Partial,SNPtest,SNPtrainSVD)

  empty_tofill=gydir[[1]]
  empty_tofill[,]=0

  
  # aggiungi i valori mancanti
  gydir[["Mpred"]]  = empty_tofill
  gydir[["oMpred"]] = empty_tofill
  gmv1[["Gpred"]]   = empty_tofill
  gmv2[["Gpred"]]   = empty_tofill
  residula[["Gpred"]] = empty_tofill
  
  nomi=names(gmv1)
  for(x in nomi) gmv1[[x]] = gmv1[[x]]+gydir[[x]]
  for(x in nomi) gmv2[[x]] = gmv2[[x]]+gydir[[x]]
  for(x in nomi) PredPmy[[x]] = PredPmy[[x]]+gydir[[x]]
  for(x in nomi) residula[[x]] = residula[[x]]+gmv2[[x]]
  

  # Return all results
  return(list(gmv1 = gmv1, gmv2 = gmv2, gydir= gydir,my = PredPmy,gmres2=residula))
}

cat(".")


# split test and train !
test_train=function(test_s) {

    cat("chose test and train")
    
    test <- c(cb[cb$sire %in% test_s, ]$animGLDS, pb[pb$sire %in% test_s, ]$animGLDS)
    training <- c(cb[!cb$sire %in% test_s, ]$animGLDS, pb[!pb$sire %in% test_s, ]$animGLDS)

    maskTest <- all$animGLDS %in% test
    maskTrain <- all$animGLDS %in% training

    # ALMENO DUE PER CAGE
    cat("################################ Removing the test with fewer than 3 individuals per group...\n")

    pbtr <- table(pb[pb$animGLDS %in% training, "PIG_BirthCG"])
    cbtr <- table(cb[cb$animGLDS %in% training, "CIG_Rep2"])
    pbtesr <- table(pb[pb$animGLDS %in% test, "PIG_BirthCG"])
    cbtesr <- table(cb[cb$animGLDS %in% test, "CIG_Rep2"])

    rmeffTR <- c(names(pbtr[pbtr < 3]), names(cbtesr[cbtesr < 3]))
    rmeffTEST <- c(names(pbtesr[pbtesr < 3]), names(cbtesr[cbtesr < 3]))

    rmtes <- c(pb[pb$PIG_BirthCG %in% rmeffTR, ]$animGLDS, cb[cb$CIG_Rep2 %in% rmeffTR, ]$animGLDS)
    rmtr <- c(pb[pb$PIG_BirthCG %in% rmeffTEST, ]$animGLDS, cb[cb$CIG_Rep2 %in% rmeffTEST, ]$animGLDS)

    training <- training[!training %in% rmtes]
    test <- test[!test %in% rmtes]

    return(list(training=training,test=test))
  
}


mainCV=function(cv) {
    
    training=cv[["training"]]
    test=cv[["test"]]
    

    cat("############################  Create fixed matrices \n")

    fixCtrain <- model.matrix(~ CIG_Rep2 + CIG_Gender, data = cb[!cb$animGLDS %in% test, ])
    fixPtrain <- model.matrix(~ PIG_BirthCG, data = pb[!pb$animGLDS %in% test, ])
    
    fixCtest <- model.matrix(~ CIG_Rep2 + CIG_Gender, data = cb[cb$animGLDS %in% test, ])
    fixPtest <- model.matrix(~ PIG_BirthCG, data = pb[pb$animGLDS %in% test, ])

    FIXEDtr <- plyr::rbind.fill(as.data.frame(fixCtrain), as.data.frame(fixPtrain))
    FIXEDtr[is.na(FIXEDtr)] <- 0
    gv_trained <- c(rep(1, ncol(fixCtrain)), rep(2, ncol(fixPtrain) - 1))

    FIXEDtest <- plyr::rbind.fill(as.data.frame(fixCtest), as.data.frame(fixPtest))
    FIXEDtest[is.na(FIXEDtest)] <- 0
    gv_test <- c(rep(1, ncol(fixCtest)), rep(2, ncol(fixPtest) - 1))

    # Prepare phenotypes
    Ytmp <- as.matrix(Y[row.names(Y) %in% training, ])
    Ytest <- as.matrix(Y[row.names(Y) %in% test, ])

    cat(",    ,   Now, making SVD for train SNP...\n")

    SNPtrain <- SNP[row.names(SNP) %in% training, ]
    SNPtrainSVD <- svd(SNPtrain)
    SNPxx <- SNPtrainSVD$u %*% diag(SNPtrainSVD$d)

    cat("And ... making SVD for train OTU...\n")
    OTUi <- pass[["OTU"]]
    OTUtrain <- OTUi[row.names(OTUi) %in% training, ]
    OTUtrainSVD <- svd(OTUtrain)
    OTUxx <- OTUtrainSVD$u %*% diag(OTUtrainSVD$d)

    cat("Now, making SVD for test SNP...\n")
    SNPtest <- SNP[row.names(SNP) %in% test, ]
    SNPtestSVD <- svd(SNPtest)
    SNPtestx <- SNPtestSVD$u %*% diag(SNPtestSVD$d)

    cat("And ... making SVD for test OTU...\n")
    OTUtest <- as.matrix(OTUi[row.names(OTUi) %in% test, ])
    
    cat("############################# cretate all models\n")
    # Ensure all inputs are matrices
    FIXEDtr <- as.matrix(FIXEDtr)
    FIXEDtest <- as.matrix(FIXEDtest)
    SNPxx <- as.matrix(SNPxx)
    SNPtestx <- as.matrix(SNPtestx)

    # Define model structures
    allmodel_P <- list(
        Gpred = list(G = list(X = SNPxx, model = "BRR"), 
                     F = list(X = FIXEDtr, model = "FIXED", common = FALSE, idColumns = gv_trained)),
        Mpred = list(M = list(X = OTUxx, model = "BRR"), 
                     F = list(X = FIXEDtr, model = "FIXED", common = FALSE, idColumns = gv_trained)),
        oMpred = list(M = list(X = OTUxx, model = "BRR")),
        GMpred = list(G = list(X = SNPxx, model = "BRR"), 
                      M = list(X = OTUxx, model = "BRR"), 
                      F = list(X = FIXEDtr, model = "FIXED", common = FALSE, idColumns = gv_trained))
    )

    GpredTR <- list(G = list(X = SNPtestx, model = "BRR"), 
                    F = list(X = FIXEDtest, model = "FIXED", common = FALSE, idColumns = gv_test))

    cat("Running main models...\n")
    
    tmp=run_pipeline(allmodel_P,Ytmp,GpredTR,OTUtest,OTUtrainSVD,SNPtestx,SNPtrainSVD,OTUxx,SNPtest)
    
    return(tmp)
}

