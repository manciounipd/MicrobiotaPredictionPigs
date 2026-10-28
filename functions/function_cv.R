
##########################################################################
#   FUNZIONE PER I PREDIREEE
############################################################################

getlr=function(comp) {
      x=comp[,1]
      y=comp[,2]
      bias=function(x,y) format((mean(x-y)/var(x)),nsmall=3)
      slope=function(x,y) format(cov(x, y) / var(x) ,nsmall=3)
      acc=function(x,y) format(cor(x,y),nsmall=3)
      return(list(acc=acc(x,y),slope=slope(x,y),bias=bias(x,y)))
}  

ready_pheno = function(db,col) {
                   yC=db[,col]
                   yC=scale(as.double(yC))
                   yC=as.matrix(yC)
                   row.names(yC) = db$animGLDS
                   #names(yC)=col
                   return(yC)
}


Multitraits = function(Y,ETA)  return(Multitrait(y=Y,intercept = TRUE,resCov=list(type="DIAG"),ETA=ETA,nIter = round(nIter*1.5), burnIn = burnIn, thin = thin,verbose=FALSE))
MTMrun=function(ETA) return(Multitraits(Y,ETA))
MTMrun2=function(Y) return(Multitraits(Y,G_W))

rungibbs=function(y,ETA)  return(BGLR(y = y , ETA = ETA, nIter = nIter*2, burnIn = burnIn, thin = thin, df = dfigr, verbose =FALSE))
run=function(ETA) rungibbs(yC,ETA)
run2=function(yC) rungibbs(yC,Gp_CB)

runPB=function(ETA) rungibbs(yP,ETA)
run2PB=function(yP) rungibbs(yP,Gp_PB)

backsolvingSVD = function(V,U,D,y) as.numeric(V %*% (crossprod(U, y) / diag(D)))
backsolvingSVD2=function(x) backsolvingSVD(VMC,UMC,DMC,x)
   
corr=function(x,db) cor(yP,db[,x],use = "pairwise.complete.obs")
rmse=function(x,db) Metrics::rmse(yP,db[,x])

projection_svd <- function(x, BGLRc, Xtest, Vtrain, effect,Y) {
    # Determine the length of n (number of samples or rows in the response)
    # Extract the effects based on the specified effect type
    if(!is.null(names(BGLRc[[x]]))) {
    eff <- BGLRc[[x]][["ETA"]][[effect]][["b"]]
    if (is.null(eff)) eff <- BGLRc[[x]][["ETA"]][[effect]][["beta"]]
    
    # If effect is null, handle the case where we compute my as a default
    if (is.null(eff) || all(is.na(eff)) ) {
            my=matrix(0,nrow=dim(Y)[1],ncol=dim(Y)[2])   
        } else {
            my <- Xtest %*% Vtrain %*% eff  # Apply projection if effect exists
        }
          #return(my)  # Return the result of the projection
    } else {
            my=matrix(0,nrow=dim(Y)[1],ncol=dim(Y)[2])
             #return(my)
}

      row.names(my)=row.names(Y)
      return(my)
  
}



obt_predictor = function(x,BGLRc,X,effect,Y) {

            eff = BGLRc[[x]]
            if(length(eff)>1){
               eff=eff[["ETA"]][[effect]]            
               if(any(names(eff) %in% "b"))   eff=eff[["b"]]
               if(any(names(eff) %in% "beta")) eff=eff[["beta"]]
               
               if (is.null(eff) || all(is.na(eff)) ) {
                     my=matrix(0,nrow=dim(Y)[1],ncol=dim(Y)[2])
               } else  {my=X%*% eff}
         
               } else {
                  my=matrix(0,nrow=dim(Y)[1],ncol=dim(Y)[2])
               }
      row.names(my)=row.names(Y)
      return(my)
}


##########################################################################
#   FUNZIONE PER I COLLAZINARE I RISULTATI
############################################################################


merge_by_row = function(x, y) {
  tmp = merge(x, y, by = "row.names")
  row.names(tmp) = tmp[, 1]
  tmp = tmp[, -1]
  return(tmp)
}

merge_by_row_l = function(l) {
  tmp = l[[1]]
  for (i in 2:length(l)) {
    tmp = merge_by_row(tmp, l[[i]])
  }
  return(tmp)
}

corr=function(x) cor(x,use="complete.obs")


mutimerged <- function(df_list,all) {
  if (length(df_list) < 2) return(df_list[[1]])
  
  merged_df <- df_list[[1]]
  
  for (i in 2:length(df_list)) {
    merged_df <- merge(merged_df, df_list[[i]], by = "row.names", all = all)
    
    # After merge, row names are now in Row.names column
    rownames(merged_df) <- merged_df$Row.names
    merged_df$Row.names <- NULL
  }
  
  return(merged_df)
}

biastutti <- function(df, abs = FALSE) {
  cols <- colnames(df)
  n <- length(cols)
  mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      diff <- df[[i]] - df[[j]]
      bias <- mean(diff, na.rm = TRUE)
      mat[i, j] <- if (abs) abs(bias) else bias
    }
  }
  
  rownames(mat) <- cols
  colnames(mat) <- cols
  
  return(mat)
}

slopeall<- function(df) {
  cols <- colnames(df)
  n <- length(cols)
  slope_mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      x <- df[[i]]
      y <- df[[j]]
      valid <- complete.cases(x, y)
      if (sum(valid) == 0) {
            slope_mat[i, j] <- NA
        } else {
            fit <- lm(y ~ x)
            slope <- coef(fit)[2]  # extract slope
            slope_mat[i, j] <- slope
      }
      
    }
  }
  
  rownames(slope_mat) <- cols  # predictor
  colnames(slope_mat) <- cols  # response
  
  return(slope_mat)
}

rmseall <- function(df) {
  cols <- colnames(df)
  n <- length(cols)
  rmse_mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      x <- df[[i]]
      y <- df[[j]]
      valid <- complete.cases(x, y)
        if (sum(valid) == 0) {
            rmse <- NA
        } else {
            rmse <- sqrt(mean((x -y)^2, na.rm = TRUE))
      }
      rmse_mat[i, j] <- rmse
    }
  }
  
  rownames(rmse_mat) <- cols  # predictor
  colnames(rmse_mat) <- cols  # response
  
  return(rmse_mat)
}


corall <- function(df) {
  cols <- colnames(df)
  n <- length(cols)
  cor_mat <- matrix(NA, n, n)
  
  for (i in 1:n) {
    for (j in 1:n) {
      x <- df[[i]]
      y <- df[[j]]
      valid <- complete.cases(x, y)
      
      if (sum(valid) == 0) {
        cor_val <- NA
      } else {
        cor_val <- cor(x[valid], y[valid])
      }
      
      cor_mat[i, j] <- cor_val
    }
  }
  
  rownames(cor_mat) <- cols
  colnames(cor_mat) <- cols
  
  return(cor_mat)
}





getmetrics=function(gpred) {
    a=cbind(metrics="Acc",reshape2::melt(corall(gpred)))
    b=cbind(metrics="Slope",reshape2::melt(slopeall(gpred)))
    d=cbind(metrics="Bias",reshape2::melt(biastutti(gpred)))
    f=cbind(metrics="RMSE",reshape2::melt(rmseall(gpred)))
    return(rbind(a,b,d,f))
}


maing=function(ii) {
   l = list(Pb_g[[ii]], Whole_g[[ii]], partial_g[[ii]], partial_g2[[ii]], gy_pb[[ii]])
   gpred = merge_by_row_l(l)
   names(gpred) = c("Pb_Si","Cb_tmt","Pb_mtm","Pred1","Pred2","gy") 
   gpred=cbind(names(ETA_W)[[ii]],getmetrics(gpred))
   return(gpred)
}

mainP=function(ii) {
   l = list(yP,  partial_g[[ii]], partial_g2[[ii]],partial_my[[ii]],partial_predv2[[ii]])
   Ypred = merge_by_row_l(l)
   names(Ypred) = c("Y","g","g2","my","altmy") 
   gpred=cbind(names(ETA_W)[[ii]],getmetrics(Ypred))
return(gpred)
}
