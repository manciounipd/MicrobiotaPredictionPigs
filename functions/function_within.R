
#combinazini=1
predict_within=function(combinazini) {
   
   test=test_l[[combinazini]]
   train=train_l[[combinazini]]
   
   
   cat("perpare the matrix:")
   cat("         do test\n")

   yCtmp=yC[train]
   yCtmp=as.matrix(yCtmp)
   row.names(yCtmp)=cb$animGLDS[train]

   yCtest=yC[test]
   yCtest=as.matrix(yCtest)
   row.names(yCtest)=cb$animGLDS[test]

   cat("          do trainig\n")
   OTU_cbtmp=OTU_cb[train,]
   OTU_cbTEST=OTU_cb[test,]

   SNP_cbtmp=SNP_cb[train,]
   SNP_cbTEST=SNP_cb[test,]

   cat("making svd....\n")

   SVD_tmp=svd(OTU_cbtmp)
   SVD_test=svd(OTU_cbTEST)

   SNPSVD_tmp=svd(SNP_cbtmp)
   SNPSVD_test=svd(SNP_cbTEST)

   cat("split populations \n")

   cb_tmp=cb[train,]
   cb_test=cb[test,]
    
   nomitest=row.names(yC)[test]
   nomitrain=row.names(yC)[train]

   cat("create test and train matrix #####\n")

   f1=formula("~CIG_Gender + CIG_Rep2")
   fixC <- model.matrix(f1, data = cb_tmp)
   row.names(fixC)=cb_tmp$animGLDS
   fixC[,1]=1

########################################
   f1=formula("~CIG_Gender + CIG_Rep2")
   fixCT <- model.matrix(f1, data = cb_test)
   row.names(fixCT)=cb_test$animGLDS
   fixCT[,1]=1

###################################
      SVDM_ctmp=SVD_tmp$u%*%diag(SVD_tmp$d) # MICROIOTA
      SVDG_ctmp=SNPSVD_tmp$u%*%diag(SNPSVD_tmp$d) # SNP
##   
      SVDG_ctest=SNPSVD_test$u%*%diag(SNPSVD_test$d) # SNP TEST
      SVD_test=SVD_test$u%*%diag(SVD_test$d)    #MICROBIOTA TEST


###############################

   cat("run  base..\n")

   Gpred  <-  list(
                F=list(X =   fixC, df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_ctmp, df = dfig, model = "BRR",saveEffects=TRUE)
                )

   oMpred <- list( M = list(X = SVDM_ctmp, df = dfig, model = "BRR",saveEffects=TRUE))

   Mpred <-  list(
                  F = list(X = fixC, df = dfig, model = "FIXED",saveEffects=TRUE),             
                  M = list(X = SVDM_ctmp, df = dfig, model = "BRR",saveEffects=TRUE))

   GMpred <- list( F = list(X = fixC, df = dfig, model = "FIXED",saveEffects=TRUE),
                  G = list(X = SVDG_ctmp,  df = dfig, model = "BRR",saveEffects=TRUE), 
                  M = list(X = SVDM_ctmp, df = dfig, model = "BRR",saveEffects=TRUE))


   GpredTEST  <- list(
                F=list(X =   fixCT, df = dfig, model = "FIXED",saveEffects=TRUE),
                G = list(X = SVDG_ctest, df = dfig, model = "BRR",saveEffects=TRUE)
                )
   ########################################################################



   ETAsMicro=list("Gpred"=Gpred,"Mpred"=Mpred,"oMpred"=oMpred,"GMpred"=GMpred)
   
   equation1 = function(x) return(rungibbs(yCtmp,ETAsMicro[[x]]))
   equation2 = function(x) return(rungibbs_eq2(myall_p[[x]],Gpred))
   equation2Test = function(x) return(rungibbs_eq2(myall_ptest[[x]],GpredTEST))

   cat("doing first equantion in training set \n")
   EQUATION1=lapply(names(ETAsMicro),equation1)
   names(EQUATION1) = names(ETAsMicro)

   cat("obtained my as test preidctor\n")
   myall_ptest=lapply(2:4,function(i) OTU_cbTEST%*%SVD_tmp$v%*%EQUATION1[[i]]$ETA$M$b)
   names(myall_ptest) = names(ETAsMicro)[2:4]
   
   cat("obtained my as trainig pseudophenotype\n")
   myall_p=lapply(2:4,function(i) SVDM_ctmp%*%EQUATION1[[i]]$ETA$M$b) # pseufo phenotype
   names(myall_p) = names(ETAsMicro)[2:4]

   cat(" solve equation 2 using my in training \n ")
   EQUATION2=lapply(names(myall_p),equation2)

   cat("making alternative prediction uses my in test\n")    
   EQUATION2test=lapply(names(myall_ptest),equation2Test)
   names(EQUATION2test) = names(myall_ptest)

   # fai uno vuuot per rempiere
 
   cat("residual in test set\n")
   residual = lapply(1:3,function(i) EQUATION2test[[i]][["y"]]-EQUATION2test[[i]][["yHat"]] ) 
   residual <- lapply(1:3, function(x) {
                          df <- as.data.frame(residual[[x]])
                         rownames(df) <- nomitest
                                 return(df)})
   
   names(residual)=names(EQUATION1)[2:4]
   cat("add empty vecotr\n")
   empty_pred=as.data.frame(rep(0,nrow(residual[[1]])))
   row.names(empty_pred)=cb$animGLDS[test]

   residual[["Gpred"]] = empty_pred

   cat("extract i predictor \n")
   cat("direct genetic effect Gpred and GMpred\n")
   #
   gdir=names(EQUATION1)[c(1,4)]
   gy_test=lapply(gdir,function(i) SNP_cbTEST%*%SNPSVD_tmp$v%*%EQUATION1[[i]]$ETA$G$b)
   names(gy_test)=gdir
   for(x in names(EQUATION1)[-c(1,4)]) gy_test[[x]]=empty_pred
   
   cat("direct microb effect Mpred\n")
   my_test =  myall_ptest
   names(my_test)=names(EQUATION1)[2:4]
   my_test[["Gpred"]]=empty_pred
   
   cat("mediated genentic effect\n")
   gm_testv1=lapply(1:3,function(i) SNP_cbTEST%*%SNPSVD_tmp$v%*%EQUATION2[[i]]$ETA$G$b) 
   for(i in 1:3) row.names(gm_testv1[[i]])=row.names(yCtest)
   gm_testv1[["Gpred"]]=empty_pred
   
   names(gm_testv1)=names(EQUATION1)[2:4]
   my_test[["Gpred"]]=empty_pred   
   #
   gm_testv2=lapply(1:3,function(i) SVDG_ctest%*%EQUATION2test[[i]]$ETA$G$b) 
   names(gm_testv2)=names(EQUATION1)[2:4]
   gm_testv2[["Gpred"]]=empty_pred   
   #
   
   cat("combine the predictor\n")

   
   g_allv1=lapply(names(EQUATION1),function(x) gy_test[[x]]+gm_testv1[[x]])
   names(g_allv1)=names(EQUATION1)
   
   g_allv2=lapply(names(EQUATION1),function(x) gy_test[[x]]+gm_testv2[[x]])
   names(g_allv2)=names(EQUATION1)
   
   my=lapply(names(EQUATION1),function(x) gy_test[[x]]+my_test[[x]])
   names(my)=names(EQUATION1)

   gv2_res=lapply(names(EQUATION1),function(x) residual[[x]]+g_allv2[[x]])
   names(gv2_res)=names(EQUATION1)
   
   return(list("gdir"=gy_test,"gv1"=g_allv1,"gv2"=g_allv2,"my"=my, "gv2_res"=gv2_res,"test"=test))


}




 get_results=function(){  

   df_list=list("Mpred"=test2,"GMpred"=test1,"oMpred"=test3,"Gpred"=gy_Gpred,
                "mGMpred"=myp1,"mGMpred2"=myp1+gdir,"mMpred"=myp2,"moMpred"=myp3,"y"=yC)

   base=df_list[[1]]   
   colnames(base)=names(df_list)[1]
   for(i in 2:length(df_list))  {
       id=row.names(base)
       base=merge(base,df_list[[i]],by="row.names")[,-1]
       row.names(base)=id       
       colnames(base)[i]=names(df_list)[i]
   }
 
   #------->  breeding values <--------#
   
getlr=function(comp) {
      x=comp[,1]
      y=comp[,2]
      bias=function(x,y) format((mean(x-y)),nsmall=3)
      slope=function(x,y) format(cov(x, y) / var(x) ,nsmall=3)
      acc=function(x,y) format(cor(x,y),nsmall=3)
      return(list(acc=acc(x,y),slope=slope(x,y),bias=bias(x,y)))
}  

# fomuila che moltipocat controlla cosa dice
   lrx=list( oMpred=getlr(merge(test3,ebvoMpred,by="row.names")[,-1]),
         Mpred= getlr(merge(test2,ebvMpred,by="row.names")[,-1]),
         GMpred=getlr(merge(test1,ebvGMpred,by="row.names")[,-1]),
#    getlr(merge(rnd,gdir,by="row.names")[,-1])
         Gpred=getlr(merge(gy_Gpred,ebvGpred,by="row.names")[,-1]))



# predizione se avessi solo il microb
onlym=c(cor(merge(ebvMpred,ebvGpred,by="row.names")[,-1])[2,1],
   cor(merge(ebvoMpred,ebvGpred,by="row.names")[,-1])[1,2])

names(onlym)=c("Mpred","oMpred")
# packed in just one varibake

rtnr=list(ypred=base,indi=onlym,lr=lrx)

return(rtnr)
}

  
  






   