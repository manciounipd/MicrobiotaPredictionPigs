setwd("/home/enrico/Work/prediction_micro_v2/")

system("ls *")

# Clear environment
rm(list = ls())

merge_by_row=function(db1,db2){
  db=merge(db1, db2, by = "row.names", all = F)
  row.names(db)=db[,1]
  db=db[,-1]
  return(db)
}


getlr=function(comp) {
      x=comp[,1]
      y=comp[,2]
      bias=function(x,y) format((mean(x-y)),nsmall=3)
      slope=function(x,y) format(cov(x, y) / var(x) ,nsmall=3)
      acc=function(x,y) format(cor(x,y),nsmall=3)
      return(list(acc=acc(x,y),slope=slope(x,y),bias=bias(x,y)))
}  


getphen=function(comp) {

      x=comp[,1]
      y=comp[,2]
      
      rmse_na <- function(true, pred) {
        idx <- complete.cases(true, pred)
        sqrt(mean((true[idx] - pred[idx])^2))
      }

      cor_na <- function(true, pred) {
        idx <- complete.cases(true, pred)
        cor(true[idx], pred[idx])
      }

      return(list(acc=cor_na(x,y),rmse=rmse_na(x,y)))
}  


collectlr <- function() {
  ITER <- 1
  collect <- list()
  predicto <- c("gv2", "gv1", "gdir")

  for (nCV in seq_along(results)) {
    for (PREDI in predicto) {
      tmp <- results[[nCV]][[PREDI]]
      nomi <- names(tmp)

      for (xi in nomi) {
        whole <- wg[[xi]]
        row.names(whole) <- row.names(yC)
        partial <- tmp[[xi]]
        merg <- data.frame(merge_by_row(partial, whole))
        collect[[ITER]] <- cbind(
          it = nCV,
          nomi = xi,
          do.call("cbind", getlr(merg)),
          pheno = pheno,
          predict = PREDI
        )
        ITER <- ITER + 1
      }
    }
  }

  require(dplyr)
  require(reshape2)

  tmp <- as.data.frame(do.call("rbind", collect)) %>%
    reshape2::melt(id = c("it", "nomi", "pheno", "predict")) %>%
    mutate(value = as.double(value))

  return(tmp)
}


collecpredphen <- function() {
  predicto <- c("gv2", "gv1", "gdir", "my", "gv2_res")
  collect <- list()
  ITER <- 1

  # Loop through each cross-validation fold and predictor
  for (nCV in seq_along(results)) {
    for (PREDI in predicto) {
      tmp <- results[[nCV]][[PREDI]]
      nomi <- names(tmp)

      for (xi in nomi) {
        partial <- tmp[[xi]]
        merg <- data.frame(merge_by_row(partial, yC))
        collect[[ITER]] <- cbind(
          it = nCV,
          nomi = xi,
          do.call("cbind", getphen(merg)),
          pheno = pheno,
          predict = PREDI
        )
        ITER <- ITER + 1
      }
    }
  }

  require(dplyr)
  require(reshape2)

  tmp <- as.data.frame(do.call("rbind", collect)) %>%
    reshape2::melt(id = c("it", "nomi", "pheno", "predict")) %>%
    mutate(value = as.double(value))

  return(tmp)
}




## faccio la stessa cosa ma per i fenomitpo
setwd("/home/enrico/Work/prediction_micro_v2/")
require("tidyverse")

#system("ls CV_withincage/*")

lrall <- data.frame()
yphen <- data.frame()

#system("ls CV_withincage")

# Set date for file naming
current_date <- "2025-11-19"
# Phenotypes to loop over
phenotypes <- c("CIG_MarketBF", "CIG_CarcDGWean_kgs",
                  "CIG_MarketLEA",
                "CIG_LiveAverageDailyGain_W2M_kgs")

system("ls *")
phenotypes=system(paste0("ls CV_withincage/pred_",current_date,"_*"),intern=TRUE)

for (pheno in phenotypes) {
  
  cat("@@@@@@@@@@@@@@@@@@@@@@ Loading analysis =>", pheno, "..")
  file_path <- pheno #paste0("CV_withincage/pred_", current_date, "_", pheno, ".Rdata")
  load(file_path) # ricordati che ci sono anche le funizone
  cat("read ! @@@@@@@@@@@@@@@@@@@@@@@@@@\n")
 # lrall=rbind(lrall,collectlr())
  yphen=rbind(yphen,collecpredphen())
  
}


# Sistema i files 

cat("finish..................\n")


cat("A) organize the results")
library(ggplot2)
require("tidyverse")

ab=yphen %>% mutate(value=as.double(value)) %>% 
      group_by(nomi,predict,pheno,variable) %>% filter(!is.na(value) ) %>%
  summarize(
    mean = mean(value),
    sd = sd(value),
    q025 = quantile(value, 0.025),
    q50  = quantile(value, 0.5),
    q975 = quantile(value, 0.975),
    .groups = "drop"
  )   


ab$nomi <- factor(ab$nomi, levels = rev(levels(factor(ab$nomi))))
ab$pheno=as.factor(ab$pheno)
ab$pheno2=ab$pheno
table(ab$variable)
levels(ab$pheno2)=c("IMF","MARB","CARCW","ADG","BF","LD","LEA")

save(file="rdata/within.RData",x=ab)


# Prepare dataset for plotting

rt=ab %>%  filter(!nomi %in% c("Gpred") & !predict %in% c("gdir"))

tmp <- ab %>%
  filter(nomi == "Gpred") %>%
  select(pheno2,variable, mean,q025,q975) %>%
  filter(!is.na(mean)) %>%
  distinct()


getwd()
dodge <- position_dodge(width = .7, preserve = "single")
dodge2 <- position_dodge(width = .7, preserve = "single")

dir.create("plot")

rt$predict=as.factor(rt$predict)
levels(rt$predict)=c("S2","S3","S4","S1")
rt$predict=as.factor(as.character(rt$predict))

rt=rt%>%filter(pheno2 %in% c("ADG","BF","IMF"))
tmp=tmp%>%filter(pheno2 %in% c("ADG","BF","IMF"))

tmp$pheno2 <- factor(tmp$pheno2, levels = c("ADG","BF","IMF"))
rt$pheno2 <- factor(rt$pheno2, levels = c("ADG","BF","IMF"))



plotta <- ggplot(rt, aes(fill = predict, y = mean, x = nomi, group = predict)) +
  geom_col(position = dodge, width = 0.6, color = "black", alpha = 0.6) +
  geom_errorbar(
    aes(ymin = q025, ymax = q975),
    position = dodge,
    width = 0.2,
    color = "black"
  ) +
  geom_hline(data = tmp, aes(yintercept = mean) ,color = "red") +
  geom_hline(data = tmp, aes(yintercept = q975), linetype = "dashed", color = "red") +  
  geom_hline(data = tmp, aes(yintercept =   q025), linetype = "dashed", color = "red") +
  facet_grid(pheno2 ~ variable) +
  geom_text(
    aes(y = q975, label = sprintf("%.2f", mean)),
    position = position_dodge2(width = .7),
    vjust = -0.3
  ) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme_bw() 


png("plot/WITHIN.png", units="in", width=8.5, height=12, res=300)
plotta+xlab("")+ylab(expression(paste("Cor(", hat(y[v]), ", ", y[v], ")")))+
labs(fill="Model")
dev.off()


dir.create("res")
tmp$nomi="Gpred"
W=plyr::rbind.fill(rt,tmp)

save(file="res/cvW.RData",x=W)
