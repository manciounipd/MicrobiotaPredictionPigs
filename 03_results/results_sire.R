# divio anche i pb con i cb ?

main_fx=function(result) {

models=c("Gpred","Mpred","oMpred","GMpred")
# Initialize an empty data frame to store correlation and RMSE results
metrics_results <- data.frame(nomi = character(), nIter=numeric(),model = character(), correlation = numeric(), rmse = numeric(), stringsAsFactors = FALSE)
for(nIter in  1:length(result)){
# Loop through each element of result
for (nomi in names(result[[nIter]])) {
    for (model in models) {

        if (!is.null(nrow(result[[nIter]][[nomi]][[model]])) ) {
            tmp_tmp <- as.data.frame(result[[nIter]][[nomi]][[model]])
            names(tmp_tmp) <- c("pb", "cb")
            tmp_tmp$RowNames <- row.names(tmp_tmp)
            tmp_tmp <- reshape2::melt(tmp_tmp, id.vars = "RowNames")
            
            Y <- as.data.frame(Y)
            names(Y) <- c("pb", "cb")
            Y$RowNames <- row.names(Y)
            Y_tmp <- reshape2::melt(Y, id.vars = "RowNames")
            
            comp_tmp <- merge(tmp_tmp, Y_tmp, by = c("RowNames", "variable"))
            comp_tmp <- comp_tmp[complete.cases(comp_tmp), ]
            comp_tmp[,"value.x"] <- as.double(comp_tmp[,"value.x"])
            comp_tmp[,"value.y"] <- as.double(comp_tmp[,"value.y"])
            
            # Calculate correlation
            corr <- cor(comp_tmp[,"value.x"], comp_tmp[,"value.y"], use = "pairwise.complete.obs")
            
            # Calculate RMSE
            rmse <- sqrt(mean((comp_tmp[,"value.x"] - comp_tmp[,"value.y"])^2))
            
            # Append the metrics result to the data frame
            metrics_results <- rbind(metrics_results, data.frame(nomi = nomi,nIter=nIter, model = model, correlation = corr, rmse = rmse))
            
            cat("Metrics for nomi:", nomi, "model:", model, "Correlation:", corr, "RMSE:", rmse, "\n")
        }
    }
}

}
    return(metrics_results)
}

# Print the metrics results data fram

require("tidyverse")
# 
setwd("/home/enrico/Work/prediction_micro_v2/")
env_folder = "/home/enrico/Work/prediction_micro_v2/CV_sire/"
data="2025-11-19"

pheno = c("CIG_MarketBF", #"CIG_CarcDGWean_kgs",
          "ActIMF",
          "CIG_LiveAverageDailyGain_W2M_kgs")

env_files <- list.files(env_folder , full.names = TRUE)

## keep only files matching phenotypes and, if multiple per phenotype,
## retain the most recent (by date encoded in the filename)
env_files <- env_files[Reduce(`|`, lapply(pheno, grepl, x = env_files))]

selected_files <- character()
for (p in pheno) {
  matches <- env_files[grepl(p, env_files)]
  if (length(matches) > 0) {
    matches <- sort(matches)
    selected_files <- c(selected_files, tail(matches, 1))
  }
}
env_files <- unique(selected_files)


env_files <- c(env_files,list.files(env_folder, pattern = data, full.names = TRUE))

all_results <- list()
# Loop through each file and load the environment

for (file in env_files) {
    cat("Loading:", file, "\n")
    load(file)
    tmp=main_fx(result)
    add_nome = sub("pred_[0-9-]+_", "", tools::file_path_sans_ext(basename(file)))
    tmp$phen=add_nome
    all_results[[file]] <- tmp
}

metrics_results = do.call(rbind, all_results)
row.names(metrics_results) = 1:nrow(metrics_results)

fao=metrics_results %>% reshape2::melt(id=(c("nomi","nIter","model","phen")))
table(fao$phen)

summario=fao %>% group_by(nomi,phen,model,variable) %>% 
     filter(!is.na(value)) %>%
  summarize(mean = mean(value),
            q025 = quantile(value, 0.025),
            q975 = quantile(value, 0.975),
            sd = sd(value)) %>% ungroup()


table(summario$phen)

#            mean_rmse = mean(rmse), 
#            sd_rmse = sd(rmse)) 

# Extract Gpred correlations per phenotype
tmp = summario %>%
  filter(model == "Gpred") %>%
  select(phen,variable,mean,q025,q975) %>% ungroup() %>%
  dplyr::distinct()

summario=summario %>% filter(!model=="Gpred" & ! nomi=="gydir")


summario$model=as.factor(summario$model)
summario$model <- factor(summario$model, levels = c("oMpred","Mpred","GMpred"))

summario$nomi=as.factor(summario$nomi)
summario$nomi= factor(summario$nomi, levels = c("my","gmv1","gmv2","gmres2","gydir"))
levels(summario$nomi)=paste0("S",1:5)

dodge <- position_dodge(width = .7, preserve = "single")
dodge2 <- position_dodge(width = .7, preserve = "single")

table(summario$phen)

summario=summario %>% filter(phen  %in% c("CIC_ActIMF","CIG_MarketBF","CIG_LiveAverageDailyGain_W2M_kgs"))
summario$phen <- as.factor(summario$phen)# levels = c("CIG_MarketBF", "CIG_LiveAverageDailyGain_W2M_kgs", "pred_2025-06-15_CIC_ActIMF"))
levels(summario$phen) <- c("IMF","ADG","BF")

tmp=tmp %>% filter(phen  %in% c("CIC_ActIMF","CIG_MarketBF","CIG_LiveAverageDailyGain_W2M_kgs"))
tmp$phen <- as.factor(tmp$phen)# levels = c("CIG_MarketBF", "CIG_LiveAverageDailyGain_W2M_kgs", "pred_2025-06-15_CIC_ActIMF"))
levels(tmp$phen) <-  c("IMF","ADG","BF")
tmp$model="Gpred"

tail(summario)

a=summario%>% 
 ggplot(aes(fill = nomi, y = mean, x = model, label = round(mean, 3))) +
        geom_col(position = dodge, width = 0.6, color = "black", alpha = 0.6) +
  geom_errorbar(
    aes(ymin = q025, ymax = q975),
    position = dodge,
    width = 0.2,
    color = "black"
  ) +
   geom_text(
    aes(y = q975, label = sprintf("%.2f", mean)),
    position = position_dodge2(width = .7),
    vjust = -0.3
  ) +
  geom_hline(data = tmp, aes(yintercept = mean), linetype = "dashed", color = "red") +
  facet_grid(phen ~ variable) +
  #scale_fill_manual(values = c("#00AFBB", "#E7B800")) +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 10, face = "italic", hjust = 1),
    axis.text.x = element_text(size = 14, face = "italic"),
    strip.text.y = element_text(size = 10, face = "italic", angle = .9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    legend.title = element_text(size = 13, face = "italic"),  # Legend title
    legend.text = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "italic")
  )


png("plot/Sire.png", units="in", width=8.5, height=12, res=300)
a
dev.off()

dir.create("res")
S=plyr::rbind.fill(tmp,summario)
save(file="res/cvSire.RData",x=S)

########## giiungi i resiudiiiiiii########
