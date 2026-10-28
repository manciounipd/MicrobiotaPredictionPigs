setwd("/home/enrico/Work/prediction_micro_v2/")

system("ls CV_multi/*")

# Clear environment
rm(list = ls())


# Set date for file naming

current_date <- "2025-11-18"
# Phenotypes to loop over

# Initialize result containers
lr1 <- list()
lr2 <- list()
lr3 <- list()
lr4 <- list()

phenotypes=system(paste0("ls CV_multi/pred_",current_date,"*"),intern=TRUE)

geneticvalue= FALSE

# Main loop over phenotypes
if (geneticvalue) {
  for (pheno in phenotypes) {
    cat("@@@@@@@@@@@@@@@@@@@@@@\nProcessing =>", pheno, "\n")

    load(pheno)  # Make sure this loads the required variables

    for (i in 1:4) row.names(gm_pb2[[i]]) <- row.names(yP)

    partial_g2 <- lapply(1:4, function(i) gm_pb2[[i]] + gy_pb[[i]])
    residuals <- lapply(1:4, function(i) Partial_eq2V2[[i]][["y"]] - Partial_eq2V2[[i]][["yHat"]])
    partial_predv2 <- lapply(1:4, function(i) residuals[[i]] + partial_g2[[i]])


      whole_PB <- lapply(1:4, function(i) data.frame(Whole = Whole_g[[i]][, 2]))

      tmp1 <- do.call("rbind", lapply(1:4, function(i) {
        paste(names(ETA_W)[i], getlr(Whole_g[[i]]))
      }))
      tmp1 <- as.data.frame(tmp1); names(tmp1) <- c("acc", "slope", "bias")

      cvlr <- function(i) getlr(merge_by_row(partial[[i]], whole[[i]]))

      create_summary_table <- function() {
        summary <- data.frame(do.call("rbind", lapply(1:4, cvlr)))
        row.names(summary) <- names(ETA_W)
        summary$phen <- pheno
        return(summary)
      }

      # Set up and compute summaries
      whole <- whole_PB; partial <- partial_g2; lr1[[pheno]] <- create_summary_table()
      whole <- Pb_g;      partial <- whole_PB;  lr2[[pheno]] <- create_summary_table()
      whole <- whole_PB;  partial <- partial_g; lr3[[pheno]] <- create_summary_table()
      whole <- Pb_g;      partial <- partial_g; lr4[[pheno]] <- create_summary_table()

      # Final table
      a1 <- cbind(do.call("rbind", lr1), predtype = "pred2")
      b1 <- cbind(do.call("rbind", lr3), predtype = "pred1")
      print(as.data.frame(rbind(a1, b1)))
  }



    cat("## Final Phenotype Prediction Evaluation ##\n")
    a1=cbind(do.call("rbind",lr1),predtype="pred2")
    b1=cbind(do.call("rbind",lr3),predtype="pred1")

    as.data.frame(rbind(a1,b1))  
}
    #Now create the table with the results or plot


cat("## Final Phenotype Prediction Evaluation ##\n")

# Store RMSE and correlations
pheno_predictions <- list()

calculate_metrics <- function(df, col_index) {
    sqrt(sum((df[, 1] - df[, col_index])^2, na.rm = TRUE) / nrow(df))
  }




for (pheno in phenotypes) {
  cat("Processing =>", pheno, "\n")
  load(pheno)  # Make sure this loads the required variables
  

  results <- list()
  for (i in 1:4) {

    predictions <- merge_by_row_l(list(yP, partial_g2[[i]], partial_g[[i]], partial_my[[i]], partial_predv2[[i]]))
    names(predictions) <- c("yP", "Partial_g2", "Partial_g", "Partial_my", "Partial_predv2")
    round=1
    rmse=list();corr=list()
    repeat {
      x=sample(1:nrow(predictions),size=round(nrow(predictions)*0.8))
      rmse[[round]]=do.call("c", lapply(2:5, function(j) calculate_metrics(predictions[x,], j)))
      corr[[round]]=cor(predictions[x,], use = "pairwise.complete.obs")[1, -1]
      round=round+1
      if(round>1000) {break}
    }
    
    rmse=t(apply(do.call("rbind",rmse),2, function(x) quantile(x, na.rm=T, probs = c(0.025,0.5,0.975))))
    corr=t(apply(do.call("rbind",corr),2, function(x) quantile(x, na.rm=T, probs = c(0.025,0.5,0.975))))
    row.names(rmse)=row.names(corr)
    metrics=rbind(data.frame(corr,"metrics"="corr"),data.frame(rmse,"metrics"="rmse"))

    metrics$model <- names(ETA_W)[i]
    metrics$phen <- pheno
    metrics$type <- rep(row.names(corr),2)
    names(metrics)[1:3] = c("lower","median","upper")
    results[[i]] <- metrics
  }

  pheno_predictions[[pheno]] <- do.call("rbind", results)
}


final_results=do.call("rbind",pheno_predictions)
row.names(final_results)=1:nrow(final_results)


require("tidyverse")

# Extract Gpred correlations per phenotype

print("Stesiiiiiiiiiiii fenotipi!!!!!!!!!!!!!!!!!!!!!!!!!")

final_results$phen=as.factor(final_results$phen)
levels(final_results$phen) = c("ADG","BF","IMF","LEA")
final_results = final_results %>% filter(phen %in% c("ADG","BF","IMF"))

final_results$type=as.factor(final_results$type)
levels(final_results$type) = c("S2","S3","S1","S4")
final_results$type=as.character(final_results$type)
final_results$type=as.factor(final_results$type)

row.names(final_results)=1:nrow(final_results)



save(file="res/cv_Between.RData",x=final_results)



tmp <- final_results %>%
  filter(model == "Gpred") %>%
  select(phen, cor) %>%
  filter(!is.na(cor)) %>%
  distinct()


final_results$model=as.factor(final_results$model)
final_results$model <- factor(final_results$model, levels = c("Gpred","oMpred","Mpred","GMpred"))
# Make the plot
dodge <- position_dodge(width = .7, preserve = "single")

plt=final_results %>%
  filter(model != "Gpred") %>%
  ggplot(aes(fill = type, y = cor, x = model, label = round(cor, 3))) +
  geom_col(position = dodge, width = 0.6, color = "black", alpha = 0.6) +
  facet_grid(phen ~ .) +
  geom_text(
    aes(y = cor, label = sprintf("%.2f", cor)),
    position = position_dodge2(width = .7),
    vjust = -0.3
  )+
  geom_hline(data = tmp, aes(yintercept = cor), linetype = "dashed", color = "red")+
    geom_hline(data = tmp, aes(yintercept = cor), linetype = "dashed", color = "red")+
      geom_hline(data = tmp, aes(yintercept = cor), linetype = "dashed", color = "red")+
   theme_bw() +
  theme(
    axis.text.y = element_text(size = 15, face = "italic", hjust = 1),
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


png("plot/ACROSS.png", units="in", width=8.5, height=12, res=300)
plt+xlab("")+ylab(expression(paste("Cor(", hat(y[v]), ", ", y[v], ")")))+
labs(fill="Model")
dev.off()




final_results %>%
  filter(!model %in% c("Gpred","Mpred") & !type %in% c("Partial_predv2")) %>%
  ggplot(aes(x = type, y = cor, color="black",fill = model, label = round(cor, 3))) +
  geom_col(position = position_dodge2()) +
  facet_grid(phen ~ .) +
  geom_text(position = position_dodge2(width = 0.9), vjust = -0.3) +
  theme_bw() +
  geom_hline(data = tmp, aes(yintercept = cor), linetype = "dashed", color = "red")




tmp <- final_results %>%
  filter(model == "Gpred") %>%
  select(phen, rmse) %>%
  filter(!is.na(rmse)) %>%
  distinct()


final_results %>%
  filter(!model %in% c("Gpred","Mpred") & !type %in% c("Partial_predv2")) %>%
  ggplot(aes(x = type, y = rmse,group=model, color=model,fill = model, label = round(rmse, 3))) +
   geom_point(aes(shape=model),size=4,na.rm=TRUE, position=position_dodge(width=0.3))+
  facet_grid(phen ~ .,scale="free") +
  geom_text(position = position_dodge2(width = 0.3), vjust = -1) +
  theme_bw() 
  geom_hline(data = tmp, aes(yintercept = rmse), linetype = "dashed", color = "red")
