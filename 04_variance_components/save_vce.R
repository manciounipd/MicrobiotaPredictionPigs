print("load vce")
source("functions/function_plot.R")
setwd("/home/enrico/Work/prediction_micro_v2/")
load("rdata/all_vceCB.RData")

traits <- names(collect1) 
components <- c("Gpred", "GMpred", "Mpred", "oMpred")

# --- helper functions ---

get_r <- function(ok) {
  tmp <- apply(ok, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975)))
  fx  <- function(x) format(x, digits=3, nsmall=3)
  m1  <- apply(tmp, 2, function(x) paste(fx(x[2]), "(", fx(x[1]), fx(x[3]), ")"))
  return(m1)
}

# clean std_to1
std_to1 <- function(x) {
  x <- as.data.frame(x)
  row_sums <- rowSums(x)
  row_sums[row_sums == 0] <- NA
  sweep(x, 1, row_sums, "/")
}

get_perc <- function(ok) {
  ok  <- std_to1(ok)
  tmp <- apply(ok, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
  fx  <- function(x) format(x, digits=3, nsmall=3)
  m1  <- apply(tmp, 2, function(x) paste0(fx(x[2]), " (", fx(x[1]), ", ", fx(x[3]), ")"))
  return(m1)
}

get_perccorr <- function(ok, v) {
  ok  <- std_to1(ok)
  tmp <- apply(ok, 2, function(x) quantile(x, probs=c(0.025, 0.5, 0.975), na.rm=TRUE))
  fx  <- function(x) format(x, digits=3, nsmall=3)
  m1  <- apply(tmp, 2, function(x) paste0(fx(x[2]*v), " (", fx(x[1]*v), ", ", fx(x[3]*v), ")"))
  return(m1)
}

# --- collect results ---

eq1=list(); eq2=list(); eq2_per=list(); eq22=list()
x=1
for (trait in traits) {
  for (comp in components) {

    m1 = t(as.data.frame(get_perc(collect1[[trait]][[comp]])))
    meM = mean(collect1[[trait]][[comp]][["M"]])
    eq1[[x]] = data.frame(Model=comp, trait=trait, m1)

    mi = collect2[[trait]][[comp]]
    if (!is.null(mi)) {
      m2   = t(as.data.frame(get_perc(mi)))
      m2a  = t(data.frame(get_r(mi)))
      m32a = t(data.frame(get_perccorr(mi, meM)))

      eq2_per[[x]] = data.frame(Model=comp, trait=trait, as.data.frame(m2))
      eq2[[x]]     = data.frame(Model=comp, trait=trait, m2a)
      eq22[[x]]    = data.frame(Model=comp, trait=trait, m32a)
    }
    x = x + 1
  }
}

# --- build databases ---

db   <- as.data.frame(dplyr::bind_rows(eq1))
db2  <- as.data.frame(dplyr::bind_rows(eq2_per))
db2w <- as.data.frame(dplyr::bind_rows(eq2))
db22 <- as.data.frame(dplyr::bind_rows(eq22))

library(ggplot2)
if (requireNamespace("webr", quietly = TRUE)) library(webr)
library(dplyr)
library(tidyverse)


#---------------------------------------------#
# --- functions for processing/plotting ---
#------------------------------------------#

GET_DB_Gpred=function(PHENO){
    TMP=db[db$trait==PHENO & db$Model=="Gpred",-1,]
    tmp=data.frame();for(i in names(TMP)[-1]) tmp=rbind(tmp,cbind(i,TMP[,i] ))    
    tmp=tmp[!is.na(tmp$V2),]
    tmp$val = as.double(do.call("c", lapply(1:nrow(tmp), function(cx) stringr::str_split_1(tmp$V2[cx], " ")[1])))
    tmp$subeffect=tmp$i
    tmp$effect=""
    return(tmp)
}

GET_DB <- function(PHENO, MODEL) {

  p2 = as.data.frame(t(db22[db22$trait==PHENO & db22$Model==MODEL, -c(1:2)]))
  p1 = t(db[db$trait==PHENO & db$Model==MODEL, -c(1:2)])
  

      pp = as.data.frame(rbind(
        cbind(type="nested", subeffect="M", effect=row.names(p2), v=p2[,1]),
        cbind(type="whole",  subeffect=row.names(p1), effect=row.names(p1), v=p1[,1])
        ))

      pp = pp[!is.na(pp$v),]

      tmp = pp[pp$type=="whole" & pp$effect!="M",]
      tmp$type = "nested"
      pp = rbind(tmp, pp)
     
     

      pp$val = as.double(do.call("c", lapply(1:nrow(pp), function(cx) stringr::str_split_1(pp$v[cx], " ")[1])))
    
      PP = pp[pp$effect!="M",]
      PPx = PP %>% filter(type=="nested")
      PPx[PPx$subeffect!="M","effect"] = ""

      df_agg = aggregate(val ~ subeffect + effect, data=PPx, sum)
      df_agg$effect2 = paste0(df_agg$val, "eff")
      df_agg[df_agg$subeffect!="M","effect2"] = ""
 
      return(df_agg)
}


ALL_EFFECTS <- c("m[y]","e[y]", "g[y]", "sys[y]", "g[m]", "res[m]", "sys[m]")


GET_DB_X_PLOT <- function(df) {

  inner <- df %>%
    dplyr::group_by(subeffect) %>%
    dplyr::summarise(val = sum(val), .groups = "drop") %>%
    dplyr::mutate(level = 1, label = subeffect)

  outer <- df %>%
    dplyr::filter(subeffect == "M") %>%
    dplyr::mutate(level = 2, label = effect)

    inner <- inner %>%
        dplyr::mutate(label2 = dplyr::case_when(
          
          label == "F"   ~ "sys[y]",
          label == "RES" ~ "e[y]",
          label == "G"   ~ "g[y]",
          label == "M"   ~ "m[y]",
          TRUE ~ label
        ))
  
  outer <- outer %>%
        dplyr::mutate(label2 = dplyr::case_when(
          effect == "G"  ~ "g[m]",
          effect == "F"  ~ "sys[m]",
          effect == "RES"~ "res[m]",
          TRUE ~ label
        ))

  donut_df <- dplyr::bind_rows(inner, outer)

  # ensure all effects are present
  donut_df <- donut_df %>%
    tidyr::complete(label2 = ALL_EFFECTS, fill = list(val = 0, level = 1))

  # fix factor order (m[y] last)
  donut_df$label2 <- factor(donut_df$label2, levels = ALL_EFFECTS)
#  donut_df$label <- factor(donut_df$label, levels=c("F","G","RES","M"))

  # recompute positions
  donut_df <- donut_df %>%
    dplyr::group_by(level) %>%
    dplyr::arrange(level, label2) %>%
    dplyr::mutate(
      ymax = cumsum(val),
      ymin = ymax - val,
      label_pos = (ymax + ymin) / 2,
      label3 = paste0(label2, " ~ '", round(val*100, 0), "%'")
    )

  return(donut_df)
}



# --- example usage ---
PHENOS=traits; names(PHENOS)=c("BF","IMF","ADG")
PHENOS=traits; names(PHENOS)=c("BF","IMF","ADG")

PLOY = list()
LEGEND <- NULL



new_name=c("M1","NO","M12","M13")
old_name=c("Gpred","oMpred","Mpred","GMpred")
t=data.frame(old_name,new_name)

for(PH in c(3,1,2)) {

  PHENO = PHENOS[PH]
  donut_df2 <- data.frame()
  L <- list()
   
  for(i in c("Gpred","Mpred","GMpred")){

    if(i=="Gpred") df_agg  = GET_DB_Gpred(PHENO) else df_agg   <- GET_DB(PHENO, i)
    xx=t[t$old_name==i,]$new_name
    donut_df =  GET_DB_X_PLOT(df_agg)
    donut_df$MODELLO <-xx 
    L[[i]] <- MAKE_PLOT(donut_df, xx)+theme(plot.margin = unit(c(0,0,0,0), "cm"))#+theme(legend.position = "none")
    donut_df2 <- dplyr::bind_rows(donut_df2, donut_df)
  }
    TMP = ggpubr::ggarrange(plotlist = c(L), ncol = 3, nrow = 1,legend="none", common.legend = TRUE)

 # if(names(PHENO)=="IMF") TMP = ggpubr::ggarrange(plotlist = c(L), ncol = 3, nrow = 1,legend="bottom", common.legend = TRUE)
 # if(!names(PHENO)=="IMF") TMP = ggpubr::ggarrange(plotlist = c(L), ncol = 3, nrow = 1,legend="none", common.legend = TRUE)
 
  # ggpubr::annotate_figure() uses gridExtra::arrangeGrob(padding = unit(0.5, "line"))
  # which adds visible whitespace; keep padding at 0 to make the layout tighter.
  TMP = gridExtra::arrangeGrob(
    TMP,
    left = ggpubr::text_grob(names(PHENO), face = "italic", size = 25),
    padding = grid::unit(0, "pt")
  ) %>%
    ggpubr::as_ggplot()
  PLOY[[PHENO]] = TMP

  TMP =L[[max(length(L))]]+theme(legend.position = "bottom") +
    guides(fill = guide_legend(order = 1),
           color = guide_legend(order = 2))
}


library(patchwork)

plot3 <- PLOY[["CIG_LiveAverageDailyGain_W2M_kgs"]]+PLOY[["CIG_MarketBF"]]
 

final_plot=ggpubr::ggarrange(plotlist=PLOY,ncol=1,hjust = -0.1, vjust = 0.5,legend="none", common.legend = TRUE,align='h')

for(i in names(PLOY)){
png(paste0("plot/Donut_Plots",i,".png"), units = "in", width = 15, height = 17/3, res = 300)
print(PLOY[[i]])
dev.off()}

# save to file
png("plot/Donut_PlotsPB.png", units = "in", width = 15, height = 17, res = 300)
print(final_plot)
dev.off()


png("plot/Donut_Plots_Legend.png", units = "in", width = 5, height = 6, res = 300)
print(TMP + theme(
  legend.position = "right",   # or "left"
  legend.direction = "vertical"
))
dev.off()


LEGGEND = L[[3]]+
  lims(x = c(0,0), y = c(0,0))+
  theme_void()+
  theme(
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.position = c(0.5,0.5),
        legend.key.size = unit(1, "cm"),
        legend.text = element_text(size =  12),
        legend.title = element_text(size = 15, face = "bold"))+
  guides(colour = guide_legend(override.aes = list(size=8)))
# Create plot

png("plot/Donut_PlotsCB.png", units="in", width=7*2.8, height=7*2.5, res=300)
PLOY_tight <- lapply(PLOY, function(p) p + theme(plot.margin = margin(0, 0, 0, 0, unit = "pt")))

legend_grob <- ggpubr::get_legend(LEGGEND + theme(legend.position = "right", legend.direction = "vertical"))

plot_grid <- cowplot::plot_grid(
  plotlist = PLOY_tight,
  ncol = 1,
  nrow = 3,
  align = "v"
)

# ggpubr combines legend with gridExtra::arrangeGrob(padding = unit(0.5, "line")).
# Build it manually with padding = 0 for a tighter result.
lwidth <- sum(legend_grob$width)
plot_all <- ggpubr::as_ggplot(
  gridExtra::arrangeGrob(
    plot_grid,
    legend_grob,
    ncol = 2,
    widths = grid::unit.c(grid::unit(1, "npc") - lwidth, lwidth),
    padding = grid::unit(0, "pt")
  )
) + theme(plot.margin = margin(0, 0, 0, 0, unit = "pt"))

print(plot_all)


dev.off()


getwd()
db=data.frame()
for(i in 1:length(eq2))    db=rbind.fill(db,eq2_per[[i]])

dir.create("table")
data.table::fwrite(file="table/tab2.csv",db,sep=",")


# Combine all rows into a single data frame
summary_table <- do.call(rbind, results)
row.names(summary_table) <- NULL

# View the final result
print(summary_table)
