#

rm(list=ls())
require("tidyverse")

setwd("/home/enrico/Work/prediction_micro_v2/")

for (f in list.files("res", pattern = "\\.RData$", full.names = TRUE)) {
  obj_name <- gsub("\\.RData$", "", basename(f))
  tmp <- load(f)            # load restituisce il nome dell'oggetto caricato
  assign(obj_name, get(tmp)) # prendi l'oggetto e assegnalo con il nuovo nome
}


ls()

row.names(cv_Between)=1:nrow(cv_Between)

W=W[,c("predict","variable","nomi","pheno2","q50","q025","q975")]

S=S[,c("nomi","variable","model","phen","mean","q025","q975")]

B=cv_Between[!is.na(cv_Between$median) ,c("type","metrics","model","phen","median","lower","upper")]


nomi=c("scenario","metrics","model","phenotype","median","lower","upper")



names(W)=nomi
names(B)=nomi
names(S)=nomi

S[is.na(S$scenario) & S$metrics=="correlation",]

W$CV="Within"
B$CV="Between"
S$CV="Sire"



names(W)
names(S)
names(B)

table(S$phenotype)
View(S %>% filter(model=="Gpred") )

all=rbind(as.data.frame(W),as.data.frame(S),as.data.frame(B))

tmp=all %>%
  filter(model %in% c("Gpred") & !scenario %in% c("S0")) %>%
  filter(!phenotype %in% "LEA") %>%
  mutate(CV = factor(CV, levels = c("Within", "Between", "Sire")))

tmp$median=apply(tmp[,c("lower","upper")],1,median)


db_=list()
for (i in unique(all$model)[1:3] ) {

gh=all %>%
  filter(model %in% c(i) & !scenario %in% c("S0")) %>%
  filter(!phenotype %in% "LEA") %>%
#  mutate(scenario = factor(scenario, levels = c("S2", "S3", "S4", "S1"))) %>%
  mutate(CV = factor(CV, levels = c("Within", "Between", "Sire")))
 gh$model=i
 db_[[i]]=gh
}
  


db=all
Gpred=db[db$model=="Gpred",c("model","metrics","phenotype","median","lower","upper","CV")]  %>% distinct()
Gpred[is.na(Gpred$median),c("median")] = apply(Gpred[is.na(Gpred$median),c("lower","upper")],1,median)
db=db[!db$model=="Gpred",]

db$scenario=gsub("S","",db$scenario)
db[db$model=="oMpred","model"]="M2"
db[db$model=="Mpred","model"]="M3"
db[db$model=="GMpred","model"]="M4"

db$d=(as.integer(db$scenario)-1)
db$d[db$d==0]=""

tmp=paste0(db$model,db$d)

x_switched <- sapply(tmp, function(s) {
  paste0(substr(s, 1, 1), substr(s, 3, 3), substr(s, 2, 2))
})

db$model_all=x_switched
db[db$scenario==1,]

table(db$CV)
head(db)

db[!db$metrics=="rmse","metrics"]="acc"
Gpred[!Gpred$metrics=="rmse","metrics"]="acc"

getwd()

db = db[!db$model_all %in% c("M2","M12","M22","M32"),]

db[db$model_all %in% c("M3"),"model_all"]="M2"
db[db$model_all %in% c("M4"),"model_all"]="M3"

db[db$model_all %in% c("M13"),"model_all"]="M12"
db[db$model_all %in% c("M14"),"model_all"]="M13"

db[db$model_all %in% c("M23"),"model_all"]="M22"
db[db$model_all %in% c("M24"),"model_all"]="M23"

db[db$model_all %in% c("M33"),"model_all"]="M32"
db[db$model_all %in% c("M34"),"model_all"]="M33"

col_blue <- "#8FBBD9"   # blu pastello
col_red  <- "#E6A0A0"   # rosso pastello

i="Sire"

for(i in unique(db$CV)) {
print(i)

PLOT=db %>% 
  filter(CV == i & metrics == "acc") %>%
  ggplot(aes(x = model_all, y = median, fill = model, label = round(median, 2))) +
  
  geom_col(position = position_dodge2()) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge2(width = 0.9),
    width = 0.25,
    color = "black"
  ) +
  
  facet_grid(
    phenotype ~ paste0("S", scenario),
    scales = "free_x"
  ) +
  
  geom_text(
    aes(y = upper, label = round(median, 2)),
    position = position_dodge2(width = 0.9),
    vjust = -0.3
  ) +
  
  scale_fill_manual(values = c(col_blue, col_red)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 19, face = "italic"),
    strip.text  = element_text(size = 19, face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  
  xlab("") +
  ylab(expression(paste("Cor(", hat(y[v]), ", ", y[v], ")"))) +
  
  geom_hline(
    data = Gpred %>% filter(CV == i & metrics == "acc"),
    aes(yintercept = median),
    color = "red"
  ) +
  
  geom_hline(
    data = Gpred %>% filter(CV == i & metrics == "acc"),
    aes(yintercept = upper),
    linetype = "dashed",
    color = "red"
  ) +
  
  geom_hline(
    data = Gpred %>% filter(CV == i & metrics == "acc"),
    aes(yintercept = lower),
    linetype = "dashed",
    color = "red"
  )

getwd()
png(paste0("plot/",i, ".png"), units = "in", width = 10.5, height = 12, res = 300)
print(PLOT)
dev.off()
}



i="Between"
for(i in unique(db$CV)) {
print(i)
PLOT= print(i)

PLOT <- db %>% 
  filter(CV == i & metrics != "acc") %>%
  
  ggplot(aes(x = model_all, y = median, color = model, label = round(median, 2))) +
  
  geom_point(
    position = position_dodge2(width = 0.9),
    size = 10
  ) +
  
  geom_errorbar(
    aes(ymin = lower, ymax = upper),
    position = position_dodge2(width = 0.9),
    width = 0.5
  ) +
  
  facet_grid(
    phenotype ~ paste0("S", scenario),
    scales = "free_x"
  ) +
  
  geom_text(
    aes(y = upper, label = round(median, 2)),
    position = position_dodge2(width = 0.9),
    color = "black",
    vjust = -0.3
  ) +
  
  scale_color_manual(values = c(col_blue, col_red)) +
  
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 19, face = "italic"),
    strip.text  = element_text(size = 19, face = "italic"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
    legend.position = "none"
  ) +
  
  xlab("") +
  ylab(expression(paste("Cor(", hat(y[v]), ", ", y[v], ")"))) +
  geom_hline(data = Gpred %>% filter(median!=0 & CV==i & !metrics=="acc"), aes(yintercept = median), color = "red")+
    geom_hline(data = Gpred %>% filter(median!=0 & CV==i & !metrics=="acc"), aes(yintercept = upper), linetype = "dashed", color = "red")+
    geom_hline(data = Gpred%>% filter(median!=0 & CV==i & !metrics=="acc"), aes(yintercept = lower), linetype = "dashed", color = "red")


png(paste0("plot/",i, "RMSE.png"), units = "in", width = 10.5, height = 12, res = 300)
print(PLOT)
dev.off()
}

# Loop over scenarios to check consistency per scenario




C=C+geom_hline(data=for_all,aes(yintercept=cor),linetype="dotted",color="red")


anova(lm(median ~ model,data=all%>% filter(scenario=="S4" & !model=="Gpred")))






# Calcolo rank per scenario × phenotype
df_rank <- all %>% filter(!model=="Gpred" ) %>%
  group_by(scenario, phenotype,CV) %>%
  mutate(rank = rank(-median, ties.method = "min")) %>%
  ungroup()

df_rank %>% filter(phenotype=="IMF" )
# Grafico linee dei ranghi


ggplot(df_rank %>% filter( model %in% c("Mpred","oMpred","GMpred")), 
aes(x = phenotype, y = rank, color=model,group = model,fill=model)) + geom_line()+
 facet_grid(scenario~CV)

  geom_line(size = 1) +
  #geom_ribbon(
  #    aes(ymin = lower, ymax = upper),alpha = 0.2
  #  )+
  geom_point(size = 3) + facet_grid(phenotype~CV)


ggplot(df_rank %>% filter(model %in% c("oMpred","Mpred","GMpred")), 
aes(x = model, y = rank, group = scenario, color =scenario)) +
  geom_line(size = 1) +
  geom_point(size = 3) + facet_grid(CV~phenotype)




all %>%
  group_by(scenario,phenotype,CV) %>%
  do({
    model <- aov(median ~ model, data = .)
    tidy <- broom::tidy(model)
    tidy
  }) %>% filter(!is.na(statistic)) %>% select(scenario, CV,phenotype,p.value) %>% arrange(p.value)


all %>%
  group_by(scenario) %>%
  do({
    model <- aov(median ~ model,, data = .)
    tidy <- broom::tidy(model)
    tidy
  }) %>% filter(!is.na(statistic)) %>% select(scenario,p.value) %>% arrange(p.value)

all %>% filter(!model=="Gpred") %>%
  group_by(model) %>%
  do({
    model <- aov(median ~ model, data = .)
    tidy <- broom::tidy(model)
    tidy
  }) %>% filter(!is.na(statistic)) %>% select(model,p.value) %>% arrange(p.value)

