rm(list=ls())




dir.create("/home/enrico/Work/prediction_micro_v2/rdata",recursive =TRUE)


cat("First Script for prepare data and create Rdata file to be easy in the next step\n")

cat("in un file Rdata ci da:
    1. le matrici genomiche
    2. il datset pulite piu id animali e sire")


cat("1. Load file with all function alredy\n")

source("functions/function.R")

cat("2. Load packages")

require(tidyverse)
require(data.table)
require(readxl)
require(eimpute)
require(reshape2)

cat("che tempo il microbiota ")
args <- commandArgs(trailingOnly = TRUE)

# Check if a string is passed
if (length(args) < 1) {
  stop("Please pass a string like '1 2 3 4'")
}

# Read the string
input_string <- args[1]

# Print raw input
cat("Input string:", input_string, "\n")

# Split into vector and convert to numeric
period  <- as.numeric(unlist(strsplit(input_string, " ")))
print(period)
# Show result
cat("Parsed values:", paste(period, collapse = ", "), "\n")

#stop()

cat("3.Load all Data and chose effect ")

data=as.data.frame(fread("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/AllData16.csv"))
data[data=="."] = NA
add=names(data)[grep("Sire",names(data))][c(1,2)]

# loop per ogni fenotipo

#In TE:
# 1) models fixed effects included sire (28 levels fitted solely for models including only microbial information), 
# 2) contemporary group (12 levels)
# 3) maternal genetic line (2 levels)
# 4R) the random effects of pen (which here was the physical group of same-sex paternal half-sibs),
# and the residuals.

#In the NU population, 
# 1) fixed effects included sire (28 levels fitted solely for models including only microbial information),
# 2) contemporary group (66 levels), 
#3)  sex (2 levels),
# 4) the random effects of the biological litter where the individual was born, plus the residual error. 

# However, according with francesco i put this effect

croseff=c("CIG_Gender","CIG_Rep2","CIG_Pen","sireGLDS")
pureff=c("PIG_Gender", "PIG_BirthCG","PIG_BirthLitterID","sireGLDS")


# Part two \n make prediction cross breed 
# phenotype / next steps create a list of model 

pbcb=read_excel(("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/models1.xlsx"))
#names(data)

y_pb=pbcb$PB
y_cb=pbcb$CB

cbrm=y_cb[!y_cb%in% names(data) ]
pnrm=y_pb[!y_pb%in% names(data) ]

pbcb=pbcb[!pbcb$CB %in% cbrm ,]
pbcb=pbcb[!pbcb$PB %in% pnrm ,]

y_pb=pbcb$PB
y_cb=pbcb$CB


data=as.data.frame(data)

cb = data[! (is.na(data[,y_cb])),c("animGLDS",croseff,y_cb)]
pb =data[! (is.na(data[,y_pb])),c("animGLDS",pureff,y_pb)]


apply(cb[,croseff],2,function(x) length(unique(x)))
apply(pb[,pureff],2,function(x) length(unique(x)))


cat("create the two kernel maltrix ...\n")

# genomic matrix 
G = f90_to_BLGR("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/geno.RData",
                "/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/map.RData")
G = qcG(G,0.05,0.95,0.95)
K=Ikernel(G)


# microbiome matrix
microb=minfo("/mnt/user_data/enrico/PROGETTI/prediction_micro/DATA_Maschoff/AbundanceTransf_20231123.RData",period)

Km = microb[[1]]
M = microb[[2]]

cat("correlation with matrix\n")

tmp1=reshape2::melt(K)
tmp2=reshape2::melt(Km)
names(tmp1) = c("id","id2","g")
names(tmp2) = c("id","id2","m")
tmp3=merge(tmp1,tmp2)
cat("correlations between two matrix: ",cor(tmp3[,3:4])[1,2],"\n")
rm(tmp1,tmp2)
cat("correlation of diagonal element")
cor(tmp3[tmp3$id==tmp3$id2,c("g","m")])
cat("correlation offdiagonal element")
cor(tmp3[!tmp3$id==tmp3$id2,c("g","m")])
#print(plot(tmp3[!tmp3$id==tmp3$id2,c("g","m")]))

#save(x=tmp3,file="rdata/comparedMvsG.R")
rm(tmp3)


anim=c(pb$animGLDS,cb$animGLDS)
animal_keep <- Reduce(intersect, list(anim, row.names(K), row.names(Km)))

cat("keep only animal with information\nthen use single step\n")

pb = pb[pb$animGLDS %in% animal_keep,]
cb = cb[cb$animGLDS %in% animal_keep,]

cat(nrow(pb),"<=>",nrow(cb),"\n")

pb_anim = pb$animGLDS
cb_anim = cb$animGLDS

tmp=table(data[data$animGLDS %in% pb_anim,]$sireGLDS)

sire=names(tmp)

pass_info_1=list(pb=pb,cb=cb,pb_anim=pb_anim,cb_anim=cb_anim,G=G,M=M,Km=Km,K=K,sire=sire)

save(x=pass_info_1,file="/home/enrico/Work/prediction_micro_v2/rdata/info_r1.R")


