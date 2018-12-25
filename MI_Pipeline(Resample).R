#command to clear console
cat("\014")  
#command to clear all obj
rm(list = ls())

############## ONLY CHANGE THESE VALUES ##################

#set working directory
setwd("/ddn/gs1/home/sankara2/Data/Pipetest3")
#ANOVA clustering FILENAME 
clusters = "ClumpResult_ForClustering_tox21-er-bla-agonist-p2_ratio.txt.dat"
meta = "protocol_meta_KRSmod.txt"

#Only need to install *ONCE*
#install.packages("dplyr")

################ END HERE ################################
#parallel 
library(foreach)
library(doParallel)
registerDoParallel(cores = 60)

#Output file prefix
library(dplyr)
x = (strsplit(clusters, ".", fixed = TRUE)[[1]][1] %>% strsplit("_", fixed = TRUE) %>% unlist()) [3:4] 
out.name = paste(x[1], x[2], sep = "_")


#START w/ ANOVA cluster .txt.dat file
#create .txt file accounting for all observations
x = colnames(read.delim(clusters, header=TRUE));
data.mat = read.delim(clusters, quote= "", header=TRUE);
data.mat = as.data.frame(sapply(data.mat, function(x) gsub("\"", "", x)))
colnames(data.mat) = x
dir.create("./Data")
setwd("./Data")
clusters = substr(clusters, 1, regexpr("dat", clusters, fixed = TRUE)[1] - 2)
write.table(data.mat, file= clusters, sep = "\t", row = FALSE)
rm(list = setdiff(ls(), c("meta", "clusters", "out.name")))


#Write-out CASRN w/ homogeneous profiles
data.mat = read.delim(clusters,  header=TRUE);
attach(data.mat)
t1 = summarise(group_by(data.mat, CAS), count= n())
t2 = summarise(group_by(data.mat, CAS, clump), count= n()) %>% subset(clump==1, select= c(CAS, count))
data.filt = data.mat[data.mat$CAS %in% intersect(t1, t2)$CAS,]
detach(data.mat)
CASRN = unique(data.filt$CAS, MARGIN= 1)
write.table(CASRN, file = "CASRN.txt", sep ="\t", row=FALSE)
rm(list = setdiff(ls(), c("meta", "clusters", "out.name")))




#Obtain unique CASRN IDs and generate pair-wise combinations
CASRN = read.delim("CASRN.txt", header= T)[[1]]
pwarray = combn(CASRN, 2)


#Reading data... NEED metadata
data.mat = read.delim(clusters,  header=TRUE);
meta.dat = subset(read.delim(paste("../", meta,sep = ""), header = TRUE),select=c(SAMPLE_DATA_TYPE, ASSAY_SD));

#Merging negative controls SD
if ("ratio" %in% meta.dat$SAMPLE_DATA_TYPE){
  data.mat$ASSAY_SD = meta.dat$ASSAY_SD[3] ;
} else{
  data.mat = merge(data.mat, meta.dat, by.x= "Sample.Data.Type", by.y= "SAMPLE_DATA_TYPE");
}

#Write subset (CASID, concentrations, responses, ASSAY_SD) for homogeneous profiles
conc.ind = grep("conc", colnames(data.mat)) [1] : (grep("conc", colnames(data.mat)) [1] + 14)
resp.ind = grep("resp", colnames(data.mat)) [1] : (grep("resp", colnames(data.mat)) [1] + 14)
CAS.ind = grep("^CAS$", colnames(data.mat))
ASSAY_SD.ind = grep("^ASSAY_SD$", colnames(data.mat))

data.sub = NULL;
for (i in 1:length(CASRN)){
  pos = grep(as.character(CASRN[i]), as.character(data.mat[,CAS.ind]))
  data.sub = rbind(data.sub, data.mat[pos, c(CAS.ind, conc.ind, resp.ind, ASSAY_SD.ind)])  
}
write.table(data.sub, paste(out.name, "_0.txt", sep = ""), sep = "\t", row = FALSE)


#Write 30 PERMUTED subsets (CASID, concentrations, responses, ASSAY_SD) for homogeneous profiles
# 0 = original conc-resp data; >0 = permuted conc-resp data
data.sub = read.delim(paste(out.name, "_0.txt", sep = ""), header = TRUE)
foreach(k= 1:30) %dopar% {
  data.sub[,17:31]= t(apply(data.sub[,17:31], MARGIN = 1, function(x) sample(x, replace = TRUE)))
  write.table(data.sub,paste(out.name, "_", k, ".txt", sep = ""), sep = "\t", row = F)
}


library(splines)
#MI calculations as a f(x)
info <- function(arg1, arg2, data){
  data.mat = data
  #Positions of chemicals
  pos1 <- grep(arg1, as.character(data.mat[,1]))
  pos2 <- grep(arg2, as.character(data.mat[,1]))
  
  
  #----------SPLINE----------
  
  c = rbind(log10(data.mat[pos1, 2:16]*10^6), log10(data.mat[pos2, 2:16]*10^6))
  range.mat = t(apply(c, 1, range, na.rm= T))
  #Concentration region of overlap
  range = c(max(range.mat[,1]), min(range.mat[,2]))
  
  #Interpolating spline estimates of responses @ 15 equally spaced concentrations
  conc.mat = matrix(NA, length(pos1) + length(pos2), 15)
  resp.mat = matrix(NA, length(pos1) + length(pos2), 15)
  for (i in 1:(length(pos1) + length(pos2))){
    
    if(i <= length(pos1)){
      ispl <- interpSpline(c[i,], data.mat[pos1[i], 17:31], na.action = na.omit)
      nspl <- predict(ispl, seq(range[1], range[2], length.out = 15))
      resp.mat[i,] = nspl$y
      conc.mat[i,] = 10^(nspl$x)
    }
    else{
      ispl <- interpSpline(c[i, ], data.mat[pos2[i -  length(pos1)], 17:31], na.action = na.omit)
      nspl <- predict(ispl, seq(range[1], range[2], length.out = 15))
      resp.mat[i,] = nspl$y
      conc.mat[i,] = 10^(nspl$x)
    }
  }
  
  
  #-----------------Data Discretization------------------
  #f1 -> marginal profile of chemical X
  #f2 -> marginal profile of chemical Y
  f1=matrix(NA, nrow=length(pos1), ncol=15);
  f2=matrix(NA, nrow=length(pos2), ncol=15); 
  
  for(i in 1:15){
    for(j in 1:length(pos1)){
      DETECTIONLIMIT = data.mat$ASSAY_SD[pos1[j]]
      response = resp.mat[j, i]
      
      f1[j, i] = trunc(response/(3 * DETECTIONLIMIT))
    }
    
    for (k in 1:length(pos2)){
      DETECTIONLIMIT = data.mat$ASSAY_SD[pos2[k]]
      response = resp.mat[k + length(pos1), i]
      

      f2[k, i] = trunc(response/(3 * DETECTIONLIMIT))
    }
    
  }
  
  
  
  #-------Estimating MI-----------
  
  #f12 -> joint profile of chemical X and Y.
  f12 = rbind(f1, f2);
  
  #Number of concentrations sampled (= 15)
  denom = ncol(f12);
  
  #distinct obs. for marginal and joint distributions
  uniquef1 = unique(f1, MARGIN=2);
  uniquef2 = unique(f2, MARGIN=2);
  uniquef12 = unique(f12, MARGIN=2);
  
  
  #evaluation of k(x, *) term and k(x, *)log2[k(x, *)]
  p1=vector(); term1=vector();
  ent1=vector();
  for(i in 1:ncol(uniquef1)){
    p1[i]    = sum(apply(f1, 2, identical, uniquef1[,i]));
    term1[i] = p1[i]*log2(p1[i])
    ent1[i] = p1[i]/denom * log2(p1[i]/denom)
  }
  
  #evaluation of k(*, y) term and k(*, y)log2[k(*, y)]
  p2=vector(); term2=vector(); 
  ent2=vector();
  for(i in 1:ncol(uniquef2)){
    p2[i]    = sum(apply(f2, 2, identical, uniquef2[,i]));
    term2[i] = p2[i]*log2(p2[i])
    ent2[i] = p2[i]/denom * log2(p2[i]/denom)
  }
  
  #evaluation of k(x,y) term and k(x,y)log2[k(x,y)]
  p12=vector(); term12=vector(); 
  for(i in 1:ncol(uniquef12)){
    p12[i]    = sum(apply(f12, 2, identical, uniquef12[,i]));
    term12[i] = p12[i]*log2(p12[i])
  }
  
  #Calculation of MI and scaled MI; n = # of concentrations
  MutInfo = log2(denom) + 1/denom*(sum(term12) - sum(term1) - sum(term2));
  MIScaled = MutInfo/min(-sum(ent1), -sum(ent2))
  
  #mean responses @ 15 equally spaced concentrations for each chemical
  chemX = apply(resp.mat[1:length(pos1),], 2, mean)
  chemY = apply(resp.mat[(length(pos1)+1):(length(pos1)+length(pos2)), ], 2, mean)
  
  #Pearson and Spearman correlation estimates
  pear = cor.test(chemX, chemY, method = "pearson")
  spear = cor.test(chemX, chemY, method = "spearman")
  
  
  return(c(MutInfo, MIScaled, -sum(ent1), -sum(ent2), ncol(uniquef1), ncol(uniquef2), pear$estimate, spear$estimate))
  
}



#-------------Main Program----------
foreach(i= 0:30) %dopar% {
  data.mat = read.delim( paste(out.name, "_", i , ".txt", sep = ""), header=TRUE);
  #print(paste(out.name, "_", i , ".txt", sep = ""))
  
  #Calculate and store CAS IDs, and MI, Pearson, and Spearman #s
  t1= vector();
  t2= vector();
  dat= matrix(data = NA, nrow = length(pwarray[1,]), ncol = 8);
  
  for (a in 1: length(pwarray[1,])){
    t1[a] = as.character(pwarray[1,a])
    t2[a] = as.character(pwarray[2,a])
    dat[a,] = try(info(t1[a], t2[a], data.mat) )
  }
  
  #Create and write dataframe with CAS IDs, and MI, Pearson, and Spearman
  data.stats = as.data.frame(cbind(t1, t2, dat))
  colnames(data.stats) <- c("CASX", "CASY", "MI",  "MIScaled", "HX", "HY", "|X|", "|Y|", "Pearson", "Spearman")
  write.table(data.stats, file= paste(out.name,"_Spline_", i ,".txt", sep = ""), sep = "\t", row = FALSE)
}





#Total distribution of permuted estimates (MI, Pearson, Spearman)
perm.MI = NULL
perm.pearson = NULL
perm.spearman = NULL

for (i in 1:30){
  dat <- read.delim(paste(out.name,"_Spline_", i ,".txt", sep = ""), header = T)
  perm.MI = append(perm.MI, dat$MI)
  perm.pearson = append(perm.pearson, dat$Pearson)
  perm.spearman = append(perm.spearman, dat$Spearman)
}
write.table(data.frame(perm.MI, perm.pearson, perm.spearman), file= "./Permuted_Est.txt", sep="\t", row.names = F)
rm(list = setdiff(ls(), c("meta", "clusters", "out.name")))



#-----Estimate Probabilities------------
dat <- read.delim(paste(out.name, "_Spline_0.txt", sep = ""), header = T)
perm <- read.delim("./Permuted_Est.txt", header = T)


finalMatrix <- foreach (i = 1:dim(dat)[1], .combine = rbind, .inorder = TRUE) %dopar% {
  MI.P = sum(perm$perm.MI > dat$MI[i]) / dim(perm)[1]
  
  if (dat$Pearson[i] < 0)
    Pearson.P = sum(perm$perm.pearson < dat$Pearson[i]) / dim(perm)[1]
  else
    Pearson.P = sum(perm$perm.pearson > dat$Pearson[i]) / dim(perm)[1]
  
  if (dat$Spearman[i] < 0)
    Spearman.P = sum(perm$perm.spearman < dat$Spearman[i]) / dim(perm)[1]
  else
    Spearman.P = sum(perm$perm.spearman > dat$Spearman[i]) / dim(perm)[1]
  
  c(MI.P, Pearson.P, Spearman.P)
}
colnames(finalMatrix) <- c("MI.P", "Pearson.P","Spearman.P")
write.table(data.frame(dat,finalMatrix), file= paste(out.name, "_Spline_F.txt", sep = ""), sep = "\t",
            row.names = FALSE)
rm(list = setdiff(ls(), c("meta", "clusters", "out.name")))

#------Plots-----
dir.create("../Graphics")
setwd("../Graphics")

#Permuted distributions
perm <- read.delim("../Data/Permuted_Est.txt", header = T)

jpeg(file="Perm_MI.jpg", width = 6, height = 4, units ="in", res = 100)
hist(perm$perm.MI, main = "Distribution of Permuted MI", xlab = "MI (bits)")
dev.off()

jpeg(file="Perm_Pearson.jpg", width = 6, height = 4, units ="in", res = 100)
hist(perm$perm.pearson, main = "Distribution of Permuted Pearson Correlations", xlab = "Pearson Correlation Coefficient (r)")
dev.off()

jpeg(file="Perm_Spearman.jpg", width = 6, height = 4, units ="in", res = 100)
hist(perm$perm.spearman, main = "Distribution of Permuted Spearman Correlations", xlab = "Spearman Correlation Coefficient (p)")
dev.off()


#Scatter plots
dat <- read.delim(paste("../Data/",out.name, "_Spline_F.txt", sep = ""), header = T)

jpeg(file="Pearson_MI_Scatter.jpg", width = 6, height = 6, units ="in", res = 200)
plot(dat$MI, dat$Pearson, cex = .1, main = "Scatterplot of MI vs. Pearson Correlations", xlab = "MI (bits)", 
     ylab = "Pearson Correlation Coefficient (r)")
dev.off()

jpeg(file="Spearman_MI_Scatter.jpg", width = 6, height = 6, units ="in", res = 200)
plot(dat$MI, dat$Spearman, cex = .1, main = "Scatterplot of MI vs. Spearman Correlations", xlab = "MI (bits)", 
     ylab = "Spearman Correlation Coefficient (p)")
dev.off()


#Sample graphs of non-linear relationships
sample_size <- function(x){
  if (x > 10)
    return(10)
  else 
    return(x)
}

data.mat = read.delim(paste("../Data/", out.name, "_0.txt", sep = ""), header = TRUE);
#Pearson
dir.create("./Pearson")
setwd("./Pearson")
temp = subset(dat, dat$MI.P < 0.05 & dat$Pearson.P > 0.05)
temp = sample_n(temp, size = sample_size(nrow(temp)))

if (nrow(temp) > 0){
  for (i in 1:nrow(temp)){
    jpeg(filename = paste("Plot_", i, ".jpg", sep = ""),  width = 10, height = 8, units = "in", res = 150)
    pos1 <- grep(temp$CASX[i], as.character(data.mat[,1]))
    pos2 <- grep(temp$CASY[i], as.character(data.mat[,1]))
    plot(1:10, 1:10, ylim=c(-300,300), xlim=c(0.0001, 100), log="x", col="white", 
         xlab= paste("conc [", expression(mu),"M]", sep = ""), ylab = "resp", main = paste(temp$CASX[i], "vs", temp$CASY[i]))
    for(k in 1:length(pos1)){
      points(data.mat[pos1[k], 2:16]*10^6, data.mat[pos1[k], 17:31], pch=19)
      lines(data.mat[pos1[k], 2:16]*10^6, data.mat[pos1[k], 17:31], lwd=2)
    }
    for(m in 1:length(pos2)){
      points(data.mat[pos2[m], 2:16]*10^6, data.mat[pos2[m], 17:31], pch=19, col="red")
      lines(data.mat[pos2[m], 2:16]*10^6, data.mat[pos2[m], 17:31], lwd=2, col="red")
    }
    
    mtext(paste("MutInfo=", signif(temp$MI[i], digits=4)), line=-2, at=0.01)
    mtext(paste("Pearson=", signif(temp$Pearson[i], digits=4)), line= -3, at=0.01)
    mtext(paste("MI p-value=", signif(temp$MI.P[i], digits=4)), line=-4, at=0.01)
    mtext(paste("Pearson p-value=", signif(temp$Pearson.P[i], digits=4)), line=-5, at=0.01)
    dev.off()
  }
}

#Spearman
dir.create("../Spearman")
setwd("../Spearman")
temp = subset(dat, dat$MI.P < 0.05 & dat$Spearman.P > 0.05)
temp = sample_n(temp, size = sample_size(nrow(temp)))

if(nrow(temp) > 0){
  for (i in 1:nrow(temp)){
    jpeg(filename = paste("Plot_", i, ".jpg", sep = ""),  width = 10, height = 8, units = "in", res = 150)
    
    pos1 <- grep(temp$CASX[i], as.character(data.mat[,1]))
    pos2 <- grep(temp$CASY[i], as.character(data.mat[,1]))
    plot(1:10, 1:10, ylim=c(-300,300), xlim=c(0.0001, 100), log="x", col="white",
         xlab= paste("conc [", expression(mu),"M]", sep = ""), ylab = "resp",  main = paste(temp$CASX[i], "vs", temp$CASY[i]))
    for(k in 1:length(pos1)){
      points(data.mat[pos1[k], 2:16]*10^6, data.mat[pos1[k], 17:31], pch=19)
      lines(data.mat[pos1[k], 2:16]*10^6, data.mat[pos1[k], 17:31], lwd=2)
    }
    for(m in 1:length(pos2)){
      points(data.mat[pos2[m], 2:16]*10^6, data.mat[pos2[m], 17:31], pch=19, col="red")
      lines(data.mat[pos2[m], 2:16]*10^6, data.mat[pos2[m], 17:31], lwd=2, col="red")
    }
    
    mtext(paste("MutInfo=", signif(temp$MI[i], digits=4)), line=-2, at=0.01)
    mtext(paste("Spearman=", signif(temp$Spearman[i], digits=4)), line= -3, at=0.01)
    mtext(paste("MI p-value=", signif(temp$MI.P[i], digits=4)), line=-4, at=0.01)
    mtext(paste("Spearman p-value=", signif(temp$Spearman.P[i], digits=4)), line=-5, at=0.01)
    dev.off()
  }
}









 