control <- trainControl("repeatedcv", repeats = 10, selectionFunction = "oneSE",returnResamp = "all", savePredictions = "final", search="grid", classProbs=T,summaryFunction = twoClassSummary)
ncompf <- min(0.75*ncol(abund_clean_env_helCLF) - 1, 0.75*ncol(abund_clean_env_helCLF))
set.seed(sample(ncompf))
ncompf <- min(0.75*ncol(abund_clean_env_helCLF) - 1, 0.75*ncol(abund_clean_env_helCLF))
grid <- expand.grid(ncomp=seq(0,ncompf,1))
samp <- sample(1:nrow(grid), 10)
grid[samp,]
fit.pls <- train(bloom ~. , data=abund_clean_env_helCLF, method="pls", metric=metric, tunegrid=samp, trControl=control,tuneLength = 20,returnResamp = "all",savePredictions = "all")


metric="ROC"
control <- trainControl("repeatedcv", repeats = 10, selectionFunction = "oneSE",returnResamp = "all", savePredictions = "final", search="grid", classProbs=T,summaryFunction = twoClassSummary)
ncompf <- min(0.75*ncol(abund_clean_env_helCLF) - 1, 0.75*ncol(abund_clean_env_helCLF))
set.seed(sample(ncompf))
ncompf <- min(0.75*ncol(abund_clean_env_helCLF) - 1, 0.75*ncol(abund_clean_env_helCLF))
grid <- expand.grid(ncomp=seq(0,ncompf,1))
samp <- sample(1:nrow(grid), 500)
ncompT <- data.frame(.ncomp = ncomp=0.75*nrow)
fit.pls <- train(bloom ~. , data=abund_clean_env_helCLF, method="pls", metric=metric, tuneGrid=ncompT, trControl=control,tuneLength = 20,returnResamp = "all",savePredictions = "all")


#pls
metric="RMSE" 
control <- trainControl("repeatedcv", repeats = 10, selectionFunction = "oneSE",search="grid")
set.seed(sample(ncompf))
ncompf <- min(0.75*nrow(abund_clean_env_helCLF) - 1, 0.75*nrow(abund_clean_env_helCLF))
grid <- expand.grid(ncomp=seq(0,ncompf,1))
samp <- sample(1:nrow(grid), 45)
ncompT <- data.frame(.ncomp = c(min(samp):max(samp)))
fit.pls <- train(bloom ~. , data=abund_clean_env_helCLF, method="pls", metric=metric, tuneGrid=ncompT, trControl=control,tuneLength = 20,returnResamp = "all",savePredictions = "all")
plot(fit.pls)
fit.pls$bestTune
plot(varImp(fit.pls), 10, main = "PLS")















### MATNEL example https://www.flutterbys.com.au/stats/tut/tut15.2.html ###
coenocline <- function(x,A0,m,r,a,g, int=T, noise=T) {
             b <- a/(a+g)
              d <- (b^a)*(1-b)^g
              cc <- (A0/d)*((((x-m)/r)+b)^a)*((1-(((x-m)/r)+b))^g)
              if (noise) {n <- A0/10; n[n<0]<-0; cc<-cc+rnorm(length(cc),0,n)}
              cc[cc<0] <- 0
              cc[is.na(cc)]<-0
              if (int) cc<-round(cc,0)
              cc
     }
   
   dummy <- function(x) {
       nms <- colnames(x)
       ff <- eval(parse(text=paste("~",paste(nms,collapse="+"))))
       mm <- model.matrix(ff,x)
       nms <- colnames(mm)
       mm <- as.matrix(mm[,-1])
       colnames(mm) <- nms[-1]
       mm
     }

   set.seed(1)
    x <- seq(0,50,l=10)
    n <- 10
    sp1<-coenocline(x=x,A0=5,m=0,r=2,a=1,g=1,int=T, noise=T)
    sp2<-coenocline(x=x,A0=70,m=7,r=30,a=1,g=1,int=T, noise=T)
    sp3<-coenocline(x=x,A0=50,m=15,r=30,a=1,g=1,int=T, noise=T)
    sp4<-coenocline(x=x,A0=7,m=25,r=20,a=0.4,g=0.1,int=T, noise=T)
    sp5<-coenocline(x=x,A0=40,m=30,r=30,a=0.6,g=0.5,int=T, noise=T)
    sp6<-coenocline(x=x,A0=15,m=35,r=15,a=0.2,g=0.3,int=T, noise=T)
    sp7<-coenocline(x=x,A0=20,m=45,r=25,a=0.5,g=0.9,int=T, noise=T)
    sp8<-coenocline(x=x,A0=5,m=45,r=5,a=1,g=1,int=T, noise=T)
    sp9<-coenocline(x=x,A0=20,m=45,r=15,a=1,g=1,int=T, noise=T)
    sp10<-coenocline(x=x,A0=30,m=50,r=5,a=1,g=1,int=T, noise=T)
    X <- cbind(sp1, sp10,sp9,sp2,sp3,sp8,sp4,sp5,sp7,sp6)
    #X<-X[c(1,10,9,2,3,8,4,5,7,6),] 
      colnames(X) <- paste("Sp",1:10,sep="")
    rownames(X) <- paste("Site", c(1,10,9,2,3,8,4,5,7,6), sep="")
    X <- X[c(1,4,5,7,8,10,9,6,3,2),]
    data <- data.frame(Sites=factor(rownames(X),levels=rownames(X)), X)
data
    
set.seed(1)
     Site <- gl(10,1,10,lab=paste('Site',1:10, sep=""))
     Y <- matrix(c(
       6.1,4.2,101325,2,
       6.7,9.2,101352,510,
       6.8,8.6,101356,546,
       7.0,7.4,101372,758,
       7.2,5.8,101384,813,
       7.5,8.4,101395,856,
       7.5,0.5,101396,854,
       7.0,11.8,101370,734,
       8.4,8.2,101347,360,
       6.2,1.5,101345,356
         ),10,4, byrow=TRUE)
     colnames(Y) <- c('pH','Slope', 'Pressure', 'Altitude')
     Substrate <- factor(c('Quartz','Shale','Shale','Shale','Shale','Quartz','Quartz','Shale','Quartz','Quartz'))
     enviro <- data.frame(Site,Y,Substrate)
enviro     

data.dist <- vegdist(wisconsin(sqrt(data[,-1])),"bray")
data.dist
dim(data.dist)

str(enviro)
enviro1 <- within(enviro, Substrate <- as.numeric(Substrate))
env.dist <- vegdist(decostand(enviro1[,-1],"standardize"),"euclid")
env.dist
dim(env.dist)

plot(data.dist, env.dist)

data.mantel <- mantel(data.dist, env.dist, perm=1000)
data.mantel

enviro1 <- dummy(enviro[,-1])

library(car)
vif(lm(1:nrow(enviro1) ~ pH+Slope+Pressure+Altitude+SubstrateShale, data=data.frame(enviro1)))
vif(lm(1:nrow(enviro1) ~ pH+Slope+Altitude+SubstrateShale, data=data.frame(enviro1)))

data.bioenv <- bioenv(wisconsin(sqrt(data[,-1])), decostand(enviro1[,-3],"standardize"))
data.bioenv







library(mlbench);
library(caret);
library(caretEnsemble);
library(tidyverse);
library(magrittr);
library(purrrlyr)
library(lubridate);
library(ape);
library(extrafont);
library(metagenomeSeq);
library(boral);
library(corrplot);
library(igraph);
library(RColorBrewer);
library(Cairo);
library(data.table);
library(plyr);
library(gridExtra);
library(ggplot2);
library(reshape2);
library(Matrix);
library(AICcmodavg);
library(compositions);
library(zCompositions);
library(rfUtilities);
library(randomForest)
library(devtools)
install_github("zdk123/SpiecEasi")
library(SpiecEasi)
library("remotes")
install.packages("corrplot")
library("corrplot")
install.packages("GUniFrac")
library("GUniFrac")
install.packages("corrr")
library("corrr")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("microbiome")
devtools::install_github("vmikk/metagMisc")
setwd("/Users/kiristern/Downloads")
vir_count<-read.csv("viral_count_for_bac.csv",row.names=1,header=T)
head(vir_count[1:5,1:5])
Bact<-read.table("Champ_ASVs_counts2.txt",row.names=1,header=T)
head(Bact[1:5,1:5])
gLV<-read.csv("gLV_table_FILT_10_2.csv", header=T)
meta<-read.csv("metadata3.csv",row.names=1,header=T)
meta_2<-read.csv("met_gLV.csv", row.names=1,header=T)
met_bac<-read.csv("met_asv_Bact.csv", row.names = 1, header=T)
met_vir<-read.csv("meta_vir.csv", row.names = 1, header=T)
taxbact<-as.matrix(read.table("taxonomy_greengenes.txt", header=T, row.names=1))
#Hellinger transformation
vir_count.no0 = vir_count[ ,colSums(vir_count)!=0, ]
vir_count_hel<-decostand(vir_count.no0, method="hellinger")
vir_pa <- decostand(vir_count,"pa")
vir_pa_10 <- vir_count[,colSums(vir_pa)/nrow(vir_pa) >= 0.10]
names(vir_pa_10)
vir_cleaned_hel <- vir_count_hel[,names(vir_pa_10)]
vir_cleaned <- vir_count.no0[,names(vir_pa_10)]
bac_count.no0 = Bact[ ,colSums(Bact)!=0, ]
bac_count_hel<-decostand(bac_count.no0, method="hellinger")
bac_pa <- decostand(Bact,"pa")
bac_pa_10 <- Bact[,colSums(bac_pa)/nrow(bac_pa) >= 0.10]
names(bac_pa_10)
bac_cleaned_hel <- bac_count_hel[,names(bac_pa_10)]
bac_cleaned <- bac_count.no0[,names(bac_pa_10)]
bac_phy<-otu_table(t(bac_cleaned),taxa_are_rows=T)
sample_infobac<-sample_data(met_bac)
TAX = tax_table(taxbact)
bac_physeq<-phyloseq(bac_phy,TAX, sample_infobac)
colnames(tax_table(bac_physeq)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cyan=subset_taxa(bac_physeq, Phylum=="p__Cyanobacteria")
#cyan_1<-phyloseq_to_df(cyan, addtax = F)
cyan_1 <- cyan %>% otu_table() %>% as.data.frame
head(cyan_1[1:5,1:5])
#cyano_2<-cyan_1 %>% remove_rownames %>% column_to_rownames(var="OTU")
head(cyano_2[1:5,1:5])
#cyano_3<-t(cyano_2)
cyano_3<-t(cyan_1)
head(cyano_3[1:5,1:5])
#Coinertia, Procrustes,Mantel (how both matrices or df are correlated)
dudi.vir <- dudi.pca(vir_cleaned_hel, scale = FALSE, scan=FALSE)
dudi.cyan <- dudi.pca(cyano_helli, scale = FALSE, scan=FALSE)
dudi.bact<-dudi.pca(bac_cleaned_hel,scale = FALSE, scan=FALSE )
coinert<-coinertia(dudi.vir, dudi.cyan,scannf = FALSE, nf = 2 )
coinert<-coinertia(dudi.vir, dudi.bact,scannf = FALSE, nf = 2 )
coinert$RV
RV.rtest(cyano_helli, vir_cleaned_hel, nrepet = 999)
RV.rtest(bac_cleaned_hel, vir_cleaned_hel, nrepet = 999)
protest(dudi.bact, dudi.vir)
dist_vir<-sqrt(vegdist(vir_cleaned, method = "bray"))
dist_bac<-sqrt(vegdist(bac_cleaned, method = "bray"))
dist_cyano<-sqrt(vegdist(cyano_3,method = "bray"))
mantel(dist_vir,dist_bac, method = "pearson")
mantel(dist_vir,dist_cyano, method = "pearson")
protest(dist_vir,dist_cyano)
protest(dist_vir,dist_bac)














