vir_count<-read.table("ASVs_counts_copy.tsv",row.names=1,header=T)
meta<-read.csv("meta_cmd.csv", header = T, row.names = 1)
meta$Sample_ID<-rownames(meta)
dim(vir_count)
vir_count2<-t(vir_count)
dim(vir_count2)
vir_count2_hel<-decostand(vir_count2, "hellinger")
head(vir_count2_hel)[1:5,1:5]
vir_count2_hel2<-as.data.frame(vir_count2_hel)
vir_count2_hel2$Sample_ID<-rownames(vir_count2_hel2)
test<-merge(meta,vir_count2_hel2, by="Sample_ID")
drop <- c("Microcystin")
test2 = test[,!(names(test) %in% drop)]
install.packages("devtools")
library(devtools)
devtools::install_github("cran/mvpart")
install_github("cran/MVPARTwrap", force=TRUE)
library(MVPARTwrap)
dim(test2)
spe.hel<-test2[,21:5364]
dim(spe.hel)
env<-test2[,1:20]
doubs.mrt <- mvpart(as.matrix(spe.hel) ~Years ,env,
                    legend=FALSE, margin=0.01, cp=0, xv="pick",
                    xval=nrow(spe.hel), xvmult=100, which=4)

colnames(test[1:25])
head(env[1:2, 1:20])

vc <- ASV_count
vc <- vir_count
meta2 <- meta

head(vc)
tvc <- t(vc)
vc_helli <- decostand(tvc, "hellinger")
dim(meta2)
dim(vc_helli)
vc_hel_keep <- vc_helli[rownames(vc_helli) %in% rownames(meta2),]
dim(vc_hel_keep)
meta_keep <- meta2[rownames(meta2) %in% rownames(vc_hel_keep),]
dim(meta_keep)
doubs.mrt <- mvpart(as.matrix(vc_hel_keep) ~Years ,meta_keep,
                    legend=FALSE, margin=0.01, cp=0, xv="pick",
                    xval=nrow(vc_hel_keep), xvmult=100, which=4)

order(spe.hel) %in% order(vc_hel_keep)

dim(spe.hel)
dim(vc_hel_keep)
    
