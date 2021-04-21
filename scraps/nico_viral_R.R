#install.packages("mlbench")
#install.packages("caret")
#install.packages("caretEnsemble")
#install.packages("tidyverse")
#install.packages("magrittr")
#install.packages("purrrlyr")
#install.packages("lubridate")
#install.packages("ape")
#install.packages("extrafont")
#install.packages("AICcmodavg")
#install.packages("zCompositions")
#install.packages("rfUtilities")
#install.packages("randomForest")


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

setwd("/Users/nico/Desktop/Collaboration/Etudiants/Kiri/")
vir_count<-read.csv("viral_count_for_bac.csv",row.names=1,header=T)
Bact<-read.table("Champ_ASVs_counts2.txt",row.names=1,header=T)
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
dist_cyano<-sqrt(vegdist(cyano_helli,method = "bray"))
mantel(dist_vir,dist_bac)
mantel(dist_vir,dist_cyano)

#Phyloseq
vir_phy<-otu_table(t(vir_count),taxa_are_rows=T)
sample_info<-sample_data(met_vir)
viral_physeq<-phyloseq(vir_phy,sample_info)
vir_filt<-filter_taxa(viral_physeq, function(x) sum(x > 3) > (0.1*length(x)), TRUE)

vir_phy2<-otu_table(t(vir_cleaned_hel),taxa_are_rows=T)
sample_infobac<-sample_data(met_vir)
vir_physeq2<-phyloseq(vir_phy2,sample_infobac)
vir_helli<-phyloseq_to_df(vir_physeq2, addtax = F)
vir_hel<-vir_helli %>% remove_rownames %>% column_to_rownames(var="OTU")
head(vir_hel[1:5,1:5])
vir_helli<-t(vir_hel)
head(vir_helli[1:5,1:5])

bac_phy<-otu_table(t(Bact),taxa_are_rows=T)
sample_infobac<-sample_data(met_bac)
TAX = tax_table(taxbact)
bac_physeq<-phyloseq(bac_phy,TAX, sample_infobac)
bac_filt<-filter_taxa(bac_physeq, function(x) sum(x > 3) > (0.1*length(x)), TRUE)
colnames(tax_table(bac_filt)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cyan=subset_taxa(bac_filt, Phylum=="p__Cyanobacteria")

#Extract cyano from transformed table using phyloseq
bac_phy<-otu_table(t(bac_cleaned_hel),taxa_are_rows=T)
sample_infobac<-sample_data(met_bac)
TAX = tax_table(taxbact)
bac_physeq<-phyloseq(bac_phy,TAX, sample_infobac)
colnames(tax_table(bac_physeq)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cyan=subset_taxa(bac_physeq, Phylum=="p__Cyanobacteria")
cyan_helli<-phyloseq_to_df(cyan, addtax = F)
cyano_hel<-cyan_helli %>% remove_rownames %>% column_to_rownames(var="OTU")
head(cyano_hel[1:5,1:5])
cyano_helli<-t(cyano_hel)
head(cyano_helli[1:5,1:5])

#cyano_hel2<-cyan_helli %>% remove_rownames %>% column_to_rownames(var="OTU")
#vir_hel2<-vir_helli %>% remove_rownames %>% column_to_rownames(var="OTU")
#cyano_hel_fin<-t(cyano_hel2)
#vir_hel_fin<-t(vir_hel2)

#corr
f <- function (x, y) rcorr(x, y, method = "spearman")
tab <- outer(vir_cleaned_hel, cyano_helli, Vectorize(f))
as.data.frame.table(tab)

#network
se.hmp2 <- spiec.easi(list(cyan, vir_physeq2), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))
dtype <- c(rep(1,ntaxa(cyan)), rep(2,ntaxa(viral_physeq)))
plot(adj2igraph(getRefit(se.hmp2)), vertex.color=dtype+1, vertex.size=9)
betaMat=as.matrix(symBeta(getOptBeta(se.hmp2)))


#CLRtransfromation
abund<-vir_count
abund.no0=abund[,colSums(abund)!=0,]
abund.count<-t(cmultRepl(t(abund.no0),method="CZM",output="p-counts"))
abund.prop.rows<-decostand(abund.count,"total",MARGIN=1)
abund.clr<-t(apply(abund.prop.rows,1,function(x){log(x)-mean(log(x))}))
abund.clr<-as.data.frame(abund.clr)

abund_bc<-Bact
abundcb.no0=abund_bc[,colSums(abund_bc)!=0,]
abundbc.count<-t(cmultRepl(t(abundcb.no0),method="CZM",output="p-counts"))
abundbc.prop.rows<-decostand(abundbc.count,"total",MARGIN=1)
abundbc.clr<-t(apply(abundbc.prop.rows,1,function(x){log(x)-mean(log(x))}))
abundbc.clr<-as.data.frame(abundbc.clr)

abundbc.clr$Samples_ID<-rownames(abund_bc)
abund.clr$Samples_ID<-rownames(abund)

vir_Bac<-merge(abundbc.clr,abund.clr,by="Samples_ID")

abund_dolicho<-abundbc.clr[,c("ASV_11","ASV_16","ASV_18","ASV_50","ASV_85","ASV_110","ASV_129","ASV_280","ASV_292","ASV_356","ASV_379","ASV_683","ASV_735","ASV_836",
"ASV_890","ASV_1111","ASV_1114","ASV_1184","ASV_1243","ASV_1397","ASV_1488","ASV_1506","ASV_1812","ASV_1831","ASV_1881","ASV_2167","ASV_2185","ASV_2303","ASV_2628","ASV_2850","ASV_2867","ASV_2985","ASV_1675","ASV_5900","ASV_3991")]
abund_microcystis<-abundbc.clr[,c("ASV_7","ASV_8","ASV_143","ASV_859","ASV_1928")]
#Merging Virus and Cyano and doing spiecEasi network (with one table this time)


