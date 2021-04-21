#R codes used
library(vegan)
library(mvpart)
library(rpart)
library(rdaTest)
library(labdsv)
library(plyr)
library(MASS)
library(phyloseq)
library(plotrix)

#data_tables

abund_table<-read.table("ASVs_counts_copy.tsv",row.names=1,check.names=FALSE,header=T)
meta_table1<-read.csv("metadata_w_cyano.csv",header=T,row.names=1,check.names=FALSE)

----------------------------------------
  
  #Test of collinearity
  
meta_table1 = log(meta_table1+1)
env.pearson<-cor(meta_table1)
round(env.pearson, 2)
corvif(meta_table1) #http://www.highstat.com/Book2/HighstatLibV6.R

#To remove the colinear variables

----------------------------------------
  
  #RDA
  
#Transform the abundance table: Hellinger transformation is performed here to avoid the use of double-zeros as indications of resemblance.
vir_abund_helli 

#Standardize the environmental variables table: the environmental data could be in different units, then a standardization is required to perform RDA
meta_table2=meta_table1
meta_table2<-decostand(meta_table1, method="standardize")

rda_a <- rda(sp.asv.rm ~ ., data=env_keep)
scrs<-scores(rda_a,display=c("sp","wa","lc","bp","cn"),scaling=2)
#Extract site data first
df_sites<-data.frame(scrs$sites)
colnames(df_sites)<-c("x","y")
#Draw sites
p<-ggplot()
p<-p+geom_point(data=df_sites,aes(x,y))
p<-p+labs(x="RDA1(x%)",y="RDA2(x%)") #x values could be taken from summary(rda_a)
multiplier <- vegan:::ordiArrowMul(scrs$biplot)
df_arrows<- scrs$biplot*multiplier
colnames(df_arrows)<-c("x","y")
df_arrows=as.data.frame(df_arrows)
p<-p+geom_segment(data=df_arrows, aes(x = 0, y = 0, xend = x, yend = y),
                  arrow = arrow(length = unit(0.2, "cm")),color="black",alpha=0.6)
p
p2<-p+geom_text(data=as.data.frame(df_arrows*1.3),aes(x, y, label = rownames(df_arrows),size=3,fontface="bold"),color="black",alpha=0.6, size=2, fontface="bold")
p2
# Draw species
df_species<- as.data.frame(scrs$species)
colnames(df_species)<-c("x","y")


----------------------------------------
  
  
  #Regression_tree
  
  # time.tree <- mvpart(as.matrix(abund_table2) ~meta_table1$variable ,meta_table1,
  #                     legend=FALSE, margin=0.01, cp=0, xv="pick",
  #                     xval=nrow(abund_table2), xvmult=100, which=4, big.pts=T, bars=F)

rpart.pca(time.tree)

----------------------------------------
  
#K-means partitioning
  
cKM <- cascadeKM(vir_abun_removed, inf.gr=2, sup.gr=10, iter = 999, criterion = 'calinski')
summary(cascadeKM(vir_abun_removed, inf.gr=2, sup.gr=10, iter = 999, criterion = 'calinski'))
#examine the output file struture 
summary(resultat.cascadeKM)

#obtain SCE and calinski index to selectc the optimal partition.
#size=number of object in a groups
resultat.cascadeKM$partition

#Graph of partition
output1 = plot(resultat.cascadeKM)
#Regroup as possible object by group
output2 = plot(resultat.cascadeKM, sortg=TRUE)

----------------------------------------
  
  #Phyloseq analysis
  
otu_table = import_biom("viral_physeq")
mapping_file = import_qiime_sample_data("")
tree_file = read_tree_greengenes("")
data = merge_phyloseq(otu_table,mapping_file,tree_file)
colnames(tax_table(Smile)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Transform to relative abundance
SmileT<- transform_sample_counts(Smile, function(x) x / sum(x) )


viral_physeq %>% otu_table()
#Transform distance to the sqrt()
#Unifrac is a metric but non-Euclidean, Bray-Curits is a semi-metric, jsd is a semi-metric.
dist = sqrt(dist)
#example:dist = sqrt(phyloseq::distance(SmileT, "jsd"))

#ordination_betadiversity_PCOA
pcoa=ordinate(SmileT, "PCoA", distance=dist)
nmds=ordinate(SmileT,"NMDS",distance=dist)

#Permanova test using adonis function
adonis(dist ~ variable, as(sample_data(data), "data.frame"))

#Dispersion analysis
betatax=betadisper(dist,data.frame(sample_data(data))$variable)
p=permutest(betatax)
p$tab

#ANOSIM test
variable_group = get_variable(data, "variable")
variable_group = anosim(dist, variable_group)

----------------------------------------
  
  #Bray-Curtis Time decay (separate analysis for littoral and pelagic)
  
Date<-read.csv("",sep=";",row.names=1,header=T) #time matrix
Bray<-read.csv("",sep=";",row.names=1,header=T) #Bray-Curtis distance matrix

Bray2<-as.vector(matrix(unlist(Bray),ncol=1))
head(Bray2)

Date2<-as.vector(matrix(unlist(Date),ncol=1))
Bray2<-subset(Bray2,Date2!="NA")

Date2<-subset(Date2,Date2!="NA")

plot(Date2,Bray2)
plot(as.factor(Date2),Bray2)

k=365 #if years (time decay analysis)

Date3 <- k*(round(Date2/k,0))
Date3
cat<-sort(unique(Date3))
data<-sapply(cat,function(x){
  Â  
  Â  Brayx<-subset(Bray2,Date3==x)
  Meanx=mean(Brayx)
  Sex=std.error(Brayx)
  return(c(Meanx,Sex))
})