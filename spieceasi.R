#https://github.com/zdk123/SpiecEasi 
library(SpiecEasi)

cyano_ps

#ensure viral ps has same samples as cyano_ps 
meta2

colnames(asv_count) <- meta$description
vir_counts <- asv_count[,(colnames(asv_count) %in% meta2$description)]
rownames(vir_counts) <-paste0("vir_", rownames(vir_counts))
mock_taxa2<- mock_taxa
rownames(mock_taxa2)<- paste0("vir_", rownames(mock_taxa2))

#add ASV count table, metadata, virTree to phyloseq table
count_phy <- otu_table(vir_counts, taxa_are_rows=T)
sample_info <- sample_data(meta2)

#add to phyloseq object
virps <- phyloseq(count_phy, sample_info, mock_taxa2)
colnames(tax_table(virps)) <- c("Kingdom", "Phylum", "Class",
                                      "Order", "Family", "Genus", "ASV")

#reorder phyloseq by chronological date
otu_table(virps) <- otu_table(virps)[,toorder]
virps_filt <- filter_taxa(virps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)
cyanops_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)

virps_filt
cyanops_filt 

#spiec easi
spie <- spiec.easi(list(virps_filt, cyanops_filt), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep("Viral",ntaxa(virps_filt)), rep("Cyanobacteria",ntaxa(cyanops_filt)))
#plot(adj2igraph(getRefit(spie)), vertex.color=dtype+1, vertex.size=9)


#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(spie)))
dim(betaMatsym)

#check positive, negative, and total edges
(p.edge =length(betaMat[betaMatsym>0])/2)
(n.edge =length(betaMat[betaMatsym<0])/2)
(tot.edge =length(betaMat[betaMatsym!=0])/2)



#get weights
bm <- symBeta(getOptBeta(spie), mode="maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
FG.ig <- adj2igraph(Matrix::drop0(getRefit(spie)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(cyanops_filt))))

#plot with weights
plot_network(FG.ig, list(virps_filt, cyanops_filt))

#postive weights only:
weights.pos <- (1-Matrix::summary(t(bm))[,3])/2
FG.ig.pos <- adj2igraph(Matrix::drop0(getRefit(spie)),
                    edge.attr=list(weight=weights.pos))
plot_network(FG.ig.pos, list(virps_filt, cyanops_filt))

write.graph(FG.ig,"spieceasi.ncol.txt",format="ncol") 
head(corr.tab <- read.table("spieceasi.ncol.txt"))

#Create a custom color scale
library(RColorBrewer)
dtype <- as.factor(c(rep(1,ntaxa(virps_filt)), rep(2,ntaxa(cyanops_filt))))
otu.id <- c(taxa_names(virps_filt), taxa_names(cyanops_filt))
# colours.df <- cbind(data.frame(otu.id), data.frame(dtype))
# 
# myColors <- brewer.pal(2,"Set1")
# names(myColors) <- levels(colours.df$dtype)
# nodeColr <- scale_colour_manual(name = "dtype",values = myColors)



#https://ramellose.github.io/networktutorials/workshop_MDA.html
library(igraph)
#Network centrality: degree centrality (ie. degree = number of connections a node has)
spiec.deg <- degree(FG.ig)
hist(spiec.deg)
range(spiec.deg)

#if the degree distribution of a network follows a power law, that network is scale-free
plaw.fit <- fit_power_law(spiec.deg) #The fit_power_law functions fits a power law to the degree distribution of the network.
plaw.fit
#The values for the fit are compared to expected values with the Kolmogorov-Smirnov test. 
#The null hypothesis for this test is that the degree distribution is drawn from a reference distribution. 
#In this case, that reference distribution is generated from a power law.

#The null hypothesis can only be rejected if the p-value of the test is below 0.05. 
#Here, the p-value is 1. Therefore, we cannot conclude that the degree distribution is drawn from a different distribution than the power-law distribution.
#Scale-free networks are networks with a degree distribution that follows a power law. 
#Our result indicates that the network may be scale-free and contains nodes with a degree far larger than the average degree. 
#While there is not that much known about the effect of scale-freeness on microbial networks, studies (e.g., Cohen et al 2001, https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.3682) 
#indicate that scale-freeness decreases the network’s sensitivity to random attacks. 
#However, we still do not know to what extent biological networks follow a power law as we have few true biological networks.
#Lima-Mendez and van Helden (2009) (https://pubs.rsc.org/en/content/articlehtml/2009/mb/b908681a) discuss some of the weaknesses of this theory.

library(ggnet)
ggnet2(FG.ig.pos,
       color = dtype, palette = c("Viral" = "#E1AF00", "Cyanobacteria" = "steelblue"), 
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = spiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 7,
       edge.size = weights.pos, edge.alpha = 0.25,
       label = otu.id, label.size = 1)+
  ggtitle("Viral and Cyanobacteria correlation network")
 # guides(color=FALSE)
  


