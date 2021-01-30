#https://github.com/zdk123/SpiecEasi 
library(SpiecEasi)
library(igraph)

cyano_ps <- subset_taxa(bact_physeq, Phylum == "p__Cyanobacteria")
doli_ps <- subset_taxa(bact_physeq, Genus == "g__Dolichospermum")
micro_ps <- subset_taxa(bact_physeq, Genus == "g__Microcystis")

#replace name so don't have to edit whole script
cyano_ps <- micro_ps

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

#taxa_names(viral_physeq) <- paste0("vir_", taxa_names(viral_physeq))
taxa_names(doli_ps) <- paste0("doli_", taxa_names(doli_ps))
taxa_names(micro_ps) <- paste0("micro_", taxa_names(micro_ps))

# colnames(tax_table(viral_physeq)) <- c("Kingdom", "Phylum", "Class",
#                                        "Order", "Family", "Genus", "ASV")

doli.ps <- prune_samples(rownames(sample_data(doli_ps)) %in% rownames(meta2), doli_ps)
micro.ps <- prune_samples(rownames(sample_data(micro_ps)) %in% rownames(meta2), micro_ps)

#reorder phyloseq by chronological date
otu_table(virps) <- otu_table(virps)[,toorder]
otu_table(doli.ps) <- otu_table(doli.ps)[,toorder]
otu_table(micro.ps) <- otu_table(micro.ps)[,toorder]


virps_filt <- filter_taxa(virps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)

cyanops_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)
# doli_filt <- filter_taxa(doli_ps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)
# micro_filt <- filter_taxa(micro_ps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)


virps_filt
cyanops_filt 

#spiec easi
spie <- spiec.easi(list(virps_filt, doli.ps, micro.ps), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))


#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(spie)))
dim(betaMatsym)

#select for cyano - viral connections only
list(name=c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps)))
length(taxa_names(virps_filt))
length(c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps)))

# vir.cyan <- betaMatsym[1:576, 577:614]
# 
# #check positive, negative, and total edges (divide by 2 because an edge is represented by 2 entries in the matrix)
# (p.edge =length(betaMatsym[betaMatsym>0])/2)
# (n.edge =length(betaMatsym[betaMatsym<0])/2)
# (tot.edge =length(betaMatsym[betaMatsym!=0])/2)
# 
# #viral-cyano connections only
# (p.edge =length(vir.cyan[vir.cyan>0]))
# (n.edge =length(vir.cyan[vir.cyan<0]))
# (tot.edge =length(vir.cyan[vir.cyan!=0]))


#get weights
bm <- symBeta(getOptBeta(spie), mode="maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
FG.ig <- adj2igraph(Matrix::drop0(getRefit(spie)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(doli.ps), taxa_names(micro.ps))))

#plot with weights
#plot_network(FG.ig, list(virps_filt, cyanops_filt))

#postive weights only:
weights.pos <- (1-Matrix::summary(t(bm))[,3])/2
FG.ig.pos <- adj2igraph(Matrix::drop0(getRefit(spie)),
                    edge.attr=list(weight=weights.pos))
#plot_network(FG.ig.pos, list(virps_filt, cyanops_filt))

write.graph(FG.ig,"spieceasi.ncol.txt",format="ncol") 
head(corr.tab <- read.table("spieceasi.ncol.txt"))


#isolate for viral-cyano interactions only
vircyn.pos <- corr.tab %>% 
  filter(across(V2, ~ !grepl('vir_', .))) %>%
  filter(across(V1, ~grepl('vir_', .))) %>%
  rename(weight = V3) %>%
  filter(weight > 0)
head(vircyn.pos)

#plot vircyn connections with weights only
vircyn.plot <- graph_from_data_frame(vircyn.pos, directed = TRUE, vertices = NULL)
plot_network(vircyn.plot)

which(as.data.frame(str_count(vircyn.pos$V2, "micro_"))=="1", arr.ind=T) #see which positions micro_ are in in col 2 of df

repdm <- list()
for (j in str_count(vircyn.pos$V2, "doli_")){
  if (j == 1){
    repdm[length(repdm)+1] <- print("Dolichospermum")
  } else {
    repdm[length(repdm)+1] <- print("Microcystis")
  }
}
repdm <- unlist(repdm)
head(repdm, n=8)

dtype <- as.factor(c(rep("Phage", length(unique(vircyn.pos[,1]))), repdm))
otu.id <- c(as.character(vircyn.pos[,1]), as.character(vircyn.pos[,2]))


#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
spiec.deg <- igraph::degree(vircyn.plot)
hist(spiec.deg)
range(spiec.deg)


#if the degree distribution of a network follows a power law, that network is scale-free
plaw.fit <- fit_power_law(spiec.deg) #The fit_power_law functions fits a power law to the degree distribution of the network.
plaw.fit
#The values for the fit are compared to expected values with the Kolmogorov-Smirnov test. 
#The null hypothesis for this test is that the degree distribution is drawn from a reference distribution. 
#In this case, that reference distribution is generated from a power law.

#The null hypothesis can only be rejected if the p-value of the test is below 0.05. 
#Here, the p-value is 0.321. Therefore, we cannot conclude that the degree distribution is drawn from a different distribution than the power-law distribution.
#Scale-free networks are networks with a degree distribution that follows a power law. 
#Our result indicates that the network may be scale-free and contains nodes with a degree far larger than the average degree. 
#While there is not that much known about the effect of scale-freeness on microbial networks, studies (e.g., Cohen et al 2001, https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.86.3682) 
#indicate that scale-freeness decreases the network’s sensitivity to random attacks. 
#However, we still do not know to what extent biological networks follow a power law as we have few true biological networks.
#Lima-Mendez and van Helden (2009) (https://pubs.rsc.org/en/content/articlehtml/2009/mb/b908681a) discuss some of the weaknesses of this theory.

library(ggnet)
ggnet2(vircyn.plot,
       color = dtype, palette = c("Phage" = "#E1AF00", "Dolichspermum" = "red", "Microcystis" = "steelblue"), 
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = spiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 7,
       edge.size = vircyn.pos[,3], edge.alpha = 0.5,
       label = otu.id, label.size = 1)+
  ggtitle("Viral and Microcystis correlation network")
 # guides(color=FALSE)

  
#Check which OTUs are part of different modules.
#https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html

#GREEDY COMMUNITY DETECTION
clusters<-cluster_fast_greedy(as.undirected(vircyn.plot))
clusters
modularity(clusters)
#modularity matrix
B = modularity_matrix(vircyn.plot, membership(clusters))
round(B[1,],5)
#membership of nodes
membership(clusters)
#number of communities
length(clusters)
#size of communities
sizes(clusters)
#crossing edges
crossing(clusters, vircyn.plot)

#plot communities without shaded regions
ggnet2(vircyn.plot,
       color = membership(clusters),
       alpha=0.75,
       shape = factor(dtype),
       shape.legend = "Type",
       node.size = spiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 7,
       edge.size = vircyn.pos[,3], edge.alpha = 0.5,
       label = otu.id, label.size = 1)+
  ggtitle("Viral and Microcystis correlation network by clusters")
# guides(color=FALSE)

#plot dendogram
plot_dendrogram(clusters)

#Check which OTUs are part of different modules.
clustersOneIndices=which(clusters$membership==1)
clustersOneOtus=clusters$names[clustersOneIndices]
clustersOneOtus
#OR
clusters[1]
