#https://github.com/zdk123/SpiecEasi 
library(SpiecEasi)
library(igraph)

viral_physeq2 <- virps3000

#add to phyloseq object
colnames(tax_table(viral_physeq2)) <- c("Kingdom", "Phylum", "Class",
                                "Order", "Family", "Genus", "ASV")


virps_filt <- filter_taxa(viral_physeq2, function(x) sum(x > 1) > (0.10*length(x)), TRUE)

virps_filt

#spiec easi
vir.spie2 <- spiec.easi(virps_filt, method='mb', nlambda=125,
                   lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05,
                                                               subsample.ratio=0.8,
                                                               seed = 1234,
                                                               ncores=4))

vir.spie2$select$stars$summary #if coming up with empty network: b/c max value of the StARS summary statistic never crosses the default threshold (0.05). fix by lowering lambda.min.ratio to explore denser networks
getStability(vir.spie2)
sum(getRefit(vir.spie2))/2

#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(vir.spie2)))

#get weights
bm2 <- symBeta(getOptBeta(vir.spie2), mode="maxabs")
diag(bm2) <- 0
weights2 <- Matrix::summary(t(bm2))[,3]
FG.ig.vir <- adj2igraph(Matrix::drop0(getRefit(vir.spie2)),
                     edge.attr=list(weight=weights2),
                     vertex.attr=list(name=taxa_names(virps_filt)))

#plot with weights
#plot_network(FG.ig, list(virps_filt, cyanops_filt))

vir.corr.tab <- igraph::as_data_frame(FG.ig.vir, what="edges")
write.csv(vir.corr.tab, "vir-vir.covar.csv")

#plot vircyn connections with weights only
virplot <- graph_from_data_frame(vir.corr.tab, directed = TRUE, vertices = NULL)


#https://ramellose.github.io/networktutorials/workshop_MDA.html
#Network centrality: degree centrality (ie. degree = number of connections a node has)
virspiec.deg <- igraph::degree(virplot)
hist(virspiec.deg)
range(virspiec.deg)


library(ggplot2)
library(ggnet)
ggnet2(virplot,
       alpha=0.75,
       #shape = factor(dtype),
       #shape.legend = "Type",
       node.size = virspiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 4,
       edge.size = abs(vir.corr.tab[,3]), edge.alpha = 0.5, edge.lty = ifelse(vir.corr.tab$weight > 0, 1, 2),
       label = colnames(vir.spie2$est$data), label.size = 1)+
  ggtitle("Viral correlation network")
# guides(color=FALSE)



#Check which OTUs are part of different modules.
#https://users.dimi.uniud.it/~massimo.franceschet/R/communities.html

#GREEDY COMMUNITY DETECTION
clusters.vir<- cluster_fast_greedy(as.undirected(virplot), weights = abs(E(virplot)$weight))

modularity(clusters.vir)

#membership of nodes
membership(clusters.vir)
#number of communities
length(clusters.vir)
#size of communities
sizes(clusters.vir)
#crossing edges
crossing(clusters.vir, virplot)

#see which edge connects two different communities
which(crossing(clusters.vir, virplot) == T)
length(which(crossing(clusters.vir, virplot) == T)) #number of cross community interactions


#plot communities without shaded regions

ggnet2(virplot,
       color = membership(clusters.vir),
       alpha=0.75,
       node.size = virspiec.deg,
       size.legend = "Degree of Centrality",
       size.cut = 3,
       edge.size = abs(vir.corr.tab[,3]), edge.alpha = 0.5, edge.lty = ifelse(vir.corr.tab$weight > 0, 1, 2),
       label = colnames(vir.spie2$est$data), label.size = 1)+
  ggtitle("Viral correlation network by clusters")
#guides(size=FALSE)




