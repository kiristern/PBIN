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

#reorder phyloseq by chronological date
otu_table(virps) <- otu_table(virps)[,toorder]
virps_filt <- filter_taxa(virps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)
cyanops_filt <- filter_taxa(cyano_ps, function(x) sum(x > 1e-5) > (0.10*length(x)), TRUE)

virps_filt
cyanops_filt 

#spiec easi
spie <- spiec.easi(list(virps_filt, cyanops_filt), method='mb', nlambda=40,
                      lambda.min.ratio=1e-2, pulsar.params = list(thresh = 0.05))

dtype <- c(rep(1,ntaxa(virps_filt)), rep(2,ntaxa(cyanops_filt)))
plot(adj2igraph(getRefit(spie)), vertex.color=dtype+1, vertex.size=9)



#http://psbweb05.psb.ugent.be/conet/microbialnetworks/spieceasi.php
betaMatsym <- as.matrix(symBeta(getOptBeta(spie)))
dim(betaMatsym)
betaMat <- getOptBeta(spie)
dim(betaMat)

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
                    edge.attr=list(weight=weights.pos),
                    vertex.attr = list(name=c(taxa_names(virps_filt), taxa_names(cyanops_filt))))
plot_network(FG.ig.pos, list(virps_filt, cyanops_filt))

write.graph(FG.ig,"spieceasi.ncol.txt",format="ncol") 
head(spea2cytoscape <- read.table("spieceasi.ncol.txt"))

