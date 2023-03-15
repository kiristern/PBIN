#https://fukamilab.github.io/BIO202/09-B-networks.html

library(igraph)
library(bipartite)

vircyn.plot2 <- vircyan.plot

web <- as_adjacency_matrix(vircyn.plot2, attr = "weight")
#bipartite::visweb(web)

web2 <- get.adjacency(vircyn.plot2, sparse=F, attr = "weight")

plotweb(sortweb(web2, sort.order = "inc"), method="normal", text.rot=90, col.low="#1F78B4", col.high="#B2DF8A", col.interaction="#A6CEE3", labsize = 0.4)
plotweb(sortweb(web2, sort.order = "inc"), method="normal", text.rot=90, col.low="#1F78B4", col.high="#B2DF8A", col.interaction="#A6CEE3", labsize = 0.3,
        bor.col.high="darkgrey", bor.col.low="darkgrey")

visweb(web)





#some try didnot work
require(ggplot2)
require(network)
require(igraph)
require(sna)
require(GGally)
require(ergm)
require(intergraph)
require(RColorBrewer)

bip_railway <- function (mymat, nodesize=9, label=F) {
  require(GGally)
  # Coords for mode "A"
  coordP<- cbind(rep(2,dim(mymat)[1]), seq(1, dim(mymat)[1])+2)
  # Coords for mode "P"
  coordA<- cbind(rep(4,dim(mymat)[2]), seq(1, dim(mymat)[2])+2)
  mylayout<- as.matrix(rbind(coordP, coordA))
  #
  # Initialize and plot the network with a railway layout.
  test.net<- bip_init_network(mymat)
  p<- GGally::ggnet2(test.net, mode=mylayout, label=label,
                     size= nodesize, label.size=nodesize/3,
                     layout.exp=1.5) +
    coord_flip()
  p
}
bip_init_network <- function (mymat, mode1="P", mode2="A") {
  require(network)
  require(GGally)
  if(!is.matrix(mymat)) mymat <- as.matrix(mymat)
  p<- dim(mymat)[1]    # Plants are rows
  a<- dim(mymat)[2]    # Animals are columns
  net<- network::network(mymat,
                         matrix.type = "bipartite",
                         ignore.eval = FALSE,
                         names.eval = "weights")
  net
  network::set.vertex.attribute(net,"mode",c(rep(mode1,p), rep(mode2,a)))
}

bip_railway(sortweb(web2, sort.order = "inc")) 
g<- bip_railway(web2, label=T)
g

