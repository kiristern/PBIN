#https://fukamilab.github.io/BIO202/09-B-networks.html

library(igraph)
library(bipartite)

vircyn.plot2 <- vircyn.plot

web <- as_adjacency_matrix(vircyn.plot2, attr = "weight")
#bipartite::visweb(web)

web2 <- get.adjacency(vircyn.plot2, sparse=F, attr = "weight")

plotweb(sortweb(web2, sort.order = "inc"), method="normal", text.rot=90, col.low="#1F78B4", col.high="#B2DF8A", col.interaction="#A6CEE3", labsize = 0.4)



