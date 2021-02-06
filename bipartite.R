#https://fukamilab.github.io/BIO202/09-B-networks.html

library(igraph)
library(bipartite)

vircyn.plot2 <- vircyn.plot

web <- as_adjacency_matrix(vircyn.plot2, attr = "weight")
bipartite::visweb(web)

web2 <- get.adjacency(vircyn.plot2, sparse=F, attr = "weight")

plotweb(sortweb(web2, sort.order = "inc"), method="normal", text.rot=90, col.low="green", col.high="blue")




# #https://rpubs.com/pjmurphy/317838
# V(vircyn.plot2)$type <- V(vircyn.plot2)$name %in% vircyn.pos[,1]
# 
# V(vircyn.plot2)$shape <- ifelse(V(vircyn.plot2)$type, "circle", "square")
# 
# V(vircyn.plot2)$color <- ifelse(V(vircyn.plot2)$type, "lightblue", "salmon")
# V(vircyn.plot2)$shape <- ifelse(V(vircyn.plot2)$type, "circle", "square")
# E(vircyn.plot2)$color <- "lightgrey"
# V(vircyn.plot2)$label.color <- "black"
# V(vircyn.plot2)$label.cex <- 0.5
# V(vircyn.plot2)$frame.color <- "grey"
# V(vircyn.plot2)$size <- 5
#  
# plot(as.undirected(vircyn.plot2), 
#      layout=layout.bipartite, 
#      vertex.size=3, 
#      vertex.label.cex=0.5,
#      asp=0.2) #asp controls how rectangular the plot is. <1 = wide, >1 = tall
# 
# coords <- layout_as_bipartite(vircyn.plot2)
# head(coords) #This matrix shows the x coordinates of your vertices in the 1st col and the y coordinates in the 2nd col, ordered according to your list with names.
# 
# #graph df
# vircyn.pos




  



