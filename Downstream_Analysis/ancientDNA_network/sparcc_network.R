args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(qgraph))
suppressMessages(library(vegan))
suppressMessages(library(MCL))
suppressMessages(library(SpiecEasi))

read.table(args[1], header =T, row.names =1, sep = "\t") -> inFile
t(inFile) -> tFile
tFile * 1000 -> forNet
args[2] -> outName


sparcc.matrix <- sparcc(forNet)
sparcc.cutoff <- 0.3
sparcc.adj <- ifelse(abs(sparcc.matrix$Cor) >= sparcc.cutoff, 1, 0)
# Add OTU names to rows and columns
rownames(sparcc.adj) <- colnames(forNet)
colnames(sparcc.adj) <- colnames(forNet)
# Build network from adjacency
sparcc.net <- graph.adjacency(sparcc.adj,
                              mode = "undirected",
                              diag = FALSE)

#To find hubs (large number of links compared to other nodes) - can be thought of as keystone species

# Use sparcc.net for the rest of the method
#sparcc.net
# Hub detection
sparcc.net.cn <- closeness(sparcc.net)
sparcc.net.bn <- betweenness(sparcc.net, normalized =T)
sparcc.net.pr <- page_rank(sparcc.net)$vector
sparcc.net.hs <- hub_score(sparcc.net)$vector
sparcc.net.tr <- transitivity(sparcc.net, type = "local")
sparcc.net.tr.global <- transitivity(sparcc.net, type = "global")

node.names <- V(sparcc.net)$name
#print(node.names)
#print(sparcc.net.tr)

#print(sparcc.net.bn)

#The net.hs is a vector of OTUs with value from 0-1, with values closer to 1 being the more likley hubs

# Sort the species based on hubbiness score
sparcc.net.hs.sort <- sort(sparcc.net.hs, decreasing = TRUE)
sparcc.net.cn.sort <- sort(sparcc.net.cn, decreasing = TRUE)
sparcc.net.bn.sort <- sort(sparcc.net.bn, decreasing = TRUE)
sparcc.net.pr.sort <- sort(sparcc.net.pr, decreasing = TRUE)
# Choose the top 5 keystone species
sparcc.net.hs.top5 <- head(sparcc.net.hs.sort, n = 5)
sparcc.net.cn.top5 <- head(sparcc.net.cn.sort, n = 5)
sparcc.net.bn.top5 <- head(sparcc.net.bn.sort, n = 5)
sparcc.net.pr.top5 <- head(sparcc.net.pr.sort, n = 5)

as.data.frame(sparcc.net.hs.top5) -> hubs.hs
hubs.hs$pop <- rep(outName, 5)
as.data.frame(sparcc.net.cn.top5) -> hubs.cn
hubs.cn$pop <- rep(outName, 5)
as.data.frame(sparcc.net.pr.top5) -> hubs.pr
hubs.pr$pop <- rep(outName, 5)
as.data.frame(sparcc.net.bn.sort) -> hubs.bn
hubs.bn$pop <- rep(outName, length(as.data.frame(sparcc.net.bn.sort)))
as.data.frame(sparcc.net.tr) -> hubs.tr
hubs.tr$taxa <- node.names

net.knn <- knn(sparcc.net, vids = V(sparcc.net))
#print(net.knn$knn)
print(sparcc.net.tr.global)
write.table(sparcc.net.tr.global, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_globalTrans.txt"), quote =F,  col.names = F, row.names = T)


write.table(hubs.hs, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_5hubs.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.cn, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_5closeness.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.bn, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_allBetween.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.pr, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_5pageRank.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.tr, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_allTrans.txt"), quote =F,  col.names = F, row.names = T)
#write.table(sparcc.net.tr.global, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_globalTrans"), quote =F,  col.names = F, row.names = T)

#Make clusters - i don't really understand the details here

# Get clusters
sparcc.wt <- walktrap.community(sparcc.net)
sparcc.ml <- multilevel.community(sparcc.net)
# Get membership of walktrap clusters
membership(sparcc.wt) -> member
modularity(sparcc.net, membership(sparcc.wt)) -> mod
write.table(mod, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_modularity.txt"), quote =F,  col.names = F, row.names = T)
#print(member)
#print(max(member))
# Get clusters using MCL method
sparcc.adj <- as_adjacency_matrix(sparcc.net)
sparcc.mc <- mcl(sparcc.adj, addLoops = TRUE)

connectVec <- vector(mode = "numeric", length = length(V(sparcc.net)$name))
nodeNames <- V(sparcc.net)$name
for(i in 1:length(V(sparcc.net)$name)){
    connectVec[i] <- length(neighbors(sparcc.net, nodeNames[i]))
}
as.data.frame(connectVec) -> connectDF
connectDF$taxa <- nodeNames
zero <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec <= 0)))]
one5 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 1 & connectDF$connectVec <= 5)))]
six10 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 6 & connectDF$connectVec <= 10)))]
eleven15 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 11 & connectDF$connectVec <= 15)))]
sixteen20 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 16 & connectDF$connectVec <= 20)))]
twentyone25 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 21 & connectDF$connectVec <= 25)))]
twentysix30 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 26 & connectDF$connectVec <= 30)))]
thirtyone35 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 31 & connectDF$connectVec <= 35)))]
thirtysix40 <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 36 & connectDF$connectVec <= 40)))]
fortyAbove <- nodeNames[as.numeric(row.names(subset(connectDF, connectDF$connectVec >= 41)))]
as.data.frame(zero) -> zeroDF
as.data.frame(one5) -> one5DF
as.data.frame(six10) -> six10DF
as.data.frame(eleven15) -> eleven15DF
as.data.frame(sixteen20) -> sixteen20DF
as.data.frame(twentyone25) -> twentyone25DF
as.data.frame(twentysix30) -> twentysix30DF
as.data.frame(thirtyone35) -> thirtyone35DF
as.data.frame(thirtysix40) -> thirtysix40DF
as.data.frame(fortyAbove) -> fortyAboveDF
zeroDF$pop <- rep(outName, length(rownames(zeroDF)))
one5DF$pop <- rep(outName, length(rownames(one5DF)))
six10DF$pop <- rep(outName, length(rownames(six10DF)))
eleven15DF$pop <- rep(outName, length(rownames(eleven15DF)))
sixteen20DF$pop <- rep(outName, length(rownames(sixteen20DF)))
twentyone25DF$pop <- rep(outName, length(rownames(twentyone25DF)))
twentysix30DF$pop <- rep(outName, length(rownames(twentysix30DF)))
thirtyone35DF$pop <- rep(outName, length(rownames(thirtyone35DF)))
thirtysix40DF$pop <- rep(outName, length(rownames(thirtysix40DF)))
fortyAboveDF$pop <- rep(outName, length(rownames(fortyAboveDF)))

write.table(zeroDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_zero.txt"), quote =F,  col.names = F, row.names = F)
write.table(one5DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_one5.txt"), quote =F,  col.names = F, row.names = F)
write.table(six10DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_six10.txt"), quote =F,  col.names = F, row.names = F)
write.table(eleven15DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_eleven15.txt"), quote =F,  col.names = F, row.names = F)
write.table(sixteen20DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_sixteen20.txt"), quote =F,  col.names = F, row.names = F)
write.table(twentyone25DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_twentyone25.txt"), quote =F,  col.names = F, row.names = F)
write.table(twentysix30DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_twentysix30.txt"), quote =F,  col.names = F, row.names = F)
write.table(thirtyone35DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_thirtyone35.txt"), quote =F,  col.names = F, row.names = F)
write.table(thirtysix40DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_thirtysix40.txt"), quote =F,  col.names = F, row.names = F)
write.table(fortyAboveDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_fortyAbove.txt"), quote =F,  col.names = F, row.names = F)

write.table(connectDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/connections_", outName, "_rawConnnects.txt"), quote =F,  col.names = F, row.names = F)

names(member[grep("^1$", member)]) -> cluster1
names(member[grep("^2$", member)]) -> cluster2
names(member[grep("^3$", member)]) -> cluster3
names(member[grep("^4$", member)]) -> cluster4
names(member[grep("^5$", member)]) -> cluster5
names(member[grep("^6$", member)]) -> cluster6
names(member[grep("^7$", member)]) -> cluster7

write.table(cluster1, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_cluster1.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster2, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_cluster2.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster3, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_cluster3.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster4, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_cluster4.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster5, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc", outName, "_cluster5.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster6, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc", outName, "_cluster6.txt"), quote =F,  col.names = F, row.names = T)
write.table(cluster7, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc", outName, "_cluster7.txt"), quote =F,  col.names = F, row.names = T)


#Basic features of the network

# Network features
sparcc.nodes <- V(sparcc.net)
sparcc.edges <- V(sparcc.net)
sparcc.node.names <- V(sparcc.net)$name
sparcc.num.nodes <- vcount(sparcc.net)
sparcc.num.edges <- ecount(sparcc.net)


# Function 2: Plot network with node size scaled to hubbiness
plot.net <- function(net, scores, outfile, title) {
  # Convert node label from names to numerical IDs.
  features <- V(net)$name
  col_ids <- seq(1, length(features))
  V(net)$name <- col_ids
  node.names <- features[V(net)]
  # Nodes' color.
  V(net)$color <- "white"
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex = 1)
  title(title, cex.main = 4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  legend("bottom", legend = labels, xpd = TRUE, ncol = 3, cex = 0.8)
  dev.off()
}

plot.net(sparcc.net, sparcc.net.hs, outfile = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_sparcc_uncolored"), title = outName)

sparcc.AP <- articulation.points(sparcc.net)
#as.data.frame(as.character(sparcc.AP)) -> sparcc.AP.DF
print(sparcc.AP)
#write.table(sparcc.AP.DF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_AP.txt"), quote =F,  col.names = F, row.names = T)

outName_degree.dist <- degree_distribution(sparcc.net, mode = "all", cumulative = F)
as.data.frame(outName_degree.dist) -> outDF
write.table(outDF, file = paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_sparcc_degree.txt"), quote =F,  col.names = F, row.names = T)

# Function 3: Plot network with clusters and node size scaled to hubbiness
plot.net.cls <- function(net, scores, cls, AP, outfile, title) {
  # Get size of clusters to find isolated nodes.
  cls_sizes <- sapply(groups(cls), length)
  # Randomly choosing node colors. Users can provide their own vector of colors.
  colors <- sample(colours(), length(cls))
  # Nodes in clusters will be color coded. Isolated nodes will be white.
  V(net)$color <- sapply(membership(cls),
                         function(x) {ifelse(cls_sizes[x]>1,
                                             colors[x], "white")})
  # Convert node label from names to numerical IDs.
  node.names <- V(net)$name
  col_ids <- seq(1, length(node.names))
  V(net)$name <- col_ids
  # To draw a halo around articulation points.
  #AP <- lapply(names(AP), function(x) x)
  #marks <- lapply(1:length(AP), function(x) which(node.names == AP[[x]]))
  # Define output image file.
  outfile <- paste(outfile, "jpg", sep=".")
  # Image properties.
  jpeg(outfile, width = 4800, height = 9200, res = 300, quality = 100)
  par(oma = c(4, 1, 1, 1))
  # Customized layout to avoid nodes overlapping.
  e <- get.edgelist(net)
  class(e) <- "numeric"
  l <- qgraph.layout.fruchtermanreingold(e, vcount=vcount(net),
                                         area=8*(vcount(net)^2),
                                         repulse.rad=(vcount(net)^3.1))
  # Main plot function.
  plot(net, vertex.size = (scores*5)+4, vertex.label.cex=0.9,
       vertex.label.color = "black",
       mark.border="black",
       #mark.groups = marks,
       mark.col = "white",
       mark.expand = 10,
       mark.shape = 1,
       layout=l)
  title(title, cex.main=4)
  # Plot legend containing OTU names.
  labels = paste(as.character(V(net)), node.names, sep = ") ")
  legend("bottom", legend = labels, xpd = TRUE, ncol = 3, cex = 0.8)
  dev.off()
}
# Execute this command after running Function 3
plot.net.cls(sparcc.net, sparcc.net.hs, sparcc.wt, sparcc.AP,
             outfile =  paste0("/Users/dave/Desktop/microbiom_R_network/ancient_calc/", outName, "_sparcc_colored"), title = outName)
