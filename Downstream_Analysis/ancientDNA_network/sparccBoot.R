args <- commandArgs(TRUE)
suppressMessages(library(ggplot2))
suppressMessages(library(igraph))
suppressMessages(library(qgraph))
suppressMessages(library(vegan))
suppressMessages(library(MCL))
suppressMessages(library(SpiecEasi))
suppressMessages(library(stringr))

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
sparcc.net.bn <- betweenness(sparcc.net)
sparcc.net.pr <- page_rank(sparcc.net)$vector
sparcc.net.hs <- hub_score(sparcc.net)$vector

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
as.data.frame(sparcc.net.bn.top5) -> hubs.bn
hubs.bn$pop <- rep(outName, 5)

sparcc.net.tr.global <- transitivity(sparcc.net, type = "global")
write.table(sparcc.net.tr.global, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_transitivity.txt"), quote =F,  col.names = F, row.names = T)

write.table(hubs.hs, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_5hubs.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.cn, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_5closeness.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.bn, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_5between.txt"), quote =F,  col.names = F, row.names = T)
write.table(hubs.pr, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_5pageRank.txt"), quote =F,  col.names = F, row.names = T)

#Make clusters - i don't really understand the details here

# Get clusters
sparcc.wt <- walktrap.community(sparcc.net)
sparcc.ml <- multilevel.community(sparcc.net)
# Get membership of walktrap clusters
membership(sparcc.wt) -> member
#print(member)
max(member) -> numClusters
write.table(numClusters, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_numClusters.txt"), quote =F,  col.names = F, row.names = T)

# Get clusters using MCL method
sparcc.adj <- as_adjacency_matrix(sparcc.net)
sparcc.mc <- mcl(sparcc.adj, addLoops = TRUE)

#Basic features of the network

# Network features
sparcc.nodes <- V(sparcc.net)
sparcc.edges <- V(sparcc.net)
sparcc.node.names <- V(sparcc.net)$name
sparcc.num.nodes <- vcount(sparcc.net)
sparcc.num.edges <- ecount(sparcc.net)

sparcc.AP <- articulation.points(sparcc.net)
as_ids(sparcc.AP) -> AP_ids
write.table(AP_ids, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_AP.txt"), quote =F,  col.names = F, row.names = T)
modularity(sparcc.net, membership(sparcc.wt)) -> mod
write.table(mod, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_modularity.txt"), quote =F,  col.names = F, row.names = T)


matches <- data.frame(nrows = numClusters, ncols = 2)
for(i in 1:numClusters){
    as.character(i) -> fpl
    matches[i,1] <- i
    which(member == fpl) -> l
    length(l) -> matches[i,2]
    #sum(str_count(member, fpl)) -> matches[i,2]
}
colnames(matches) <- c("cluster", "numberMatches")
write.table(matches, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_cluster_Matches.txt"), quote =F,  col.names = T, row.names = F)
as.data.frame(as.character(member)) -> m
names(member) -> m$taxa
colnames(m) <- c("cluster", "taxa")
write.table(m, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_membership.txt"), quote =F,  col.names = T, row.names = F)

as.data.frame(summary(m$cluster)) -> temp1
colnames(temp1) <- c("clusterSize")
sort(temp1$clusterSize, decreasing = T) -> temp2
print(grep("TRUE",cumsum(temp2) > 0.5 * sum(temp2))[1]) -> fifty
print(grep("TRUE",cumsum(temp2) > 0.75 * sum(temp2))[1]) -> sevenfive
print(grep("TRUE",cumsum(temp2) > 0.90 * sum(temp2))[1]) -> ninety

write.table(fifty, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_50.txt"), quote =F,  col.names = T, row.names = F)
write.table(sevenfive, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_75.txt"), quote =F,  col.names = T, row.names = F)
write.table(ninety, file = paste0("/Users/dave/Desktop/microbiom_R_network/boot/", outName, "_90.txt"), quote =F,  col.names = T, row.names = F)
