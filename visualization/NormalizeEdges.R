#########################################################################################################################
#########################################################################################################################
###### PROJECT:        DBNs
###### NAME:           NormalizeEdges.R
###### AUTHOR:         Daniel Ruiz-Perez, PhD Student
###### AFFILIATION:    Florida International University
###### 
###### DESCRIPTION:    This file normalizes the weight of the files in a folder by updating them
#########################################################################################################################
#########################################################################################################################



library(stringr)
library(igraph)
options("scipen"=100, "digits"=4)

folder="/Users/elyonyi/Documents/GitHub/longitudinal_microbiome_analysis_public/visualization/"
setwd(folder)

files = list.files(path = folder, pattern = ".graphml")

for (f in files){
  
  nameOfNetwork = paste(folder,f,sep="/")
  abundances = unlist(read.csv('meanAbundanceInfantFilteredAlignment.txt',sep="\t"))
  
  networkO = readChar(nameOfNetwork, file.info(nameOfNetwork)$size)
  network= unlist(strsplit(x =networkO, "\r\n"))
  network= unlist(strsplit(x =network, "\n"))
  
  ################## EDGE WIDTH
  weights = str_extract(network[grep("key=\"key_weight",network)], "\\-*\\d+\\.+\\d+")
  weights = weights[!is.na(weights)]
  
  ########### WE NORMALIZE THE WEIGHTS BASED ON w_i \average{x_i} / \sum_j  abs(w_j) \average{x_j}
  
  library(igraph)
  edgeList = igraph::get.edgelist(igraph::read.graph(nameOfNetwork, format = "graphml"))
  edgeAttr = igraph::get.edge.attribute(igraph::read.graph(nameOfNetwork, format = "graphml"))
  nodeAttr = igraph::get.vertex.attribute(igraph::read.graph(nameOfNetwork, format = "graphml"))
  
  #This function returns the abundance of a node based on the index of the incoming node
  abundanceOfTaxa <- function(nameOfIncomingNode,abundances){
    # we look for the abundance of that parent node
    abundance = abundances[which(names(abundances) %in% substr(nameOfIncomingNode,1,nchar(nameOfIncomingNode)-3))] 
    if (length(abundance)==0){
      #the node doesn't have abundance, it is clinical. We make it to be an average of the other abundances
      abundance = mean(abundances)
    }
    return (abundance)
  }
  
  normalizedWeights =c()
  #Iterate over all posible children. to normalize incoming edges to each of them
  for (i in 1:length(nodeAttr[[1]])){
    indexOfIncomingEdges = which(edgeList[,2]==i)
    indexOfOutgoingEdges = edgeList[indexOfIncomingEdges,1]
    
    #We need to calculate the sum of the multiplication of each weight times the mean abundance for all incoming edges of this children
    sum = 0
    for (j in 1:length(indexOfIncomingEdges)){
      nameOfIncomingNode = nodeAttr$id[indexOfOutgoingEdges[j]]
      abundance = abundanceOfTaxa(nameOfIncomingNode,abundances)
      sum = sum + abs(edgeAttr$weight[indexOfIncomingEdges[j]]*abundance)
    }
    
    #Iterate over all incoming edges for this children
    for (j in 1:length(indexOfIncomingEdges)){
      nameOfIncomingNode = nodeAttr$id[indexOfOutgoingEdges[j]]
      abundance = abundanceOfTaxa(nameOfIncomingNode,abundances)
      normalizedWeights = as.vector(c(normalizedWeights,edgeAttr$weight[indexOfIncomingEdges[j]]*abundance/sum))
    }
  }
  
  #we update the weights in the file
  for(i in 1:length(weights)){
    networkO = gsub(weights[i], normalizedWeights[i], networkO)
  }
  
  #save it and take care of the null character that appears
  writeChar(networkO,nameOfNetwork,nchars = nchar(networkO))
  r = readBin(nameOfNetwork, raw(), file.info(nameOfNetwork)$size)
  r[r==as.raw(0)] = as.raw(0x20) ## replace with 0x20 = <space>
  writeBin(r, nameOfNetwork)

  # if (f %in% "infant_gut_microbiota_dbn_sample_alignment_filtered_sr3d_dbnIntraBoot.graphml")
  #   break
  
}

