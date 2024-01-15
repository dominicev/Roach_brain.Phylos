
library(ape)
library(phangorn)
library(phytools)
library(adephylo)
library(fitdistrplus)
library(TreeSearch)
library(reshape2)
library(dplyr)
library(MonoPhy)
library(data.table)
library(foreach)
library(doParallel)
library(TreeTools)
library(castor)
library(devtools)

source_url("https://raw.githubusercontent.com/dominicev/Roach_brain_Phylos/main/createNclusterFile.r")
source_url("https://raw.githubusercontent.com/dominicev/Roach_brain_Phylos/main/deltaGLS.r")
source_url("https://raw.githubusercontent.com/dominicev/Roach_brain_Phylos/main/phylo_outlier_trimming.r")


###Unused packages###
#library(rlist)
    #Unused packages determined by...
    #funchir::stale_package_check('B:/OneDrive - University of Illinois - Urbana/Science/R programs/Roach_brain_Phylos/tree functions.R')


#use this to load this file into an r book
#source_url("https://raw.githubusercontent.com/dominicev/Roach_brain_Phylos/main/tree functions.R")
 
##Reminder: Explanation of subsetting lists/dataframes in R: https://adv-r.hadley.nz/subsetting.html


####rooting a tree from a list of possible outgroups####
###
###This function will root a tree given a list of possible outgroups. This is useful when rooting a series of incomplete gene trees
###
rootIncompleteTree<-function(phylo, outGroupList) {
  i<-0 #initialize iterator
  repeat{
    i=i+1 #increment iterator
    max=length(outGroupList)
    if(outGroupList[i] %in% phylo$tip.label | i>max ){ #it will stop the repeat loop if the outgroup tip appears in the tree or if it has gone through all the outgroups and no tip appears in the tree
      outgroup2Root=i
      break}}
  ifelse(i<max,
         return(root(phylo,outGroupList[outgroup2Root], resolve.root=TRUE )),
         return(midpoint.root(phylo))) #if there is an outgroup to root at then this returns the appropriately rooted tree, otherwise it does a midpoint root
  
}


####filter outlier and intruder tips in a tree by a single clade####

filterByOneClade<-function(phylo,cladeControlList,cladeControlListStrict,taxaInTree=phylo$tip.label ){
  
  # example    cladeControlList=taxonList$Blaberidae_control
  # example    cladeControlListStrict=sppInTree$Blaberidae_control
  # example    taxaInTree = sppInTree$Taxon
  
  ##get only the species that are in the tree from the control clade list
  taxonLs<-cladeControlList[taxonList[[2]] %in%  phylo$tip.label]
  
  ##replace unknowns (?) with NA
  controlBooleans<-unlist(rapply(as.list(taxonLs),function(x) ifelse(x=="?",NA,x), how = "replace"))
  
  ##make lists of species in each clade
  sppInClade<-taxonList$Taxon[controlBooleans>0]
  sppOutsideClade<-taxonList$Taxon[controlBooleans<1]
  unplacedTaxa<-taxonList$Taxon[taxonLs=="?"]
  
  ##make monophyclade table
  cladeTable<-data.frame(taxaInTree, cladeControlListStrict) 
  
  ##Assess monophyly
  AMObject<-AssessMonophyly(rootIncompleteTree(phylo, outgroupsOrdered),cladeTable , outlierlevel = .75)
  outliers<-AMObject$cladeControlListStrict$OutlierTips$`1`
  intruders<-AMObject$cladeControlListStrict$IntruderTips$`1`
  
  
  ##filter results based on unplaced taxa
  outliersFiltered<-outliers[!(outliers %in% unplacedTaxa)]#this is the list of outliers with the unplacedTaxa removed
  intrudersFiltered<-intruders[!(intruders %in% unplacedTaxa)]#this is the list of intruders with the unplacedTaxa removed
  #c(intruders,"",intrudersFiltered)
  #c(outliers,"",outliersFiltered)
  
  
  ##filter tree
  return(c(drop.tip(drop.tip(unroot(phylo), outliersFiltered), intrudersFiltered), outliersFiltered, intrudersFiltered))
}



####this is a wrapper function for the above. it can be applied to a tree file and looped over multiple tree files####
####
####taxonList is a data.frame with the form found in "gene tree filtering.rmd"/ G:\My Drive\Projects\PHY_Blab_2\bioinformatics\PhyBlab2 full analysis\taxonList.csv
####
####outgroupsOrdered is just a list of taxa for prioritized rooting


filterTaxaFromTrees<-function(treeFileName, taxonList, outgroupsOrdered){
  
  #Import data
  phylo<-read.tree(treeFileName)
  
  #list the species in the tree
  sppInTree<-filter(taxonList, taxonList[[2]] %in%  phylo$tip.label)
  
  #make a list to loop over
  allControlVariablesList<-data.frame(taxonList$Nyctiboridae_control,taxonList$Blaberidae_control, taxonList$Oxyhaloninae_control,taxonList$Solumblattodea_Control, taxonList$Tutricablattae_control, taxonList$Blattodea_control, taxonList$Mantodea_control, taxonList$Dictyoptera_control)
  
  allControlVariablesStrictList<-data.frame(sppInTree$Nyctiboridae_control,sppInTree$Blaberidae_control, sppInTree$Oxyhaloninae_control,sppInTree$Solumblattodea_Control, sppInTree$Tutricablattae_control, sppInTree$Blattodea_control, sppInTree$Mantodea_control, sppInTree$Dictyoptera_control)
  
  print(treeFileName)
  newPhylo<-phylo
  
  for(i in 1:length(allControlVariablesList)){
    sppInTree<-filter(taxonList, taxonList[[2]] %in%  newPhylo$tip.label)
    allControlVariablesStrictList<-data.frame(sppInTree$Nyctiboridae_control,sppInTree$Blaberidae_control, sppInTree$Oxyhaloninae_control,sppInTree$Solumblattodea_Control, sppInTree$Tutricablattae_control, sppInTree$Blattodea_control, sppInTree$Mantodea_control, sppInTree$Dictyoptera_control)
    
    #print(i)
    
    newPhylo<-filterByOneClade(newPhylo, allControlVariablesList[[i]], allControlVariablesStrictList[[i]], sppInTree$Taxon)[1]
    
  }
  print(paste(length(phylo$tip.label)-length(newPhylo$tip.label),"out of",length(phylo$tip.label), "total tips dropped"))
  
  write.tree(newPhylo,paste(strsplit(treeFileName, split=".", fixed = TRUE)[[1]][[1]],".filteredSp.tre", sep=""))
}


filterTaxaStatsByTree<-function(treeFileName, taxonList, outgroupsOrdered){
  
  #Import data
  phylo<-read.tree(treeFileName)
  
  #list the species in the tree
  sppInTree<-filter(taxonList, taxonList[[2]] %in%  phylo$tip.label)
  
  #make a list to loop over
  allControlVariablesList<-data.frame(taxonList$Nyctiboridae_control,taxonList$Blaberidae_control, taxonList$Oxyhaloninae_control,taxonList$Solumblattodea_Control, taxonList$Tutricablattae_control, taxonList$Blattodea_control, taxonList$Mantodea_control, taxonList$Dictyoptera_control)
  
  allControlVariablesStrictList<-data.frame(sppInTree$Nyctiboridae_control,sppInTree$Blaberidae_control, sppInTree$Oxyhaloninae_control,sppInTree$Solumblattodea_Control, sppInTree$Tutricablattae_control, sppInTree$Blattodea_control, sppInTree$Mantodea_control, sppInTree$Dictyoptera_control)
  
  print(treeFileName)
  newPhylo<-phylo
  
  for(i in 1:length(allControlVariablesList)){
    sppInTree<-filter(taxonList, taxonList[[2]] %in%  newPhylo$tip.label)
    allControlVariablesStrictList<-data.frame(sppInTree$Nyctiboridae_control,sppInTree$Blaberidae_control, sppInTree$Oxyhaloninae_control,sppInTree$Solumblattodea_Control, sppInTree$Tutricablattae_control, sppInTree$Blattodea_control, sppInTree$Mantodea_control, sppInTree$Dictyoptera_control)
    
    #print(i)
    
    newResults<-filterByOneClade(newPhylo, allControlVariablesList[[i]], allControlVariablesStrictList[[i]], sppInTree$Taxon)
    newPhylo<-newResults[[1]]
    
  }
  print(paste(length(phylo$tip.label)-length(newPhylo$tip.label),"out of",length(phylo$tip.label), "total tips dropped"))
  
  #write.tree(newPhylo,paste(strsplit(treeFileName, split=".", fixed = TRUE)[[1]][[1]],".filteredSp.tre", sep=""))
  #
  c(strsplit(treeFileName, split=".", fixed = TRUE)[[1]][[1]], newResults[2], newResults[3])
  
}

######Move one node up on a tree######

nodeMRCA<-function(phylo, nodeID){
  ifelse(
    is.numeric(try(as.numeric(nodeID) , silent=TRUE)), nodeID<-nodeID, nodeID<-as.character(match(nodeID,phylo$tip.label ))
    )  
                   
  sister<-getSisters(phylo, nodeID)
  return(getMRCA(phylo, unlist(c(sister, nodeID)) ))
  
  }



#####getDescendantTipNames: return the names of descendants (tips) of a given node number#####
getDescendantTipNames<-function(phylo, nodeNum){unname(slice(as.data.frame(phylo$tip.label), unlist(Descendants(phylo, nodeNum, "tips"))))
}


#####collapseUnsupportedNodes: collapse nodes in a tree based on their node label (i.e., bootstrap support)
collapseUnsupportedNodes<-function(phy, cutoff = 75){
  CollapseNode(phy, length(phy$tip.label)+match(TRUE, phy$node.label<cutoff))
  
}


#####Add node support to a tree from another tree

concatNodeLabels<-function(finalTree, sourceTree){
  #get a list of all the nodes to be annotated
  finalTreeNodes<-findOrderToTraverseHighestTree(c(finalTree))
  
  for(nodeID in finalTreeNodes[-1]){
    
    #the corresponding node in the concatenation tree
    correspondingNodeID<-try(getMRCA(sourceTree, unlist(getDescendantTipNames(finalTree, nodeID))), silent=TRUE)
    #print(paste(nodeID, correspondingNodeID))
    #check that the node is monophyletic in the corresponding tree
    takeNodeValueQ<-checkLowerTierTaxaMonophyly(c(finalTree, sourceTree),  
                                                unlist(getDescendantTipNames(finalTree,nodeID))
                                                , 2)
    ifelse(takeNodeValueQ, 
           finalTree$node.label[[nodeID-Nnode(finalTree)]]<-paste(finalTree$node.label[[nodeID-Nnode(finalTree)]], #original node support value
                                                                  sourceTree$node.label[[correspondingNodeID-Ntip(sourceTree)]], sep="|") #new node support value
           , NA)
    
  }
}



###The above - modified by GPT4
concatNodeLabels <- function(finalTree, sourceTree) {
  # Get a list of all the nodes to be annotated
  finalTreeNodes <- findOrderToTraverseHighestTree(c(finalTree))
  
  for(nodeID in finalTreeNodes[-1]) {
    # The corresponding node in the concatenation tree
    correspondingNodeID <- try(getMRCA(sourceTree, unlist(getDescendantTipNames(finalTree, nodeID))), silent = TRUE)
    
    # Check that the node is monophyletic in the corresponding tree
    takeNodeValueQ <- checkLowerTierTaxaMonophyly(c(finalTree, sourceTree), unlist(getDescendantTipNames(finalTree, nodeID)), 2)
    
    ifelse(takeNodeValueQ, 
           finalTree$node.label[[nodeID - Nnode(finalTree)]] <- paste(finalTree$node.label[[nodeID - Nnode(finalTree)]], # Original node support value
                                                                      sourceTree$node.label[[correspondingNodeID - Ntip(sourceTree)]], sep = "|") # New node support value
           , NA)
  }
  
  return(finalTree) # Return the modified tree
}


###the above written from scratch by GPT 4
concatNodeLabels3 <- function(tree1, tree2) {
  # Iterate through the nodes of tree1
  for (node in (Ntip(tree1) + 1):(Ntip(tree1) + Nnode(tree1))) {
    # Get the descendant tips for the current node in tree1
    tips1 <- tree1$tip.label[getDescendants(tree1, node)]
    
    # Find the corresponding node in tree2
    corresponding_node <- NA
    for (node2 in (Ntip(tree2) + 1):(Ntip(tree2) + Nnode(tree2))) {
      tips2 <- tree2$tip.label[getDescendants(tree2, node2)]
      if (all(tips1 %in% tips2) && all(tips2 %in% tips1)) {
        corresponding_node <- node2
        break
      }
    }
    
    # If a corresponding node is found, concatenate the node labels
    if (!is.na(corresponding_node)) {
      tree1$node.label[node - Ntip(tree1)] <- paste(tree1$node.label[node - Ntip(tree1)], tree2$node.label[corresponding_node - Ntip(tree2)], sep = "|")
    }
  }
  
  return(tree1) # Return the modified tree1
}



###the above written by GPT 4 to be faster
concatNodeLabelsFAST <- function(tree1, tree2) {
  # Precompute descendants for both trees
  descendants1 <- lapply((Ntip(tree1) + 1):(Ntip(tree1) + Nnode(tree1)), function(node) tree1$tip.label[getDescendants(tree1, node)])
  descendants2 <- lapply((Ntip(tree2) + 1):(Ntip(tree2) + Nnode(tree2)), function(node) tree2$tip.label[getDescendants(tree2, node)])
  
  # Iterate through the nodes of tree1
  for (node in seq_along(descendants1)) {
    tips1 <- descendants1[[node]]
    
    # Find the corresponding node in tree2
    corresponding_node <- which(sapply(descendants2, function(tips2) all(tips1 %in% tips2) && all(tips2 %in% tips1)))
    
    node1 <- node + Ntip(tree1)
    
    # If a corresponding node is found, concatenate the node labels
    if (length(corresponding_node) > 0) {
      node2 <- corresponding_node + Ntip(tree2)
      tree1$node.label[node1 - Ntip(tree1)] <- paste(tree1$node.label[node1 - Ntip(tree1)], tree2$node.label[node2 - Ntip(tree2)], sep = "|")
    } else {
      # If no corresponding node is found, add "NA"
      tree1$node.label[node1 - Ntip(tree1)] <- paste(tree1$node.label[node1 - Ntip(tree1)], "NA", sep = "|")
    }
  }
  
  return(tree1) # Return the modified tree1
}





###This is a custom function that evaluates input that looks like Mathematica's Table function and provides a similar output style, but in R!
###Note: This isn't useful, I will never use this :-)

tableW <- function(expr, var_list) {
  # Extract the variable, start, end, and increment
  var_name <- var_list[[1]]
  start <- var_list[[2]]
  end <- var_list[[3]]
  increment <- ifelse(length(var_list) > 3, var_list[[4]], 1)
  
  # Create a vector to store the results
  result <- vector("list", length = (end - start) / increment + 1)
  
  # Evaluate the expression for each value of the variable
  for (i in seq(start, end, by = increment)) {
    assign(var_name, i, envir = .GlobalEnv)
    result[[i - start + 1]] <- eval(parse(text = expr))
  }
  
  # Return the result as a matrix or array if possible
  return(do.call(rbind, result))
}

# Example usage
#expr <- "x^2"
#var_list <- list("x", 1, 5)
#result <- tableW(expr, var_list)
#print(result)




###Match tips between two trees with different tip sets. GPT-4
reciprocalTreeTrim <- function(tree1, tree2) {
  # Find the common tips
  common_tips <- intersect(tree1$tip.label, tree2$tip.label)
  
  # Find the tips to drop from each tree
  drop_tips_tree1 <- setdiff(tree1$tip.label, common_tips)
  drop_tips_tree2 <- setdiff(tree2$tip.label, common_tips)
  
  # Drop the extra tips from each tree
  if (length(drop_tips_tree1) > 0) {
    tree1 <- drop.tip(tree1, drop_tips_tree1)
  }
  if (length(drop_tips_tree2) > 0) {
    tree2 <- drop.tip(tree2, drop_tips_tree2)
  }
  
  # Return the trimmed trees
  return(list(tree1, tree2))
}



# Code generated on 2023-09-01 by GPT-4
reciprocalTreeTrimMulti <- function(tree1, tree2) {
  # Initialize variables
  is_multiPhylo_tree1 <- inherits(tree1, "multiPhylo")
  is_multiPhylo_tree2 <- inherits(tree2, "multiPhylo")
  
  # Function to drop tips from a single tree
  dropTipsFromTree <- function(tree, common_tips) {
    drop_tips <- setdiff(tree$tip.label, common_tips)
    if (length(drop_tips) > 0 && length(drop_tips) != length(tree$tip.label)) {
      return(drop.tip(tree, drop_tips))
    } else if (length(drop_tips) == length(tree$tip.label)) {
      return(NULL)
    }
    return(tree)
  }
  
  # Function to get common tips between a single tree and a multiPhylo object
  getCommonTipsMulti <- function(tree, multi_tree) {
    multi_tips <- unique(unlist(lapply(multi_tree, function(x) x$tip.label)))
    return(intersect(tree$tip.label, multi_tips))
  }
  
  # Handle each combination of multiPhylo and single phylo objects
  if (is_multiPhylo_tree1 && is_multiPhylo_tree2) {
    # Handle the case where both are multiPhylo objects (not covered in this example)
  } else if (is_multiPhylo_tree1) {
    common_tips <- getCommonTipsMulti(tree2, tree1)
    tree1 <- lapply(tree1, function(t) dropTipsFromTree(t, common_tips))
    tree1 <- Filter(Negate(is.null), tree1)
    tree2 <- dropTipsFromTree(tree2, common_tips)
  } else if (is_multiPhylo_tree2) {
    common_tips <- getCommonTipsMulti(tree1, tree2)
    tree2 <- lapply(tree2, function(t) dropTipsFromTree(t, common_tips))
    tree2 <- Filter(Negate(is.null), tree2)
    tree1 <- dropTipsFromTree(tree1, common_tips)
  } else {
    common_tips <- intersect(tree1$tip.label, tree2$tip.label)
    tree1 <- dropTipsFromTree(tree1, common_tips)
    tree2 <- dropTipsFromTree(tree2, common_tips)
  }
  
  # Return the trimmed trees
  return(list(tree1, tree2))
}



#this function is a wrapper for TreeSearch's SPRMoves. From GPT4
rSPRNew <- function(tree, N) {
  # Perform up to N SPR moves
  for (i in 1:N) {
    # Perform a single SPR move
    tree <- SPR(tree)
    
    # Force the tree to be binary
    tree <- multi2di(collapse.singles(tree))
    
    # Optional: Check if the tree is binary (should always be true)
    #if (!is.binary(tree)) {
    #  stop("Tree is not binary after move", i)
    #}
  }
  
  return(tree)
}



##A function to calculate external to internal branch length differences##
leafiness<-function(phy, method="median"){
  if(class(phy)=="phylo"){tree=phy}else{tree=read.tree(phy)}
  
  allD=dist.nodes(tree)
  internalDs<-allD[(length(tree[[4]])+1):(length(tree[[2]])+1),(length(tree[[4]])+1):(length(tree[[2]])+1)]
  #totalDs=allD[1:length(tree[[4]]), 1:length(tree[[4]])]
  
  allLeafdists<-foreach(i=1:length(tree[[4]]))%do%{
    
    #fir=tree[[4]][[i]]
    sec=getSisters(tree,tree[[4]][[i]],mode="number" )[[1]]
    mrca=getMRCA(tree,c(i, sec)  )#this returns the node number at the base of the specified tip by finding the shared node between the tip and it's sister taxon
    dist.nodes(tree)[i,mrca]
  }
  switch(method,  
         mean=mean(unlist(allLeafdists))/mean(internalDs), #the ratio of external to internal branch lengths
         median=median(unlist(allLeafdists))/median(internalDs)
  )
}


leafiness2 <- function(phy, method="median") {
  if(class(phy) == "phylo") { tree = phy } else { tree = read.tree(phy) }
  allD = dist.nodes(tree)
  len4 = length(tree[[4]])
  len2 = length(tree[[2]])
  internalDs <- allD[(len4 + 1):(len2 + 1), (len4 + 1):(len2 + 1)]
  
  allLeafdists <- foreach(i = 1:len4) %do% {
    sec = getSisters(tree, tree[[4]][[i]], mode="number")[[1]]
    mrca = getMRCA(tree, c(i, sec))
    dist.nodes(tree)[i, mrca]
  }
  
  leaf_dists = unlist(allLeafdists)
  switch(method,  
         mean = mean(leaf_dists) / mean(internalDs),
         median = median(leaf_dists) / median(internalDs)
  )
}



##A function to calculate internal to external branch length ratio##
#aka the inverse of the above function

steminess<-function(phy, method="median"){
  1/leafiness(phy, method)
}




#This function takes a multiPhylo object as input, compares all the node support values and 
#returns a single phylo object that has the highest of all the node values (across multiple
#trees) on each node. Made by GPT4



getHighestNodeValues <- function(multi_trees) {
  # Check if the input is a multiPhylo object
  if (!inherits(multi_trees, "multiPhylo")) {
    stop("Input must be a multiPhylo object")
  }
  
  # Check that all trees have the same topology
  base_edges <- multi_trees[[1]]$edge
  if (any(sapply(multi_trees[-1], function(tree) !identical(tree$edge, base_edges)))) {
    stop("All trees must have the same topology")
  }
  
  # Use the first tree as the base
  result_tree <- multi_trees[[1]]
  
  # Iterate through the nodes
  for (node in 1:Nnode(result_tree)) {
    # Get the node values for the current node across all trees
    node_values <- sapply(multi_trees, function(tree) as.numeric(tree$node.label[node]))
    
    # Set the node value in the result tree to the highest value
    result_tree$node.label[node] <- max(node_values, na.rm = TRUE)
  }
  
  return(result_tree)
}

                                    
