library(ape)
library(phangorn)
library(phytools)
library(adephylo)
library(TreeDist)
library(dplyr)
library(TreeTools)
source("G:\\My Drive\\Science\\R programs\\tree functions.r")

#NOTE:Call the package "AwareEnsemble" so you can name the functions "AwareSupport"

#####getDescendantTipNames: return the names of descendants (tips) of a given node number#####
getDescendantTipNames<-function(phylo, nodeNum){unname(slice(as.data.frame(phylo$tip.label), unlist(Descendants(phylo, nodeNum, "tips"))))
}

#####transferSingleNodeSupport: Transfer node support value from one tree to another#####
transferSingleNodeSupport<-function(lowerTeirTree, 
                                    nodeNumber,#this value should be greater than the number of tips but lower than the number of tips +number of nodes
                                    higherTeirTree){
  #get bipartition node value in one tree
  descTipIDs<-Descendants(lowerTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(lowerTeirTree, nodeNumber)
  #descTipNames
  nodeLabel<-lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]
  #find same node in next tree
  higherTreeCorrespondingNode<-getMRCA(higherTeirTree, unlist(as.data.frame(descTipNames)))
  #compare the node support values between the two trees
  lowerTierNodeValue<-lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]
  higherTierNodeValue<-higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]
  
  #replace the node values if the value is >95%
  ifelse(
    as.numeric(lowerTierNodeValue)>=95, 
    {higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<-lowerTierNodeValue
    print(paste(paste(lowerTierNodeValue, " replaces (vs.1) ..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-length(higherTeirTree$tip.label))))
    },{
      #print("Status check 3")
      
      #If the node value is <95% then we need to ensure that the same bipartition is present in both trees...if not we keep the node support value in the second tree
      
      #get the new taxa in the higher tier tree
      
      higherTierTaxa<-setdiff(higherTeirTree$tip.label, lowerTeirTree$tip.label)#a list of all the new taxa added in the higher tier
      
      #compare the tip taxa in the two clades
      #print(higherTreeCorrespondingNode)
      higherDescTipNames<-getDescendantTipNames(higherTeirTree,higherTreeCorrespondingNode )
      
      isTheLowerTierNodeValueValid<-{length(
        intersect(
          setdiff(higherDescTipNames, descTipNames)[[1]], higherTierTaxa
        )
      )==length(
        setdiff(higherDescTipNames, descTipNames)[[1]]
      )} #determines if the additional taxa are all higher tier taxa. If TRUE I should take the value from the lower tier, if FALSE I should take the value from the higher tier.
      
      ifelse(is.na(as.numeric(lowerTierNodeValue)),{print(paste("Node ", as.character(nodeNumber), " is NA. Retaining higher teir value." ))
        isTheLowerTierNodeValueValid<-FALSE}
        , isTheLowerTierNodeValueValid<-isTheLowerTierNodeValueValid)#this is an extra step to check that the node values isn't NA. If it is NA then we force the value FALSE, if it isn't NA then we retain the previous value of TRUE or FALSE
      
      ifelse(isTheLowerTierNodeValueValid,
             {print(paste(paste(lowerTierNodeValue, " replaces (vs.2)..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-length(higherTeirTree$tip.label))))
               higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<-lowerTierNodeValue},{
                 print(paste(paste("Node ",as.character(nodeNumber)), " not the same in both trees. Retaining higher node (",as.character(higherTreeCorrespondingNode),") label:",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]) ))
                 higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<- higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]})}
  )
  
  return(higherTeirTree$node.label)
}

#####transferAllNodeSupportValues: Looping the above over all nodes and support values#####


transferAllNodeSupportValues<-function(lowerTeirTree, higherTeirTree) {
  for (i in (length(lowerTeirTree$tip.label)):(length(lowerTeirTree$tip.label)+lowerTeirTree$Nnode))
  {higherTeirTree$node.label<-transferSingleNodeSupport(lowerTeirTree, i, higherTeirTree )}
  return(higherTeirTree)
}

#####transferSingleNodeSupportAndStatus: Transfer node support value AND status from one tree to another#####

transferSingleNodeSupportAndStatus<-function(
  lowerTeirTree, 
  nodeNumber,#this value should be greater than the number of tips but lower than the number of tips +number of nodes
  higherTeirTree, 
  defaultLowerStatus#A string signifying the status of the lower tier node values
  
){
  #get bipartition node value in one tree
  descTipIDs<-Descendants(lowerTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(lowerTeirTree, nodeNumber)
  #descTipNames
  nodeLabelAndStatus<-lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]
  
  lowerTierNodeValue<-ifelse(grepl("/",nodeLabelAndStatus),
                             strsplit(lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)], "/")[[1]][[1]] #If there is a status in the node label then split it and only store the label
                             , lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]) #if there is only a node value then no need to split the label
  
  #set the lower status for that node
  statusLower<-ifelse(grepl("/",nodeLabelAndStatus), strsplit(nodeLabelAndStatus, "/")[[1]][[2]],defaultLowerStatus )
  
  #find same node in next tree
  higherTreeCorrespondingNode<-getMRCA(higherTeirTree, unlist(as.data.frame(descTipNames)))
  #compare the node support values between the two trees
  
  higherTierNodeValue<-ifelse(grepl("/",higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]),strsplit(higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)], "/")[[1]][[1]] ,higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)])
  
  #replace the node values if the value is >95%
  ifelse(
    as.numeric(lowerTierNodeValue)>=95, 
    {print(paste(paste(lowerTierNodeValue, " replaces (vs.1) ..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-length(higherTeirTree$tip.label))))
      higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<-paste(lowerTierNodeValue, "/", statusLower,sep="")
    },{
      
      #If the node value is <95% then we need to ensure that the same bipartition is present in both trees...if not we keep the node support value in the second tree
      
      #get the new taxa in the higher tier tree
      
      higherTierTaxa<-setdiff(higherTeirTree$tip.label, lowerTeirTree$tip.label)#a list of all the new taxa added in the higher tier
      
      #compare the tip taxa in the two clades
      #print(higherTreeCorrespondingNode)
      higherDescTipNames<-getDescendantTipNames(higherTeirTree,higherTreeCorrespondingNode )
      
      isTheLowerTierNodeValueValid<-{length(
        intersect(
          setdiff(higherDescTipNames, descTipNames)[[1]], higherTierTaxa
        )
      )==length(
        setdiff(higherDescTipNames, descTipNames)[[1]]
      )} #determines if the additional taxa are all higher tier taxa. If TRUE I should take the value from the lower tier, if FALSE I should take the value from the higher tier.
      
      ifelse(is.na(as.numeric(lowerTierNodeValue)),{print(paste("Node ", as.character(nodeNumber), " is NA. Retaining higher teir value." ))
        isTheLowerTierNodeValueValid<-FALSE}
        , isTheLowerTierNodeValueValid<-isTheLowerTierNodeValueValid)#this is an extra step to check that the node values isn't NA. If it is NA then we force the value FALSE, if it isn't NA then we retain the previous value of TRUE or FALSE
      
      ifelse(isTheLowerTierNodeValueValid,
             {print(paste(paste(lowerTierNodeValue, " replaces (vs.2)..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-length(higherTeirTree$tip.label))))
               higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<-paste(lowerTierNodeValue, "/",statusLower,sep="")},
             {print(paste(paste("Node ",as.character(nodeNumber)), " not the same in both trees. Retaining higher node (",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]),") label:",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]) ))
               #higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<- higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]
             })}
  )
  
  return(higherTeirTree$node.label)
}
#####checkUnlabelledNodes: additional function to put a label on unlabeled nodes#####
checkUnlabelledNodes<-function(nodeLabel,defaultHigherStatus ){ifelse(grepl("/",nodeLabel),nodeLabel
                                                                      , ifelse(nodeLabel=="","",paste(nodeLabel,"/" , defaultHigherStatus, sep="")))}
#####transferAllNodeSupportAndStatus: Looping the above over all nodes and support values and add status#####
transferAllNodeSupportAndStatus<-function(lowerTeirTree, higherTeirTree,defaultLowerStatus,defaultHigherStatus) {
  for (i in (length(lowerTeirTree$tip.label)):(length(lowerTeirTree$tip.label)+lowerTeirTree$Nnode))
  {higherTeirTree$node.label<-transferSingleNodeSupportAndStatus(lowerTeirTree, i, higherTeirTree, defaultLowerStatus)
  }
  higherTeirTree$node.label<-unlist(
    lapply(higherTeirTree$node.label,checkUnlabelledNodes, defaultHigherStatus)
  )
  return(higherTeirTree)
}
###########################################################################################################
#####NOTE##################################################################################################
#The above functions go in order from the most conservative (most exclusive) tree to the largest (least exclusive) tree. However, it appears that many of the nodes are incorrectly labeled. A good example of this is Blaberoidea. Blaberoidea is supported in the ABA tree with 100% support but there is no Blaberoidea (SS or SL) node with 100% support (check CDD tree with modified support).
#Below I will try and accomplish the same goal going in the reverse order (starting from the largest tree and working backwards) to see if it doesn't result in a more expected result.
###########################################################################################################

#dependency function
doesLowerCladeHaveAValidSupportLevel<-function(allTierTrees, lowerTierPart, leftoverTaxa){
  binary<-try(
  {
    nodeLabelLower<-allTierTrees[[lowerTierPart]]$node.label[[ getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-length(allTierTrees[[lowerTierPart]]$tip.label)]]
  as.numeric(nodeLabelLower)>=95
  }
  , silent=TRUE)
return(isTRUE(binary))
  }

#dependency function
checkLowerTierTaxaMonophyly<-function(allTierTrees, currentCladeTaxa, lowerTierPart){binary<-try(
  {
    leftoverTaxa<<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
    ifelse(length(leftoverTaxa)>=2, 
    is.monophyletic(
      allTierTrees[[lowerTierPart]],leftoverTaxa), FALSE)
 }# if there are two or more taxa, check if the tips are monophyletic in the lower tier tree #Step 1Aiv
  , silent=TRUE)
return(isTRUE(binary))
}

#dependency function
useWhichNodeLabel<-function(tier, allTierTrees, lowerTierPart, nodeNum, tierLabelHigher, tierLabelLower, leftoverTaxa){
  #currentCladeTaxa<-getDescendantTipNames(allTierTrees[[1]],nodeNum)
  
  
  answer<-switch(tier, 
         higher = paste(allTierTrees[[1]]$node.label[[nodeNum-length(allTierTrees[[1]]$tip.label)]],"/",tierLabelHigher,sep=""),
         lower = paste(allTierTrees[[lowerTierPart]]$node.label[[getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa))-length(allTierTrees[[lowerTierPart]]$tip.label)]],"/",tierLabelLower,sep=""))
  return(answer)
}

#dependency function
checkIfTierShouldbeLowered<-function(allTierTrees, lowerTierPart, currentCladeTaxa){
totalTiers<<-length(allTierTrees)
higherTierTaxa<-unlist(allTierTrees[[1]]$tip.label)
ifelse(lowerTierPart==totalTiers, {print("Reached bottom tier on node.")
  return(FALSE)},{
    print("Checking if lower(-1) taxa are appropriate....")
    ifelse({
      isTRUE(checkLowerTierTaxaMonophyly(allTierTrees, currentCladeTaxa, lowerTierPart+1)) && length(intersect(leftoverTaxa,unlist(allTierTrees[[lowerTierPart+1]]$tip.label) ))>=2    
      
      }
    
    
    , {
    return(TRUE)
    }
    , {print("Lower tier (-1) taxa not appropriate (non-monophyly or taxa not present)")
      return(FALSE)
    })
  })
}



#list some tippy clades function
findTippiestClades<-function(highestTeirTree){
  numTaxaHigh<-length(highestTeirTree$tip.label)
  sisterList<-list(1:numTaxaHigh-1) #create a dummy list to put the values into (for tippiest nodes)
  tipCladeMRCAS<-list(1:numTaxaHigh-1)#create a dummy list to put the values into (for tip-Clade pairs)
  for(i in 2:numTaxaHigh){
    #find the sister taxon to each tip in the tree
    sister<-
      getSisters(highestTeirTree,highestTeirTree$tip.label[[i]])
    
    #store the node ID of the MRCA for each tip-tip pair. For each tip-clade pair, store NA, which we will remove later
    sisterList[[i-1]]<-ifelse(sister<=numTaxaHigh,getMRCA(highestTeirTree,
                                                          c(highestTeirTree$tip.label[[i]],     
                                                            highestTeirTree$tip.label[[sister]])), NA)
  }
  
  
  #find the tippiest clades and some other clades
  for(i in 2:numTaxaHigh){
    #find the sister taxon to each tip in the tree
    sister<-
      getSisters(highestTeirTree,highestTeirTree$tip.label[[i]])
    
    #store the node ID of the MRCA for each tip-tip pair. For each tip-clade pair, store NA, which we will remove later
    tipCladeMRCAS[[i-1]]<-ifelse(sister>numTaxaHigh,sister, NA)
  }
  
  sisterList<-unlist(sisterList)
  tipCladeMRCAS<-unlist(tipCladeMRCAS)
  tippiestNodes<<-unique(sisterList[!is.na(sisterList)])#this is a list of the nodeIDs I will first transfer between tiers
  tipCladeMRCAS<<-unique(tipCladeMRCAS[!is.na(tipCladeMRCAS)])#this is a list of some of the other nodes I will transfer, but later on
  }


#######for a single node, find which tree has the valid node support value#####
evaluateNodeSupport<-function(
                              allTierTrees, #this object should have the trees in order from highest tier to lowest tier (i.e. most inclusive to least inclusive)
                              tierLabels, #the names of the three trees in the same order as the trees
                              nodeNumber, #the number of the higher tier node we are evaluating
                              nodeIterator, #an arbitrary sequential iterator used to keep track of where to store the labeled node IDs
                              verbose = FALSE,
                              developmentMode=FALSE
){
  if(verbose==TRUE) {print(paste("Starting on highest tier node", as.character(nodeNumber)))}
  highestTeirTree<-allTierTrees[[1]]
  numTaxaHigh<-length(highestTeirTree$tip.label)
  
  #Read a single node
  descTipIDs<-Descendants(highestTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(highestTeirTree, nodeNumber)
  #Quality check 1
  if(verbose==TRUE) {print(ifelse(
    {length(setdiff(as.list(as.data.frame(highestTeirTree$tip.label)), as.list(descTipNames)))==0}, "Error: This is the root", "Print:QC1 passed" ))
  }
  

  currentNode<-nodeNumber
  currentCladeTaxa<-getDescendantTipNames(highestTeirTree,currentNode)
  lowerTierPart<-2 #this is just the initial point for lowerTierPart. It will be incremented below under the right circumstances
  leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
  if(verbose==TRUE){print(paste("Lower teir part set to initial value: ", as.character(lowerTierPart)))}
  
  cladeLabeled<-FALSE #set this to FALSE so that when this is evaluated I can exit the loop when the clade is correctly labeled
  
  
  
  while(cladeLabeled==FALSE){
    higherTierTaxa<-setdiff(highestTeirTree$tip.label,allTierTrees[[lowerTierPart]]$tip.label)
    
    ##precalculate some conditionals to make the below code a little easier to read.
    areCladeTaxaAllHigher<-length(
      unlist(currentCladeTaxa) #how many total tips are in the clade considered
    )-(length(
      intersect(unlist(currentCladeTaxa), higherTierTaxa) #how many of those tips are higher tier taxa
    ))<=1
    
    doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart, leftoverTaxa)##DOES THIS NEED TO BE HERE?
    
    always<-"NEVER"#This object should never change in value. This is only to allow the while chunk below to run properly. It will be exited when the outer loop is exited, or when a "break" occurs.
    while(always=="NEVER"){
      
          ifelse(areCladeTaxaAllHigher
           #Step 1Ai-iii
           ,#if all the taxa are higher tier taxa
           {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[1]], tierLabels[[lowerTierPart]], leftoverTaxa)
           if(verbose==TRUE){print("Use higher tier label (cond. 1)")}
           cladeLabeled<-TRUE},
           #if NOT... 
           #
           ifelse(checkLowerTierTaxaMonophyly(allTierTrees, currentCladeTaxa, lowerTierPart), #label: 4
                  ifelse(checkIfTierShouldbeLowered(allTierTrees, lowerTierPart, currentCladeTaxa), {
                    lowerTierPart<-lowerTierPart+1
                    #lowerTierPart
                    if(verbose==TRUE){print(paste("Tier lowered by 1. Now = ",as.character(lowerTierPart),". Restarting loop."))}
                    break} #if the tier needs to be lowered then we must escape the loop and start over again
                         , 
                     
                    {
                      if(verbose==TRUE){print(paste("Working on tier", as.character(tierLabels[[lowerTierPart]])))}
                      
                           ifelse({leftoverTaxa<- intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
                             check<-doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart,leftoverTaxa )
                             
                           check},{raiseNode<-TRUE
                             while(raiseNode==TRUE){
                                       ifelse(
                               is.element(getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa)), finishedNodes[lowerTierPart][[1]][!is.na(  finishedNodes[lowerTierPart][[1]])] ) #check if the node on the lower tier tree has been used in a label already
                                    ,{ print("Label used previously. Moving up one tier and rechecking.")         #if so, change the tier to a higher one
                                      lowerTierPart<-lowerTierPart-1
                                      leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
                                      ifelse(doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart, leftoverTaxa)||(lowerTierPart==1), 
                                            print("Higher label valid."),
                                            {print("Higher (+1) label invalid. Raising tier again.")
                                              lowerTierPart<-lowerTierPart-1
                                               leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
                                              }
                                             )
                                      
                                      } 
                               , #if not, continue with this label/tier
                                    { print("QC Check 3: Good")
                                      finishedNodes[[lowerTierPart]][[nodeIterator]]<<-getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa))
                                    highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("lower", allTierTrees, lowerTierPart, currentNode, tierLabels[[1]], tierLabels[[lowerTierPart]], leftoverTaxa)
                                    print(paste("Using",as.character(lowerTierPart),"tier label"))
                                    cladeLabeled<-TRUE
                                    raiseNode<-FALSE #exit the loop of raising the tier level
                                    })
                               print("QC Check 4: Good")
                             }#END OF while raiseNode==TRUE
                           print("QC Check 5: Good")
                             },
                                  
                                  {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[1]], tierLabels[[lowerTierPart]], leftovertaxa)
                                  print("Use higher tier label (cond. 3)")
                                  cladeLabeled<-TRUE}
                           )}
                    )#label 4: if true
                  ,
                         {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[1]], tierLabels[[lowerTierPart]], leftoverTaxa)
                         print("Use higher tier label (cond. 2)")
                         cladeLabeled<-TRUE}#label 4: if false
                    )
            )
           
    break}
  }
  output<-unlist(
    c(
    highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]],
    unlist(lowerTierPart[[1]]), 
    
    {
      mrca<-getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa))
      ifelse(is.null(mrca), NA, mrca)
    }
  ))
  return(output)
  }


#Determine the order of nodes to evaluate on the higher tree
findOrderToTraverseHighestTree<-function(allTierTrees){
  findTippiestClades(allTierTrees[[1]])
  
  highestTierNodesOrdered<-unlist(as.list(matrix(data = NA, nrow =((length(allTierTrees[[1]]$tip.label)+length(allTierTrees[[1]]$node.label))-(length(allTierTrees[[1]]$tip.label))), ncol=1 )))
  
  #add the root to the list
  highestTierNodesOrdered[[1]]<-length(allTierTrees[[1]]$tip.label)+1
  
  highestTeirTree<-allTierTrees[[1]]
  iterator<-1
  ancestralNode<-"Initial"
  
  for(i in 1:length(tippiestNodes)){#this loop will go through each tippy node once,
    #print(as.character(i))
    iterator<-iterator+1 #increment the iterator so the next node values goes in the next spot on the list
    currentNode<-tippiestNodes[[i]]#store the tippy node as the first node
    ancestralNode<-currentNode
    while(!is.element(ancestralNode,highestTierNodesOrdered )){
      #print(as.character(currentNode))
      highestTierNodesOrdered[[iterator]]<-currentNode #add the next node to the list. This actually starts at position 2 because we have previously defined position 1 to be the root node.
      iterator<-iterator+1 #increment the iterator so the next node values goes in the next spot on the list
      ancestralNode<-getMRCA(highestTeirTree, c(getSisters(highestTeirTree, currentNode), currentNode)) #find the MRCA of the current node (tippy node if this is the 1st iteration) and it's sister clade
      currentNode<-ancestralNode #store that as the new current node
      #ifelse(currentNode==highestTierNodesOrdered[[1]], print("Root reached. Restarting"), print("..."))     
      
    } #if the current node is already on the list of ordered nodes (highestTierNodesOrdered) then we will move onto the next tippy node
  }
  #returned<-highestTierNodesOrdered[!is.na(highestTierNodesOrdered)]
  return(highestTierNodesOrdered[!is.na(highestTierNodesOrdered)])
  
  
  
}

#Final wrapper function that evaluates the whole tree.
totalAwareSupport<-function(allTierTrees, tierNames, verbose=FALSE){
  highestTierTree<<-allTierTrees[[1]]
  nodeList<-findOrderToTraverseHighestTree(allTierTrees) #find the order to traverse the highest (most inclusive) tree
  nodeIterator<-0
  
  finishedNodes<<-data.frame(matrix(NA,nrow = length(nodeList),ncol = length(allTierTrees))) #this is a high level variable that will be used to keep track of which nodes have already been tested
  
  for(j in nodeList[-1]){
    
    nodeIterator<-nodeIterator+1
    ifelse(is.element(j, finishedNodes[[1]]), print("Error: Annotating higher tier tree node twice"), print("QC Check 2: Good."))
    
    answer<-evaluateNodeSupport(allTierTrees, tierNames,j,nodeIterator, verbose)
    
    highestTierTree$node.label[[j-length(highestTierTree$tip.label)]]<-answer[[1]]
    finishedNodes[[1]][[nodeIterator]]<-j
    finishedNodes[[as.numeric(answer[[2]])  ]][[nodeIterator]]<-answer[[3]]

  }
  return(list(highestTierTree, finishedNodes))
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
           graftedTree$node.label[[nodeID-length(finalTree$node.label)]]<-paste(finalTree$node.label[[nodeID-length(finalTree$node.label)]], #original node support value
                                                                                sourceTree$node.label[[correspondingNodeID-length(sourceTree$tip.label)]], sep="|") #new node support value
           , NA)
    
  }
}






####October 2022 NOTE: There is a problem with numbering of the clades (I think)
#### and annotation of support values. This can be seen in the Blaberidae clade with Hypnosphaeria.











