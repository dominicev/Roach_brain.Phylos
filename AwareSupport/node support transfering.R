library(ape)
library(phangorn)
library(phytools)
library(adephylo)
library(TreeDist)
library(dplyr)
library(TreeTools)
library(RRphylo)
library(devtools)
#devtools::source_url("https://github.com/dominicev/Roach_brain_Phylos/blob/main/tree%20functions.R")
source("B:/OneDrive - University of Illinois - Urbana/Science/R programs/Roach_brain_Phylos/tree functions.R")
#options(show.error.locations = TRUE)
library(ips)



#####getDescendantTipNames: return the names of descendants (tips) of a given node number#####
getDescendantTipNames<-function(phylo, nodeNum){unname(slice(as.data.frame(phylo$tip.label), unlist(Descendants(phylo, nodeNum, "tips"))))
}

#####transferSingleNodeSupport: Transfer node support value from one tree to another#####
transferSingleNodeSupport<-function(lowerTeirTree, 
                                    nodeNumber,#this value should be greater than the number of tips but lower than the number of tips +number of nodes
                                    higherTeirTree, supportCutOff = 95){
  #get bipartition node value in one tree
  descTipIDs<-Descendants(lowerTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(lowerTeirTree, nodeNumber)
  #descTipNames
  nodeLabel<-lowerTeirTree$node.label[nodeNumber-Ntip(lowerTeirTree)]
  #find same node in next tree
  higherTreeCorrespondingNode<-getMRCA(higherTeirTree, unlist(as.data.frame(descTipNames)))
  #compare the node support values between the two trees
  lowerTierNodeValue<-lowerTeirTree$node.label[nodeNumber-Ntip(lowerTeirTree)]
  higherTierNodeValue<-higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]
  
  #replace the node values if the value is >cutOff%
  ifelse(
    as.numeric(lowerTierNodeValue)>=supportCutOff, 
    {higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]<-lowerTierNodeValue
    print(paste(paste(lowerTierNodeValue, " replaces (vs.1) ..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-Ntip(higherTeirTree))))
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
             {print(paste(paste(lowerTierNodeValue, " replaces (vs.2)..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-Ntip(higherTeirTree))))
               higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<-lowerTierNodeValue},{
                 print(paste(paste("Node ",as.character(nodeNumber)), " not the same in both trees. Retaining higher node (",as.character(higherTreeCorrespondingNode),") label:",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]) ))
                 higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]<- higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]})}
  )
  
  return(higherTeirTree$node.label)
}

#####transferAllNodeSupportValues: Looping the above over all nodes and support values#####


transferAllNodeSupportValues<-function(lowerTeirTree, higherTeirTree, supportCutOff = 95) {
  for (i in (Ntip(lowerTeirTree)):(Ntip(lowerTeirTree)+lowerTeirTree$Nnode))
  {higherTeirTree$node.label<-transferSingleNodeSupport(lowerTeirTree, i, higherTeirTree , supportCutOff)}
  return(higherTeirTree)
}

#####transferSingleNodeSupportAndStatus: Transfer node support value AND status from one tree to another#####

transferSingleNodeSupportAndStatus<-function(
  lowerTeirTree, 
  nodeNumber,#this value should be greater than the number of tips but lower than the number of tips +number of nodes
  higherTeirTree, 
  defaultLowerStatus,#A string signifying the status of the lower tier node values
  supportCutOff = 95 #The bootstrap value you consider to be "high support" that you want to be carried through to next tier
){
  #get bipartition node value in one tree
  descTipIDs<-Descendants(lowerTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(lowerTeirTree, nodeNumber)
  #descTipNames
  nodeLabelAndStatus<-lowerTeirTree$node.label[nodeNumber-Ntip(lowerTeirTree)]
  
  lowerTierNodeValue<-ifelse(grepl("/",nodeLabelAndStatus),
                             strsplit(lowerTeirTree$node.label[nodeNumber-Ntip(lowerTeirTree)], "/")[[1]][[1]] #If there is a status in the node label then split it and only store the label
                             , lowerTeirTree$node.label[nodeNumber-Ntip(lowerTeirTree)]) #if there is only a node value then no need to split the label
  
  #set the lower status for that node
  statusLower<-ifelse(grepl("/",nodeLabelAndStatus), strsplit(nodeLabelAndStatus, "/")[[1]][[2]],defaultLowerStatus )
  
  #find same node in next tree
  higherTreeCorrespondingNode<-getMRCA(higherTeirTree, unlist(as.data.frame(descTipNames)))
  #compare the node support values between the two trees
  
  higherTierNodeValue<-ifelse(grepl("/",higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]),strsplit(higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)], "/")[[1]][[1]] ,higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)])
  
  #replace the node values if the value is >supportCutOff%
  ifelse(
    as.numeric(lowerTierNodeValue)>=supportCutOff, 
    {print(paste(paste(lowerTierNodeValue, " replaces (vs.1) ..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-Ntip(higherTeirTree))))
      higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]<-paste(lowerTierNodeValue, "/", statusLower,sep="")
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
             {print(paste(paste(lowerTierNodeValue, " replaces (vs.2)..."), higherTierNodeValue, " on node ", nodeNumber, "/", as.character(higherTreeCorrespondingNode-Ntip(higherTeirTree))))
               higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]<-paste(lowerTierNodeValue, "/",statusLower,sep="")},
             {print(paste(paste("Node ",as.character(nodeNumber)), " not the same in both trees. Retaining higher node (",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]),") label:",as.character(higherTeirTree$node.label[higherTreeCorrespondingNode-Ntip(higherTeirTree)]) ))
               #higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]<- higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]
             })}
  )
  
  return(higherTeirTree$node.label)
}

#####checkUnlabelledNodes: additional function to put a label on unlabeled nodes#####
checkUnlabelledNodes<-function(nodeLabel,defaultHigherStatus ){ifelse(grepl("/",nodeLabel),nodeLabel
                                                                      , ifelse(nodeLabel=="","",paste(nodeLabel,"/" , defaultHigherStatus, sep="")))}


#####transferAllNodeSupportAndStatus: Looping the above over all nodes and support values and add status#####
transferAllNodeSupportAndStatus<-function(lowerTeirTree, higherTeirTree,defaultLowerStatus,defaultHigherStatus, supportCutOff) {
  for (i in (Ntip(lowerTeirTree)):(Ntip(lowerTeirTree)+lowerTeirTree$Nnode))
  {higherTeirTree$node.label<-transferSingleNodeSupportAndStatus(lowerTeirTree, i, higherTeirTree, defaultLowerStatus, supportCutOff)
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
doesLowerCladeHaveAValidSupportLevel<-function(allTierTrees, lowerTierPart, leftoverTaxa, supportCutOff=95){
  binary<-try(
  {
    nodeLabelLower<-allTierTrees[[lowerTierPart]]$node.label[[ getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-Ntip(allTierTrees[[lowerTierPart]])]]
  as.numeric(nodeLabelLower)>=supportCutOff
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
  
  #browser()
  answer<-switch(tier, 
         higher = {
           
           newNodeNum=getMRCA(allTierTrees[[lowerTierPart-1]], 
                          intersect(unlist(getDescendantTipNames(allTierTrees[[1]], nodeNum)), allTierTrees[[lowerTierPart-1]]$tip.label)) #added this so that we can get the node number in the tier we are using ###This line might be causing problems 
           
           value=paste(
           allTierTrees[[lowerTierPart-1]]$node.label[[newNodeNum-Ntip(allTierTrees[[lowerTierPart-1]])]]
           ,"/",tierLabelHigher,sep="")
           value},#changed 1 to lowerTierPart-1 here because it shouldn't default to the highest tier. It should go with the next highest
         
         lower = paste(
           
           allTierTrees[[lowerTierPart]]$node.label[[ #the node label we're going to take
             ifelse(
               getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-Ntip(allTierTrees[[lowerTierPart]])==1, 1, getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-Ntip(allTierTrees[[lowerTierPart]]))#NOTE: SHOULD THIS be -1?
             
             ]],"/",tierLabelLower
           
           ,sep="")
         
         )#end of switch
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
  numTaxaHigh<-Ntip(highestTeirTree)
  sisterList<-list(1:numTaxaHigh-1) #create a dummy list to put the values into (for tippiest nodes)
  tipCladeMRCAS<-list(1:numTaxaHigh-1)#create a dummy list to put the values into (for tip-Clade pairs)
  
  for(i in 2:numTaxaHigh){
    #find the sister taxon to each tip in the tree
    sister<-
      getSisters(highestTeirTree,highestTeirTree$tip.label[[i]])
    
    #store the node ID of the MRCA for each tip-tip pair. For each tip-clade pair, store NA, which we will remove later
    sisterList[[i-1]]<-ifelse(all(sister <=numTaxaHigh),getMRCA(highestTeirTree,
                                                          c(highestTeirTree$tip.label[[i]],     
                                                            highestTeirTree$tip.label[[sister[[1]]  ]])), NA)
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


####This section (bipartition check) was added Dec 2022 and then functionalized Feb 2023
#determine if the current node lacks clade taxa in both halves of the descending bipartition. If TRUE, then the tier is an invalid source of node support

areCladeTaxaLackingFromAnyBiPartitions<-function(allTierTrees, currentNode, lowerTierPart){
    highestTeirTree<-allTierTrees[[1]]
    
    currentBipartition<-Children(highestTeirTree, currentNode)
    
    answer0<-c(1:length(currentBipartition))
    
    for(child in 1:length(currentBipartition)){
      
      currentChildTaxa<-getDescendantTipNames(highestTeirTree,currentBipartition[[child]])
      leftoverTaxa0<-intersect(unlist(currentChildTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
      
      answer0[child]<-length(intersect(unlist(currentChildTaxa),leftoverTaxa0 ))<1 #If the answer is "TRUE/1" then the taxa are all missing from the lower tier in one of the child clades
    }
    
    
    
    allAnswers0<-union(answer0, answer0)
    finalAnswer<-ifelse(
      (length(allAnswers0)>1#If it's greater than 1 then at least one of them is true and the other is false
       ||allAnswers0[[1]])
      , TRUE, FALSE)
    
    return(finalAnswer)
  
}


####This section lowestTierCladeAppearsIn was added feb 2023. 
####This function finds the lowest tier appropriate to get the node support values for the clade in question
####It checks both monophyly of the clade at all tiers and the bipartition presence of taxa within the clade


lowestTierCladeAppearsIn<-function(allTierTrees,currentNode, currentCladeTaxa){
  
    clTTMfun<-function(x){checkLowerTierTaxaMonophyly(allTierTrees, currentCladeTaxa, x)}# shortened version of checkLowerTierTaxaMonophyly that works in sapply easily
    aCTLFABPfun<-function(x){areCladeTaxaLackingFromAnyBiPartitions(allTierTrees, currentNode, x)}
    
    presentInTierAns<-sapply(2:length(allTierTrees), clTTMfun)
    incompleteInTierAns<-sapply(2:length(allTierTrees), aCTLFABPfun)
    
    
    combinedAnswer<-as.logical(((!incompleteInTierAns)*presentInTierAns))
    
    pie<-c()
    i<-0
    #ans<-4
    for(ans in length(allTierTrees):2){
      i<-i+1
      pie[[i]]<-ifelse(combinedAnswer[[ans-1]],ans, 1)##The FALSE condition used to be 0 here but I don't think that's right. It should be 1 because this is going to  correspond to the lowertierPart. However, if I find it needs to be 0 then I should make an extra conditional outside of this function where I deal with the 0
    }
    allPies<<-rev(sort(unlist(pie)))
    return (max(unlist(pie)))
}

findNodeID<-function(phy, taxonList){findMRCA(phy, taxonList)-Ntip(phy)}
findNodeLabel<-function(phy, taxonList){phy$node.label[[findMRCA(phy, taxonList)-Ntip(phy)]]}


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#######for a single node, find which tree has the valid node support value#####
###############################################################################

##NOTE: This is v2 of this function. V1 is in a different file.

evaluateNodeSupport<-function(
    allTierTrees, #this object should have the trees in order from highest tier to lowest tier (i.e. most inclusive to least inclusive)
    tierLabels, #the names of the three trees in the same order as the trees
    nodeNumber, #the number of the higher tier node we are evaluating
    nodeIterator,  #an arbitrary sequential iterator used to keep track of where to store the labeled node IDs
    supportCutOff = 95,
    verbose = FALSE,
    developmentMode=FALSE
){
  #browser()
  
  if(verbose==TRUE) {print(paste("Starting on highest tier tree, node ", as.character(nodeNumber)))}
  
  highestTeirTree<-allTierTrees[[1]]
  numTaxaHigh<-Ntip(highestTeirTree)
  
  #Read a single node
  
  #descTipIDs<-Descendants(highestTeirTree, nodeNumber, "tips")
  descTipNames<-getDescendantTipNames(highestTeirTree, nodeNumber)
  
  #Quality check 1
  if(developmentMode==TRUE) {print(ifelse(
    {length(setdiff(as.list(as.data.frame(highestTeirTree$tip.label)), as.list(descTipNames)))==0}, "Error: This is the root", "Print:QC1 passed" ))
  }
  
  currentNode<-nodeNumber
  currentCladeTaxa<-getDescendantTipNames(highestTeirTree,currentNode)
  
  lowerTierPart<-lowestTierCladeAppearsIn(allTierTrees,currentNode, currentCladeTaxa) #this was originally defined as "2" and then allowed to be adjusted using the workflow below. However, I think this is unnecessary.
  
  
  leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label) #this is a list of all the taxa that are in the clade of interest in the lowered tier
  if(developmentMode==TRUE){print(paste("Lower teir part set to initial value: ", as.character(lowerTierPart)))}
  
  cladeLabeled<-FALSE #set this to FALSE so that when this is evaluated I can exit the loop when the clade is correctly labeled
  
  
  while(cladeLabeled==FALSE){
    
    
    ###Check if this node label has been used before
    if(is.element(getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa)), finishedNodes[lowerTierPart][[1]]
                  #[!is.na(  finishedNodes[lowerTierPart][[1]])]
    )
    ){ 
      ###If true, go to the next valid node and repeat
      allPies<-allPies[allPies != lowerTierPart]
      lowerTierPart<-ifelse(allPies[[1]]==0, 1, allPies[[1]])
      
      if(developmentMode==TRUE){print(paste("Now evaluating tier: ", as.character(lowerTierPart)))}
      
      
    }else{
      ###If false, then label the node
      { if(developmentMode==TRUE) {print("Target value found...using label on node")} 
        
        #record that we're labeling this node
        finishedNodes[[lowerTierPart]][[nodeIterator]]<-getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa))
        
        #old method, removed feb 7 2023
        #highestTeirTree$node.label[[currentNode-Ntip(highestTeirTree)]]<-useWhichNodeLabel("lower", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa)
        
        #browser()
        
        ###7 Feb 2023: There's a problem I don't quite understand, where the node IDs for the lowertier trees and the highest tier trees are not being indexed in the same way. So, I will use a workaround where, if lowerTierPart==1 I will simply not change the node value, but add on the label to it. While this doesn't solve the fundamental source of the problem, it shouldn't cause an issue with other datasets. So, it may (or may not) be a general solution to the bug
        
        if(lowerTierPart>1){
        highestTeirTree$node.label[[findNodeID(highestTeirTree, leftoverTaxa)]]<-paste(
          allTierTrees[[lowerTierPart]]$node.label[[
            findNodeID(allTierTrees[[lowerTierPart]], leftoverTaxa)#this part of the finished nodes list, which was just defined above, should have the nodeID for the label we want to use
          ]], tierLabels[[lowerTierPart]], sep="/")}else{
                                          
            highestTeirTree$node.label[[findNodeID(highestTeirTree, leftoverTaxa)]]<- paste(highestTeirTree$node.label[[findNodeID(highestTeirTree, leftoverTaxa)]], tierLabels[[1]], sep="/")
            }
        
        
        
        
        
        
        if(verbose==TRUE) {print(
          paste(
            paste(
              paste("Using tier ",as.character(lowerTierPart)," for label on node")
              , as.character(currentNode)) #added print line feb 6 2023 to track proper node labelling
                                 ,as.character(highestTeirTree$node.label[[currentNode-Ntip(highestTeirTree)]]), sep = "->"))} 
        
        cladeLabeled<-TRUE
        raiseNode<-FALSE #exit the loop of raising the tier level
      }#End of QC Check 3 part
    }###end of check if label has been used before
    
    
    
  }
  output<-unlist(
    c(
      highestTeirTree$node.label[[currentNode-Ntip(highestTeirTree)]],
      unlist(lowerTierPart[[1]]), 
      
      {
        mrca<-getMRCA(allTierTrees[[lowerTierPart]],intersect(leftoverTaxa, allTierTrees[[lowerTierPart]]$tip.label))
        ifelse(is.null(mrca), NA, mrca)
      }
    ))
  #output
  return(output)
}







#Determine the order of nodes to evaluate on the higher tree
findOrderToTraverseHighestTree<-function(allTierTrees){
  findTippiestClades(allTierTrees[[1]])
  
  highestTierNodesOrdered<-unlist(as.list(matrix(data = NA, nrow =((Ntip(allTierTrees[[1]])+Nnode(allTierTrees[[1]]))-(Ntip(allTierTrees[[1]]))), ncol=1 )))
  
  #add the root to the list
  highestTierNodesOrdered[[1]]<-Ntip(allTierTrees[[1]])+1
  
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



#allTierTrees<-allRandomizedTierTrees

#Final wrapper function that evaluates the whole tree.
totalAwareSupport<-function(allTierTrees, tierNames, supportCutOff = 95, verbose=FALSE, developmentMode=FALSE){
  highestTierTree0<<-allTierTrees[[1]]
  nodeList<-findOrderToTraverseHighestTree(allTierTrees) #find the order to traverse the highest (most inclusive) tree
  nodeIterator<-0
  
  finishedNodes<<-data.frame(matrix(NA,nrow = length(nodeList),ncol = length(allTierTrees))) #this is a high level variable that will be used to keep track of which nodes have already been tested

    for(j in nodeList[-1]){#This is omitting the first element because the first element is the root. 
         
    nodeIterator<-nodeIterator+1
    ifelse(is.element(j, finishedNodes[[1]]), print("Error: Annotating higher tier tree node twice"), if(developmentMode==TRUE) {print("QC Check 2: Good")} else print("...") )#this checks the highest tier tree (the one we are putting support values onto) and checking that it hasn't already been annotated
    
    answer<-evaluateNodeSupport(allTierTrees, tierNames, j, nodeIterator, supportCutOff, verbose)
    #browser()
    
    if(developmentMode==TRUE) {print(paste("Replacing ", as.character(highestTierTree0$node.label[[j-Ntip(highestTeirTree)]]), " with ", as.character(answer[[1]]), sep = ""))}
    
    highestTierTree0$node.label[[j-Ntip(highestTeirTree)]]<-answer[[1]]
    finishedNodes[[1]][[nodeIterator]]<-j #this keeps track of the highest tier nodes annotated
    finishedNodes[[as.numeric(answer[[2]])  ]][[nodeIterator]]<-answer[[3]] #this keeps track which lower tier's nodes have been used

  }
  return(list(highestTierTree0, finishedNodes))
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



#####Generate some random test trees

listNodes<-function(phy){
  listOfTips<-phy$tip.label
  
  allNodes<-data.frame(matrix(NA,nrow = 1,ncol = (length(listOfTips)+phy$Nnode) ))
  
  for(tip in 1:(length(listOfTips)+phy$Nnode)){
    allNodes[[tip]]<-nodeMRCA(phy, tip)
  }
  return(union(unlist(allNodes),unlist(allNodes)))
}

range01 <- function(x){(x-min(x))/(max(x)-min(x))}#make a function to scale a list between 0 and 1



#make a randomized tree from scratch, with randomized bootstrap values
makeRandomizedTreeWithRandSupport<-function(nTips, nTiers){
  
  BackBoneTips<-nTips
  
  #generate topology
  tierBackbone<-multi2di(root(
    rtree(BackBoneTips, rooted = FALSE, tip.label=1:BackBoneTips, br=1)
    ,outgroup="1", resolve.root=TRUE)
    )
  
  #generate node labels
  allNodes<-listNodes(tierBackbone)
  
  possibleNodeValues<-round((range01(rgamma(1000, shape = 2, rate = 2))*100-100)*-1 ,0)#draw random values from a gamma distribution and scale them between 0 and 100, with the highest probability draws being near 100. 
  
  for(nodes in 1:(length(allNodes)) ){
    tierBackbone$node.label[[nodes]]<-sample(possibleNodeValues, 1)#take some of those randomly drawn node labels and put them on the tree
  }
  tierBackbone$edge.length<-replace(tierBackbone$edge.length, tierBackbone$edge.length == 0, 1)
  
  #replace non-root tip names with a numeric code
  tierBackbone$tip.label<- replace(tierBackbone$tip.label,  tierBackbone$tip.label>1, seq(10^nTiers, nTips*(10^nTiers)-1, 10^nTiers))
  
  tierBackbone$node.label[[1]]<-100
  
  plot.phylo(tierBackbone, show.node.label = TRUE)
  return(tierBackbone)
  
}



#A function that collapses unsupported nodes and then randomly resolves them again
collapseAndReResolveRandomly<-function(phy, cutoff = 75){
  
  apple<-collapseUnsupportedEdges(phy, value = "node.label", cutoff)
  #apple<-CollapseNode(phy, length(phy$tip.label)+match(TRUE, phy$node.label<cutoff))
  
  pear<-multi2di(apple)
  
  possibleNodeValues<-round((range01(rgamma(1000, shape = 2, rate = 2))*100-100)*-1 ,0)#draw random values from a gamma distribution and scale them between 0 and 100, with the highest probability draws being near 100.
  

  pear$node.label<-replace(pear$node.label,pear$node.label=="", sample(possibleNodeValues, sum(pear$node.label=="")      ))
  
  pear$edge.length<-rep.int(1, length(pear$edge.length))
  
  
  return(pear)
}


#Now make a function that randomly adds tips to an existing tree and choosing a random node label for it

addRandomTaxa<-function(phy, nTips, rootingTaxon = "1", cutoff = 75){
###NOTE: nTips can be no larger than the number of tips in phy
  
  
  previousTipNames<-as.numeric(phy$tip.label)
  
  
  a<-floor(runif(100, min=2, max=(nTips+Ntip(phy))))
  newTipLabels<-setdiff(a,phy$tip.label )
  
  newTipNames<-previousTipNames
  newTipList<-sample(newTipLabels, nTips)
  
  #first, randomly put tips onto the tree
  newTree<-phy
  
  for(n in 1:nTips){
  possibleSpots<-2:(Ntip(newTree)+newTree$Nnode)
  newTree<-add.tips(newTree,newTipList[[ n ]], sample(possibleSpots, 1), edge.length = 1)
  }
  newTree<-root.phylo(newTree, rootingTaxon, resolve.root = TRUE)
  
    #generate node labels
    
  possibleNodeValues<-round((range01(rgamma(1000, shape = 2, rate = 2))*100-100)*-1 ,0)#draw random values from a gamma distribution and scale them between 0 and 100, with the highest probability draws being near 100. 
    
    newTree2<-multi2di(newTree)
    newTree2$edge.length<-replace(newTree2$edge.length, newTree2$edge.length == 0, 1)
    allNodes<-listNodes(newTree2)
    
    newTree2$node.label[sapply(newTree2$node.label, is.null)] <- ""
    newTree2$node.label[sapply(newTree2$node.label, is.na)] <- ""
    newTree2$node.label<-replace(newTree2$node.label,newTree2$node.label=="", sample(possibleNodeValues, sum(newTree2$node.label=="")))
    
    newTree2$edge.length<-replace(newTree2$edge.length, newTree2$edge.length == NA, 1)
    
    
#    for(nodes in 1:(length(newTree2$node.label)) ){
#      newTree2$node.label[[nodes]]<-sample(possibleNodeValues, 1)#take some of those randomly drawn node labels and put them on the tree
#    }
    
    #newTree2<-multi2di(newTree2)
    plot.phylo(newTree2, show.node.label = TRUE)
    return(newTree2)
    
}

    
#now combine the randomized tree functions into a single one for more easily making randomized test data NOTE: This will only make a 4 tiered scenario

makeRandomizedTiers<-function(nTips, nTaxa2Add, cutoff){
  
  #Backbone tier
  randStartPhy<-makeRandomizedTreeWithRandSupport(nTips, 4)
  randStartPhy2<-collapseAndReResolveRandomly(randStartPhy, cutoff)

  
  #Skeleton tier
  tier2Phy<-addRandomTaxa(randStartPhy2, nTaxa2Add, rootingTaxon = "1", cutoff)
  tier2Phy2<-collapseAndReResolveRandomly(tier2Phy, cutoff)
  
  #Flesh tier
  tier3Phy<-addRandomTaxa(tier2Phy2, nTaxa2Add, rootingTaxon = "1", cutoff)
  tier3Phy2<-collapseAndReResolveRandomly(tier3Phy, cutoff)
  
  #Whole tier
  tier4Phy<-addRandomTaxa(tier3Phy2, nTaxa2Add, rootingTaxon = "1", cutoff)
  
  
  #combine
  allRandomizedTierTrees<-c(tier4Phy, tier3Phy, tier2Phy, randStartPhy)
  
  plot.phylo(tier4Phy, show.node.label = TRUE)
  
  
  return(allRandomizedTierTrees)
}









