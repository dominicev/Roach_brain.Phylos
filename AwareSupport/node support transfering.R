library(ape)
library(phangorn)
library(phytools)
library(adephylo)
library(TreeDist)
library(dplyr)
library(TreeTools)
library(RRphylo)
devtools::source_url("https://github.com/dominicev/Roach_brain_Phylos/blob/main/tree%20functions.R")
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
  nodeLabel<-lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]
  #find same node in next tree
  higherTreeCorrespondingNode<-getMRCA(higherTeirTree, unlist(as.data.frame(descTipNames)))
  #compare the node support values between the two trees
  lowerTierNodeValue<-lowerTeirTree$node.label[nodeNumber-length(lowerTeirTree$tip.label)]
  higherTierNodeValue<-higherTeirTree$node.label[higherTreeCorrespondingNode-length(higherTeirTree$tip.label)]
  
  #replace the node values if the value is >cutOff%
  ifelse(
    as.numeric(lowerTierNodeValue)>=supportCutOff, 
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


transferAllNodeSupportValues<-function(lowerTeirTree, higherTeirTree, supportCutOff = 95) {
  for (i in (length(lowerTeirTree$tip.label)):(length(lowerTeirTree$tip.label)+lowerTeirTree$Nnode))
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
  
  #replace the node values if the value is >supportCutOff%
  ifelse(
    as.numeric(lowerTierNodeValue)>=supportCutOff, 
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
transferAllNodeSupportAndStatus<-function(lowerTeirTree, higherTeirTree,defaultLowerStatus,defaultHigherStatus, supportCutOff) {
  for (i in (length(lowerTeirTree$tip.label)):(length(lowerTeirTree$tip.label)+lowerTeirTree$Nnode))
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
    nodeLabelLower<-allTierTrees[[lowerTierPart]]$node.label[[ getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-length(allTierTrees[[lowerTierPart]]$tip.label)]]
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
#

useWhichNodeLabel<-function(tier, allTierTrees, lowerTierPart, nodeNum, tierLabelHigher, tierLabelLower, leftoverTaxa){
  #currentCladeTaxa<-getDescendantTipNames(allTierTrees[[1]],nodeNum)
  
  answer<-switch(tier, 
         higher = {
           #browser()
           newNodeNum=getMRCA(allTierTrees[[lowerTierPart-1]], 
                          intersect(unlist(getDescendantTipNames(allTierTrees[[1]], nodeNum)), allTierTrees[[lowerTierPart-1]]$tip.label)) #added this so that we can get the node number in the tier we are using ###This line might be causing problems 
           
           value=paste(
           allTierTrees[[lowerTierPart-1]]$node.label[[newNodeNum-length(allTierTrees[[lowerTierPart-1]]$tip.label)]]
           ,"/",tierLabelHigher,sep="")
           value},#changed 1 to lowerTierPart-1 here because it shouldn't default to the highest tier. It should go with the next highest
         lower = paste(allTierTrees[[lowerTierPart]]$node.label[[ getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-length(allTierTrees[[lowerTierPart]]$tip.label)]],"/",tierLabelLower,sep=""))
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





###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
#######for a single node, find which tree has the valid node support value#####
###############################################################################

evaluateNodeSupport<-function(
                              allTierTrees, #this object should have the trees in order from highest tier to lowest tier (i.e. most inclusive to least inclusive)
                              tierLabels, #the names of the three trees in the same order as the trees
                              nodeNumber, #the number of the higher tier node we are evaluating
                              nodeIterator,  #an arbitrary sequential iterator used to keep track of where to store the labeled node IDs
                              supportCutOff = 95,
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
  leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)#this is a list of all the taxa that are in the clade of interest in the lowered tier
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
    
    ####This section (bipartition check) was added Dec 2022.
    
      currentBipartition<-Children(highestTeirTree, currentNode)
      #child<-2
      
      answer0<-c(1:length(currentBipartition))
        
      for(child in 1:length(currentBipartition)){
        
        currentChildTaxa<-getDescendantTipNames(highestTeirTree,currentBipartition[[child]])
        leftoverTaxa0<-intersect(unlist(currentChildTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
        
        answer0[child]<-length(intersect(unlist(currentChildTaxa),leftoverTaxa0 ))<1 #If the answer is "TRUE/1" then the taxa are all missing from the lower tier in one of the child clades
        }
      allAnswers0<-union(answer0, answer0)
      areCladeTaxaLackingFromAnyBiPartitions<-ifelse(
        (length(allAnswers0)>1#If it's greater than 1 then at least one of them is true and the other is false
        ||allAnswers0[[1]])
        , TRUE, FALSE)#since TRUE/1 means the taxa are missing from the bipartition, then BOTH bipartitions must be FALSE in order for the node on the lower tier to be valid
        
      
    #doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart, leftoverTaxa, supportCutOff)##DOES THIS NEED TO BE HERE?
    
    
    always<-"NEVER"#This object should never change in value. This is only to allow the while chunk below to run properly. It will be exited when the outer loop is exited, or when a "break" occurs.
    while(always=="NEVER"){
      
          ifelse(areCladeTaxaAllHigher||areCladeTaxaLackingFromAnyBiPartitions
           #Step 1Ai-iii
           ,#if all the taxa are higher tier taxa
           {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa)
           if(verbose==TRUE){print("Use higher tier label (cond. 1)")}
           #browser()
           cladeLabeled<-TRUE
           }, #if NOT... 
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
                             check<-doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart,leftoverTaxa, supportCutOff )
                             
                           check},{raiseNode<-TRUE
                             while(raiseNode==TRUE){
                                       ifelse(
                               is.element(getMRCA(allTierTrees[[lowerTierPart]],unlist(leftoverTaxa)), finishedNodes[lowerTierPart][[1]][!is.na(  finishedNodes[lowerTierPart][[1]])] ) #check if the node on the lower tier tree has been used in a label already
                                    ,{ print("Label used previously. Moving up one tier and rechecking.")         #if so, change the tier to a higher one
                                      lowerTierPart<-lowerTierPart-1
                                      leftoverTaxa<-intersect(unlist(currentCladeTaxa),  allTierTrees[[lowerTierPart]]$tip.label)
                                      ifelse(doesLowerCladeHaveAValidSupportLevel(allTierTrees, lowerTierPart, leftoverTaxa, supportCutOff)||(lowerTierPart==1), 
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
                                    highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("lower", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa)
                                    print(paste("Using",as.character(lowerTierPart),"tier label"))
                                    cladeLabeled<-TRUE
                                    raiseNode<-FALSE #exit the loop of raising the tier level
                                    })
                               print("QC Check 4: Good")
                             }#END OF while raiseNode==TRUE
                           print("QC Check 5: Good")
                             },
                                  
                                  {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa)#leftovertaxa to leftoverTaxa typo here and tierLabels[[1]] to tierLabels[[lowerTierPart-1]] Dec 28 2022
                                  print("Use higher tier label (cond. 3)")
                                  cladeLabeled<-TRUE}
                           )}
                    )#label 4: if true
                  ,
                         {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa) # tierLabels[[1]] to tierLabels[[lowerTierPart-1]] Dec 28 2022
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
totalAwareSupport<-function(allTierTrees, tierNames, supportCutOff = 95, verbose=FALSE){
  highestTierTree<<-allTierTrees[[1]]
  nodeList<-findOrderToTraverseHighestTree(allTierTrees) #find the order to traverse the highest (most inclusive) tree
  nodeIterator<-0
  
  finishedNodes<<-data.frame(matrix(NA,nrow = length(nodeList),ncol = length(allTierTrees))) #this is a high level variable that will be used to keep track of which nodes have already been tested

    for(j in nodeList[-1]){ #This is omitting the first element because the first element is the root. 
    nodeIterator<-nodeIterator+1
    ifelse(is.element(j, finishedNodes[[1]]), print("Error: Annotating higher tier tree node twice"), print("QC Check 2: Good."))#this checks the highest tier tree (the one we are putting support values onto) and checking that it hasn't already been annotated
    
    answer<-evaluateNodeSupport(allTierTrees, tierNames, j, nodeIterator, supportCutOff, verbose)
    
    highestTierTree$node.label[[j-length(highestTierTree$tip.label)]]<-answer[[1]]
    finishedNodes[[1]][[nodeIterator]]<-j #this keeps track of the highest tier nodes annotated
    finishedNodes[[as.numeric(answer[[2]])  ]][[nodeIterator]]<-answer[[3]] #this keeps track which lower tier's nodes have been used

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
           finalTree$node.label[[nodeID-length(finalTree$node.label)]]<-paste(finalTree$node.label[[nodeID-length(finalTree$node.label)]], #original node support value
                                                                                sourceTree$node.label[[correspondingNodeID-length(sourceTree$tip.label)]], sep="|") #new node support value
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
makeRandomizedTreeWithRandSupport<-function(nTips){
  
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
  
  plot.phylo(tierBackbone, show.node.label = TRUE)
  return(tierBackbone)
  
  
  }

#A function that collapses unsupported nodes and then randomly resolves them again

collapseAndReResolveRandomly<-function(phy, cutoff = 75){
  apple<-collapseUnsupportedEdges(phy, value = "node.label", cutoff)
  pear<-multi2di(apple, random = TRUE, equiprob = TRUE)
  
  possibleNodeValues<-round((range01(rgamma(1000, shape = 2, rate = 2))*100-100)*-1 ,0)#draw random values from a gamma distribution and scale them between 0 and 100, with the highest probability draws being near 100. 
  
  pear$node.label<-replace(pear$node.label,pear$node.label=="", sample(possibleNodeValues, pear$Nnode      ))
  
  return(pear)
}



#Now make a function that randomly adds tips to an existing tree and choosing a random node label for it

addRandomTaxa<-function(phy, nTips, factor, rootingTaxon = "1"){
###NOTE: nTips can be no larger than the number of tips in phy
  
  
  previousTipNames<-as.numeric(phy$tip.label)
  newTipNames<-previousTipNames*factor
  newTipList<-sample(newTipNames, length(newTipNames))
  
  #first, randomly put tips onto the tree
  newTree<-phy
  
  for(n in 1:nTips){
  possibleSpots<-2:(length(newTree$tip.label)+newTree$Nnode)
  newTree<-add.tips(newTree,newTipList[[ n ]], sample(possibleSpots, 1), edge.length = 1)
  }
  newTree<-root.phylo(newTree, rootingTaxon, resolve.root = TRUE)
  
    #generate node labels
    
    
    possibleNodeValues<-round((range01(rgamma(1000, shape = 2, rate = 2))*100-100)*-1 ,0)#draw random values from a gamma distribution and scale them between 0 and 100, with the highest probability draws being near 100. 
    
    newTree2<-multi2di(newTree)
    newTree2$edge.length<-replace(newTree2$edge.length, newTree2$edge.length == 0, 1)
    allNodes<-listNodes(newTree2)
    for(nodes in 1:(length(newTree2$node.label)) ){
      newTree2$node.label[[nodes]]<-sample(possibleNodeValues, 1)#take some of those randomly drawn node labels and put them on the tree
    }
    
    #newTree2<-multi2di(newTree2)
    plot.phylo(newTree2, show.node.label = TRUE)
    return(newTree2)
    
}

    
#now combine the randomized tree functions into a single one for more easily making randomized test data NOTE: This will only make a 4 tiered scenario

makeRandomizedTiers<-function(nTips, nTaxa2Add, cutoff ){
  
  #Backbone tier
  randStartPhy<-makeRandomizedTreeWithRandSupport(nTips)
  randStartPhy2<-collapseAndReResolveRandomly(randStartPhy, cutoff)
  
  #Skeleton tier
  tier2Phy<-addRandomTaxa(randStartPhy2, nTaxa2Add, 100, rootingTaxon = "1")
  tier2Phy2<-collapseAndReResolveRandomly(tier2Phy, cutoff)
  
  #Flesh tier
  tier3Phy<-addRandomTaxa(tier2Phy2, nTaxa2Add, 100, rootingTaxon = "1")
  tier3Phy2<-collapseAndReResolveRandomly(tier3Phy, cutoff)
  
  #Whole tier
  tier4Phy<-addRandomTaxa(tier3Phy2, nTaxa2Add, 100, rootingTaxon = "1")
  
  
  #combine
  allRandomizedTierTrees<-c(tier4Phy, tier3Phy, tier2Phy, randStartPhy)
  
  plot.phylo(tier4Phy, show.node.label = TRUE)
  
  
  return(allRandomizedTierTrees)
}









