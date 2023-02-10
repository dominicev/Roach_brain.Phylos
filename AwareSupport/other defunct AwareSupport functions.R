
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


#####doesLowerCladeHaveAValidSupportLevel: returns T/F after checking if the clade on the lower tier is above the support cutoff
doesLowerCladeHaveAValidSupportLevel<-function(allTierTrees, lowerTierPart, leftoverTaxa, supportCutOff=95){
  binary<-try(
    {
      nodeLabelLower<-allTierTrees[[lowerTierPart]]$node.label[[ getMRCA(allTierTrees[[lowerTierPart]],leftoverTaxa)-Ntip(allTierTrees[[lowerTierPart]])]]
      as.numeric(nodeLabelLower)>=supportCutOff
    }
    , silent=TRUE)
  return(isTRUE(binary))
}




###########################################################################################################
#####NOTE##################################################################################################
#The above functions go in order from the most conservative (most exclusive) tree to the largest (least exclusive) tree. However, it appears that many of the nodes are incorrectly labeled. A good example of this is Blaberoidea. Blaberoidea is supported in the ABA tree with 100% support but there is no Blaberoidea (SS or SL) node with 100% support (check CDD tree with modified support).
#Below I will try and accomplish the same goal going in the reverse order (starting from the largest tree and working backwards) to see if it doesn't result in a more expected result.
###########################################################################################################