
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
    higherTierTaxa<-setdiff(highestTeirTree$tip.label,allTierTrees[[lowerTierPart]]$tip.label) #all higher tier taxa in the WHOLE tree
    
    ##precalculate some conditionals to make the below code a little easier to read.
    areCladeTaxaAllHigher<-length(
      unlist(currentCladeTaxa) #how many total tips are in the clade considered
    )-(length(
      intersect(unlist(currentCladeTaxa), higherTierTaxa) #how many of those tips are higher tier taxa
    ))<=1
    
    ####This section (bipartition check) was added Dec 2022.
    
    
    aCTLFABP<-areCladeTaxaLackingFromAnyBiPartitions(allTierTrees, currentNode, lowerTierPart)#since TRUE/1 means the taxa are missing from the bipartition, then BOTH bipartitions must be FALSE in order for the node on the lower tier to be valid
    
    
    
    
    always<-"NEVER"#This object should never change in value. This is only to allow the while chunk below to run properly. It will be exited when the outer loop is exited, or when a "break" occurs.
    while(always=="NEVER"){
      
      ifelse(areCladeTaxaAllHigher||aCTLFABP
             #Step 1Ai-iii
             ,#if all the taxa are higher tier taxa
             {highestTeirTree$node.label[[currentNode-length(highestTeirTree$tip.label)]]<-useWhichNodeLabel("higher", allTierTrees, lowerTierPart, currentNode, tierLabels[[lowerTierPart-1]], tierLabels[[lowerTierPart]], leftoverTaxa)
             if(verbose==TRUE){print("Use higher tier label (cond. 1)")}
             #browser()
             cladeLabeled<-TRUE
             }, #if NOT... 
             #
             ifelse(checkLowerTierTaxaMonophyly(allTierTrees, currentCladeTaxa, lowerTierPart), #label: 4
                    
                    ###add something here
                    ###we need to check the earliest tier in which the clade arose and use that support value
                    
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











