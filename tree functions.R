
library(ape)
library(phangorn)
library(phytools)
library(adephylo)
library(fitdistrplus)
library(reshape2)
library(dplyr)
library(rlist)
library(MonoPhy)
library(data.table)
library(foreach)
library(doParallel)
library(TreeTools)

#source("B:\OneDrive - University of Illinois - Urbana\Science\R programs\Roach_brain_Phylos\tree functions.R")
 #use this to load this file into an r book

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

###removing taxa that have extremely large pairwise distances
###
###This will compare the pairwise tip distances in the whole tree to those generated under an estimated gamma distribution and remove those that appear to have numerous cases (relative to the total number of taxa) where they are outside the alpha (default 95%) quantile of that distribution.
###
#troubelshooting
#setwd("G:\\My Drive\\Projects\\PHY_Blab_2\\bioinformatics\\PhyBlab2 full analysis\\loci\\3 - gene trees\\Full gene tree analysis 28-12-21")


#treeFileNames<-list.files(pattern=".treefile")
#alpha = 0.95
#relAbund = "Default"
#randomSeed = 2345


trimPairTipOutliers<-function(tree, alpha = 0.95, relAbund = "Default", randomSeed = 2345 ) {
      ##Import tree
      phylo<-tree
      #plot.phylo(phylo)
      
      #Define the number of times the taxon name has to appear in list of long branch taxa to be considered an true long branch taxa
      if(relAbund=="Default"){relAbund<-(length(phylo[[4]])/4)}
      
      #calculate root to tip branch length distances
      dists<-distTips(phylo, tips = "all", method = "patristic")
      
      #reshape dists into a list
      df <- melt(as.matrix(dists), varnames = c("row", "col"))
      allDists<-df$value
      #hist(allDists)
      
      #fit gamma distribution to distances
      set.seed(randomSeed)
      fit.gamma = fitdist(allDists[allDists>0], distr = "gamma", method = "mle")
      #fit.gamma[1]
      
      #Calculate outliers based on random sampling of gamma distribution
      randomGammaValues<-rgamma(9999, fit.gamma$estimate[[1]], rate = fit.gamma$estimate[[2]] )
      
      #shortBranches<-quantile(randomGammaValues, c(.001))
      longBranches<-quantile(randomGammaValues, c(alpha))
      
      rowNames<-df$row; colNames<-df$col
      longBranchTaxa<-unlist(c(rowNames[df$value>longBranches],colNames[df$value>longBranches]))
      frequencies<-sort(table(longBranchTaxa), decreasing = TRUE)
      strictLBTaxa<-dimnames(frequencies[frequencies>relAbund])
      #delete outlier tips and export
      
      ifelse(length(strictLBTaxa)>0,phylo<-drop.tip(  phylo,strictLBTaxa$longBranchTaxa ), 0)#redefines the phylogeny
      
      #write.tree(phylo, paste(strsplit(files, ".", fixed=TRUE)[[1]][[1]],".tiptrimmed.tre", sep=""))
      print(paste(length(strictLBTaxa$longBranchTaxa), "tips dropped from tree"))
      
      return(phylo)
      
      
    
}

###removing taxa that have extremely large pairwise distances
###
###This will compare the pairwise tip distances in the whole tree to those generated under an estimated gamma distribution and remove those that appear to have numerous cases (relative to the total number of taxa) where they are outside the alpha (default 95%) quantile of that distribution.
###


####this is the same as the above function but it operates on a single tree at a time

trimRootTipOutliers<-function(tree, alpha = 0.99, method = "Default", randomSeed = 2345) {

    ##Import tree
    phylo<-tree
    #plot.phylo(phylo)
    
    #root tree
    ifelse(method == "Default", phylo<-midpoint.root(phylo), phylo<-rootIncompleteTree(phylo, outgroupsOrdered))
    
    #calculate root to tip branch length distances
    dists=distRoot(phylo, tips = "all", method = "patristic")
    
    #fit gamma distribution to distances
    set.seed(randomSeed)
    fit.gamma = fitdist(dists, distr = "gamma", method = "mle")
    fit.gamma[1]
    
    #Calculate outliers based on random sampling of gamma distribution
    randomGammaValues<-rgamma(9999, fit.gamma$estimate[[1]], rate = fit.gamma$estimate[[2]] )
    
    #shortBranches<-quantile(randomGammaValues, c(.001))
    longBranches<-quantile(randomGammaValues, c(alpha))
    
    #delete outlier tips and export
    { 
      longBranchTaxa<-c()
      j=0;k=0
      for(value in dists){
        j=j+1;
        name=names(dists[j])
        if(value>longBranches[[1]]){k=k+1;
        longBranchTaxa[k]<-name}
      }
    }###gets the taxon names of outliers
    longBranchTaxa
    ifelse(length(longBranchTaxa)>0,phylo<-drop.tip(phylo,longBranchTaxa ), 0)#redefines the phylogeny
    
    #write.tree(phylo, paste(strsplit(files, ".", fixed=TRUE)[[1]][[1]],".tiptrimmed.tre", sep=""))
    print(paste(length(longBranchTaxa), "tips dropped from tree"))
    return(phylo)
  
  
}


###recursive tip trimming function


recursiveTrimOutlierTaxa<-function(tree){
  phylo<-tree
  repeat{
    previousTaxaLength<-length(phylo$tip.label)
    phylo<-trimRootTipOutliers(trimPairTipOutliers(phylo))
    currentTaxonLength<-length(phylo$tip.label)
    if(currentTaxonLength==previousTaxaLength){break()}
  }
  return(phylo)
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



######deltaGLS (gene likelihood scores) test######

#deltaGLS was used in Shen et al. (2017) to show that contentious relationships can be driven by outlier genes.
#These functions will allow calculation of these values from a set of genes given 2-4 alternative hypotheses

#input...
#clusterFileName is a string specifying the name of the cluster file. Unless it's in your working directory, you need to specify the whole directory+file name
#treeHypotheses must be a multi.phylo object containing your different topological hypotheses. It will only work if this has between 2 and 4 trees in it. Any more and you need to do multiple runs.
#genesDirectory the directory where your alignment files are. These must be in fasta format and have the file extension ".fasta"
#GTs optional


deltaGLS<-function(clusterFileName, treeHypotheses, genesDirectory, GTs = c()){
  #read cluster file, which must have 4 clusters and one list of taxa to ignore. See example file for format.
  allClusters<-read.delim(clusterFileName, header=FALSE, sep = "=")
  cluster1<<-c(strsplit(allClusters[1, 1], " +")[[1]][2], lapply(strsplit(allClusters[1, 2], "\t+"), trimws ) )
  cluster2<<-c(strsplit(allClusters[2, 1], " +")[[1]][2], lapply(strsplit(allClusters[2, 2], "\t+"), trimws ) )
  cluster3<<-c(strsplit(allClusters[3, 1], " +")[[1]][2], lapply(strsplit(allClusters[3, 2], "\t+"), trimws ) )
  cluster4<<-c(strsplit(allClusters[4, 1], " +")[[1]][2], lapply(strsplit(allClusters[4, 2], "\t+"), trimws ) )
  
  taxaToRemove<-c(strsplit(allClusters[5, 1], " +")[[1]][2], lapply(strsplit(allClusters[5, 2], "\t+"), trimws ) )
  
  #if no checkGeneValidity file is available then run that subscript. Otherwise, import the list of valid genes
  ifelse(file.exists("validGenes.csv"),{
    #if the file exists then import it
    print("Importing validGenes.csv")
    validGenes<-unlist(as.list(read.csv("validGenes.csv"))[[2]])
  },{
    #if not, then find the valid files...which takes sometime
    print("Finding which alignments are valid (this will take a while...particularly if you have more than 200 alignments)")
    validGenes<-foreach(i = 1:length(GTs), .combine = "c", .packages=c('ape', 'phangorn'))%do%{
      checkGeneValidity(GTs[[i]],genesDirectory )}
    write.csv(validGenes, file="validGenes.csv")
    
  })
  
  print(paste(
    "Your analysis will have", length(validGenes), "genes in it.") )
  
   print("WARNING: If it stops on one for any more than a few seconds there may be a problem. 
    The pml() function (calculates lnL in Phangorn) doesn't handle very complex alignments (i.e., long, with messy sequences, possible paralogy).")
  
  print("You may need to remove these from your dataset if the error occurs... ... ...")
  print("sorry.")
  
  #now calculate the gene likelihood scores (GLs)  
  startTime<-Sys.time()
  
  sampleGLs<<-calculateAllGLS(unlist(validGenes), treeHypotheses,genesDirectory )
  
  endTime<-Sys.time()
  timeLength<-endTime-startTime
  paste("Done after",timeLength[[1]],attr(timeLength, "units"))
  
  #export deltaGLS values
  
  print("Exporting results...")
  print("Note: deltaLnL values are positive if the first hypothesis is better than the second hypothesis. If it's negative, then it's the reverse. Hypotheses are ordered as given in the file name. e.g. in deltaGLS.Hyp1.vs.Hyp2.csv positive values indicate hypothesis 1 (the first tree in your input file) had a higher lnL than hypothesis 2 (the second tree).")
  
  write.csv(sampleGLs[1]-sampleGLs[2], file=paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[2]),"csv", sep=".")) #if this is positive then the second term has a better lnL
  
  if(length(treeHypotheses)==3){
    write.csv(sampleGLs[1]-sampleGLs[3], file=paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[3]),"csv", sep="."))
    write.csv(sampleGLs[2]-sampleGLs[3], file=paste("deltaGLS",colnames(sampleGLs[2]), "vs",colnames(sampleGLs[3]),"csv", sep="."))
  }
  
  if(length(treeHypotheses)==4){
    write.csv(sampleGLs[1]-sampleGLs[3], file=paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[3]),"csv", sep=".")) 
    
    write.csv(sampleGLs[1]-sampleGLs[4], file=paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[4]),"csv", sep=".")) 
    
    write.csv(sampleGLs[2]-sampleGLs[3], file=paste("deltaGLS",colnames(sampleGLs[2]), "vs",colnames(sampleGLs[3]),"csv", sep="."))
    
    write.csv(sampleGLs[2]-sampleGLs[4], file=paste("deltaGLS",colnames(sampleGLs[2]), "vs",colnames(sampleGLs[4]),"csv", sep="."))
    
    write.csv(sampleGLs[3]-sampleGLs[4], file=paste("deltaGLS",colnames(sampleGLs[3]), "vs",colnames(sampleGLs[4]),"csv", sep="."))
    
  }
  
  plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[2]),"csv", sep=".")) #if this is positive then the second term has a better lnL
  
  if(length(treeHypotheses)>=3){
    plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[3]),"csv", sep="."))
    plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[2]), "vs",colnames(sampleGLs[3]),"csv", sep="."))
  }
  
  if(length(treeHypotheses)==4){
    plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[1]), "vs",colnames(sampleGLs[4]),"csv", sep="."))
    plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[2]), "vs",colnames(sampleGLs[4]),"csv", sep="."))
    plotDeltaGls(paste("deltaGLS",colnames(sampleGLs[3]), "vs",colnames(sampleGLs[4]),"csv", sep="."))
  }
  

}


#dependency function for finding alignments that are valid for the hypotheses in question
#
checkGeneValidity<-function(fileName, directory){
  aGT<-read.dna(paste(directory, "\\", fileName, sep=""), format = "fasta")
  {temp=if(
    length(intersect(cluster1[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
    length(intersect(cluster2[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
    length(intersect(cluster3[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
    length(intersect(cluster4[[2]], attr(aGT, "dimnames")[[1]]))>=1){
    fileName}else{NULL}#if true, gene is appropriate for the test
  }
  return(temp)
}

#dependency function for calculatng a single set of GLS (gene likelihood scores)
singleGLS<-function(geneName,treeHypotheses, directory ){
  
  #load the gene data and convert it to a phyDat object
  geneData<-as.phyDat(read.dna(paste(directory, "\\", 
                                     paste(    strsplit(    geneName, split="\\.")[[1]][1]    , "fasta", sep=".")
                                     , sep=""), format = "fasta"))
  
  #find taxa that are in the alignment but not the tree and drop them from the tree. If there are taxa in the alignment that aren't in the tree then this will not work.
  validTaxa<-intersect(names(geneData), treeHypotheses[[1]]$tip.label)
  taxaToDrop<-setdiff(treeHypotheses[[1]]$tip.label, validTaxa)
  
  
  #calculate the lnL of the gene data given the tree.
  {lnLs<-c()
    for(tree in 1:length(treeHypotheses)){
      lnL<-pml(drop.tip(treeHypotheses[[tree]], taxaToDrop), data=geneData, bf="empirical", model = "GTR", site.rate = "gamma")
      lnLs[tree]<-lnL$logLik
    }
    lnLs
  }
}

#dependency function which loops the above over a set of genes

calculateAllGLS<-function(genes,treeHypotheses, directory){
  
  allGls<-foreach(gene = 1:length(genes),.combine =  'rbind' , .packages=c("ape", "phangorn"), .export = "singleGLS")%do%
    {
      print(paste(gene," - processing ...",genes[[gene]]))
      singleGLS(genes[[gene]],treeHypotheses, directory )
    }
  
  allGls<-as.data.frame(allGls)
  rownames(allGls)<-genes
  colnames(allGls)<-unlist(foreach(i=1:length(treeHypotheses))%do%{paste("Hyp", i, sep="")})
  return(allGls)
}




##dependency function for plotting GLS
##
plotDeltaGls<-function(fileNAme){

{deltaGLsValues<-read.csv(fileNAme)
  title<-strsplit(fileNAme, split="\\.")[[1]][c(2:4)]
  
  outfileName<- paste(paste(title[[1]], title[[2]], title[[3]], sep=""), "svg", sep=".")
  
  colnames(deltaGLsValues)<-c('loci', 'deltalnL')
  
  sortedData<-deltaGLsValues[order(deltaGLsValues$'deltalnL'),]}
orders<-1:length(sortedData$'deltalnL')##add the sorting order as a factor for ease of plotting
plotData<-data.frame(orders, sortedData$'deltalnL')

plot<-{ggplot(sortedData, aes(x=orders, y=deltalnL)) + geom_point(size=.5, color="blue")+scale_color_brewer(palette="Dark2") + theme_minimal()+geom_hline(yintercept=0)}
plot

ggsave(file=outfileName, plot=plot)

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


