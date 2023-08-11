
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

###Unused packages###
#library(rlist)
    #Unused packages determined by...
    #funchir::stale_package_check('B:/OneDrive - University of Illinois - Urbana/Science/R programs/Roach_brain_Phylos/tree functions.R')

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


# Median Absolute Deviation (MAD) Method for dropping extremely outliers; FROM GPT4
remove_mad_outliers <- function(x, multiplier = 3.5) {
  mad <- mad(x, constant = 1)
  med <- median(x)
  
  # Identify outliers
  outliers <- abs(x - med) / mad > multiplier
  
  # Return data without outliers
  return(x[!outliers])
}

###removing taxa that have extremely large pairwise distances###
###This will compare the pairwise tip distances in the whole tree to those generated under an estimated gamma distribution and remove those that appear to have numerous cases (relative to the total number of taxa) where they are outside the alpha (default 95%) quantile of that distribution.
###
#troubelshooting
#setwd("G:\\My Drive\\Projects\\PHY_Blab_2\\bioinformatics\\PhyBlab2 full analysis\\loci\\3 - gene trees\\Full gene tree analysis 28-12-21")


#treeFileNames<-list.files(pattern=".treefile")
#alpha = 0.95
#relAbund = "Default"
#randomSeed = 2345


trimPairTipOutliers<-function(tree, alpha = 0.95, relAbund = "Default", randomSeed = 2345, fitdistsMethod = "mle" ) {
      ##Import tree
      phylo<-tree
      #plot.phylo(phylo)
      
      #Define the number of times the taxon name has to appear in list of long branch taxa to be considered a true long branch taxa
      if(relAbund=="Default"){relAbund<-(length(phylo[[4]])/4)}
      
      #calculate root to tip branch length distances
      #dists<-distTips(phylo, tips = "all", method = "patristic") #old method/slow
      dists<-get_all_pairwise_distances(phylo, only_clades = phylo$tip.label)
      
      rownames(dists)<-phylo$tip.label
      colnames(dists)<-phylo$tip.label
      dists<-as.dist(dists)

      #reshape dists into a list
      df <- melt(as.matrix(dists), varnames = c("row", "col"))
      allDists<-df$value
      #hist(allDists)
      
      #fit gamma distribution to distances
      set.seed(randomSeed)
      fit.gamma = fitdist(allDists[allDists>0], distr = "gamma", method = fitdistsMethod)
      #fit.gamma[1]
      
      #Calculate outliers based on random sampling of gamma distribution
      randomGammaValues<-rgamma(9999, fit.gamma$estimate[[1]], rate = fit.gamma$estimate[[2]] )
      
      #shortBranches<-quantile(randomGammaValues, c(.001))
      longBranches<-quantile(randomGammaValues, c(alpha))
      
      
      ####SECTION BELOW WORKS BUT MIGHT BE TOO SLOW
      #rowNames<-df$row; colNames<-df$col
      #longBranchTaxa<-unlist(c(rowNames[df$value>longBranches],colNames[df$value>longBranches]))
      #frequencies<-sort(table(longBranchTaxa), decreasing = TRUE)
      #strictLBTaxa<-dimnames(frequencies[frequencies>relAbund])
      
      #delete outlier tips and export
      
      #ifelse(length(strictLBTaxa)>0,phylo<-drop.tip(  phylo,strictLBTaxa$longBranchTaxa ), 0)#redefines the phylogeny
      
      #####SECTION BELOW IS THE REWRITTEN PART FROM GPT4
      longBranchIndices <- which(df$value > longBranches)
      longBranchTaxa <- unlist(c(df$row[longBranchIndices], df$col[longBranchIndices]))
      
      frequencies <- sort(table(longBranchTaxa), decreasing = TRUE)
      strictLBTaxa <- dimnames(frequencies[frequencies > relAbund])
      
      # Delete outlier tips and redefine the phylogeny
      ifelse(length(longBranchTaxa) > 0, phylo<-drop.tip(phylo, strictLBTaxa$longBranchTaxa), 0)
      
      
      #write.tree(phylo, paste(strsplit(files, ".", fixed=TRUE)[[1]][[1]],".tiptrimmed.tre", sep=""))
      print(paste(length(strictLBTaxa$longBranchTaxa), "tips dropped from tree"))
      
      return(phylo)
      
      
    
}

###removing taxa that have extremely large root-to-tip distances###
###This will compare the pairwise tip distances in the whole tree to those generated under an estimated gamma distribution and remove those that appear to have numerous cases (relative to the total number of taxa) where they are outside the alpha (default 95%) quantile of that distribution.
###
### rooting option will take a list of rooting taxa in the order of priority (i.e., the input for rootIncompleteTree)

#outgroupTaxonList<-rooting<-outgroupList
####this is the same as the above function but it operates on a single tree at a time
trimRootTipOutliers<-function(tree, alpha = 0.99, rooting = "Default", randomSeed = 2345, fitdistsMethod = "mle") {

    ##Import tree
    phylo<-tree
    outgroupTaxonList<-rooting
  
    #root tree
    ifelse(rooting == "Default", phylo<-midpoint.root(phylo), phylo<-rootIncompleteTree(phylo, outgroupTaxonList))
    
    #calculate root to tip branch length distances
    #dists=distRoot(phylo, tips = "all", method = "patristic") #old method/slow
    dists<-get_all_distances_to_root(phylo)[1: length(phylo$tip.label)]
    names(dists)<-phylo$tip.label
    
    #fit gamma distribution to distances
    set.seed(randomSeed)
    
    # Remove outliers in dists so the function does a better job of dropping most on the first run
    #dists_no_outliers <- remove_mad_outliers(dists)
    
    fit.gamma = fitdist(dists, distr = "gamma", method = fitdistsMethod)
    
    #Calculate outliers based on random sampling of gamma distribution
    randomGammaValues<-rgamma(9999, fit.gamma$estimate[[1]], rate = fit.gamma$estimate[[2]] )
    
    longBranches<-quantile(randomGammaValues, c(alpha))
    
    #delete outlier tips and export
    #{ 
    #  longBranchTaxa<-c()
    #  j=0;k=0
    #  for(value in dists){
    #    j=j+1;
    #    name=names(dists[j])
    #    if(value>longBranches[[1]]){k=k+1;
    #    longBranchTaxa[k]<-name}
    #  }
    
    #}###gets the taxon names of outliers
    
    ###ABOVE is the old method, which works but is slow. Below is the new method which was written by GPT4
    longBranchIndices <- which(dists > longBranches[[1]])
    longBranchTaxa <- names(dists[longBranchIndices])
    
    #longBranchTaxa
   ifelse(length(longBranchTaxa) > 0, phylo<-drop.tip(phylo, longBranchTaxa), 0) #redefines the phylogeny
    
    #write.tree(phylo, paste(strsplit(files, ".", fixed=TRUE)[[1]][[1]],".tiptrimmed.tre", sep=""))
    print(paste(length(longBranchTaxa), "tips dropped from tree"))
    return(phylo)
  
  
}


###recursive tip trimming function


recursiveTrimOutlierTaxa<-function(tree, alpha = 0.999, rooting = "Default", randomSeed = 2345, fitdistsMethod = "mle"){
  phylo<-tree
  
  repeat{
    previousTaxaLength<-length(phylo$tip.label)
    phylo<-trimRootTipOutliers(
      trimPairTipOutliers(phylo, alpha = alpha, relAbund = "Default", randomSeed = randomSeed, fitdistsMethod = fitdistsMethod )
      , alpha = alpha, rooting = rooting, randomSeed= randomSeed, fitdistsMethod = fitdistsMethod)
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
expr <- "x^2"
var_list <- list("x", 1, 5)
result <- tableW(expr, var_list)
print(result)




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

                                    