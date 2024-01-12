# Median Absolute Deviation (MAD) Method for dropping extreme outliers; FROM GPT4
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