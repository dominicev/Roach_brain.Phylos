library(ggplot2)
library(svglite)

######deltaGLS (gene likelihood scores) test######

#deltaGLS was used in Shen et al. (2017) to show that contentious relationships can be driven by outlier genes.
#These functions will allow calculation of these values from a set of genes given 2-4 alternative hypotheses

#input...
#clusterFileName is a string specifying the name of the cluster file. Unless it's in your working directory, you need to specify the whole directory+file name
#treeHypotheses must be a multi.phylo object containing your different topological hypotheses. It will only work if this has between 2 and 4 trees in it. Any more and you need to do multiple runs.
#genesDirectory the directory where your alignment files are. These must be in fasta format and have the file extension ".fasta"
#GTs optional

deltaGLS<-function(clusterFileName, treeHypotheses, genesDirectory, GTs = c(), QtestModels = FALSE){
  #read cluster file, which must have 4 clusters and one list of taxa to ignore. See example file for format.
  allClusters<-read.delim(clusterFileName, header=FALSE, sep = "=")
  #cluster1<<-c(strsplit(allClusters[1, 1], " +")[[1]], lapply(strsplit(allClusters[1, 2], ", "), trimws ) )
  #cluster2<<-c(strsplit(allClusters[2, 1], " +")[[1]], lapply(strsplit(allClusters[2, 2], ", "), trimws ) )
  #cluster3<<-c(strsplit(allClusters[3, 1], " +")[[1]], lapply(strsplit(allClusters[3, 2], ", "), trimws ) )
  #cluster4<<-c(strsplit(allClusters[4, 1], " +")[[1]], lapply(strsplit(allClusters[4, 2], ", "), trimws ) )
  
  trimTrailingSemicolon <- function(x) {
    sub(";$", "", x)  # Removes the trailing semicolon if it exists
  }
  processCluster <- function(clusterRow) {
    c(strsplit(allClusters[clusterRow, 1], " +")[[1]], 
      lapply(strsplit(trimTrailingSemicolon(allClusters[clusterRow, 2]), ", "), trimws))
  }
  
  cluster1 <<- processCluster(1)
  cluster2 <<- processCluster(2)
  cluster3 <<- processCluster(3)
  cluster4 <<- processCluster(4)
  
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
      print(paste(as.character(i), GTs[[i]], sep = ". "))
      checkGeneValidity(GTs[[i]],genesDirectory )
      }
    print("Done checking...writing validGene.csv file")
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
  
  sampleGLs<<-calculateAllGLS(unlist(validGenes), treeHypotheses,genesDirectory, QtestModels )
  
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
#checkGeneValidity<-function(fileName, directory){
#  aGT<-read.dna(paste(directory, "\\", fileName, sep=""), format = "fasta")
#  {temp=if(
#    length(intersect(cluster1[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
#    length(intersect(cluster2[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
#    length(intersect(cluster3[[2]], attr(aGT, "dimnames")[[1]]))>=1 && 
#    length(intersect(cluster4[[2]], attr(aGT, "dimnames")[[1]]))>=1){
#    fileName}else{NULL}#if true, gene is appropriate for the test
#  }
#  return(temp)
#}


###current newest version below

checkGeneValidity <- function(fileName, directory){
  # Function to check if a sequence is DNA
  isDNA <- function(sequence) {
    all(toupper(sequence) %in% c("A", "C", "G", "T", "N", "-", "?"))
  }
  
  # Read the file content
  filePath <- paste(directory, "/", fileName, sep="")
  fileContent <- readLines(filePath, n = 2)[2]  # Read only the first line to check the sequence type
  #print(paste("File:", fileName, "First line:", fileContent))
  
  # Determine if the file is DNA or AA based on the first sequence
  if (isDNA(fileContent)) {
    #print("Detected as DNA")
    aGT <- read.phyDat(filePath, format = "fasta", type = "AA")
  } else {
    #print("Detected as Protein")
    aGT <- read.phyDat(filePath, format = "fasta", type = "AA")
  }
  
  # Check if the gene is valid
  isValid <- length(intersect(cluster1[[2]], names(aGT))) >= 1 &&
    length(intersect(cluster2[[2]], names(aGT))) >= 1 &&
    length(intersect(cluster3[[2]], names(aGT))) >= 1 &&
    length(intersect(cluster4[[2]], names(aGT))) >= 1
  
  #print(paste("Is valid:", isValid))
  result <- if(isValid) fileName else NULL
  #print(paste("Result for", fileName, ":", result))
  return(result)
}



#OLD VERSION#dependency function for calculatng a single set of GLS (gene likelihood scores)
#singleGLS<-function(geneName,treeHypotheses, directory ){
#  
#  #load the gene data and convert it to a phyDat object
#  geneData<-as.phyDat(read.dna(paste(directory, "\\", 
#                                     paste(    strsplit(    geneName, split="\\.")[[1]][1]    , "fasta", sep=".")
#                                     , sep=""), format = "fasta"))
#  
#  #find taxa that are in the alignment but not the tree and drop them from the tree. If there are taxa in the alignment that aren't in the tree then this will not work.
#  validTaxa<-intersect(names(geneData), treeHypotheses[[1]]$tip.label)
#  taxaToDrop<-setdiff(treeHypotheses[[1]]$tip.label, validTaxa)
#  
# 
#  #calculate the lnL of the gene data given the tree.
#  {lnLs<-c()
#    for(tree in 1:length(treeHypotheses)){
#      lnL<-pml(drop.tip(treeHypotheses[[tree]], taxaToDrop), data=geneData, bf="empirical", model = "GTR", site.rate = "gamma")
#      lnLs[tree]<-lnL$logLik
#    }
#    lnLs
#  }
#}



singleGLS <- function(geneName, treeHypotheses, directory, QtestModels = FALSE) {
  
  isDNA <- function(sequence) {
    all(toupper(sequence) %in% c("A", "C", "G", "T", "N", "-", "?"))
  }
  
  # Read the file content
  geneFilePath <- paste(directory, "/", geneName, sep="")
  fileContent <- readLines(geneFilePath, n = 2)[2]  # Read only the first line to check the sequence type
  #print(paste("File:", fileName, "First line:", fileContent))

  
  # Convert gene data to phyDat object
  if (isDNA(fileContent)) {
    geneData <- read.phyDat(geneFilePath, format = "fasta", type = "DNA")
  } else {
    geneData <- read.phyDat(geneFilePath, format = "fasta", type = "AA")
  }
  
  # Find taxa that are in the alignment but not in the tree and drop them from the tree
  validTaxa <- intersect(names(geneData), treeHypotheses[[1]]$tip.label)
  taxaToDrop <- setdiff(treeHypotheses[[1]]$tip.label, validTaxa)
  
  consensusTree<-consensus(treeHypotheses, check.labels = TRUE)
  
  if (QtestModels) {
    tryCatch({
      if (isDNA(fileContent)) {
        # Find the best model for nucleotide data
        modelTestResult <- modelTest(geneData, tree=drop.tip(consensusTree, taxaToDrop), model = c("GTR", "SYM"), multicore = TRUE, mc.cores = 3)
      } else {
        # Find the best model for amino acid data
        attr(geneData, "type") <- "AA"
        modelTestResult <- modelTest(geneData, tree=drop.tip(consensusTree, taxaToDrop), model = c("JTT", "WAG", "LG", "Dayhoff"))
      }
      bestModel <- modelTestResult$bestModel
    }, error = function(e) {
      # Fallback to predetermined model in case of an error
      if (isDNA(fileContent)) {
        bestModel <- "GTR"
      } else {
        bestModel <- "Dayhoff"
      }
      print("Using generic model due to error in model testing")
    })
  } else {
    # Default to predetermined models without testing
    if (isDNA(fileContent)) {
      bestModel <- "GTR"
    } else {
      bestModel <- "Dayhoff"
    }
    print("Using generic model as model testing is skipped")
  }
  
  
  
  # Loop through each tree hypothesis
  lnLs <- c()
  for (tree in 1:length(treeHypotheses)) {
    # Calculate likelihood with the best model or predetermined model
    lnL <- pml(drop.tip(treeHypotheses[[tree]], taxaToDrop), data = geneData, model = bestModel)
    lnLs[tree] <- lnL$logLik
  }
  
  
  
  
  lnLs
}

# Usage example
# singleGLS("gene_name.fasta", list_of_tree_hypotheses, "path/to/directory")




#dependency function which loops the above over a set of genes

calculateAllGLS<-function(genes,treeHypotheses, directory, QtestModels = FALSE){
  
  allGls<-foreach(gene = 1:length(genes),.combine =  'rbind' , .packages=c("ape", "phangorn"), .export = "singleGLS")%do%
    {
      print(paste(gene," - processing ...",genes[[gene]]))
      singleGLS(genes[[gene]],treeHypotheses, directory, QtestModels  )
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
