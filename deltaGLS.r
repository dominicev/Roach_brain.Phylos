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



#######Plot DeltaGLS values against locus stats############################
#allStats = a dataframe where the first column (named "LocusName") are the names of the loci (without file suffix) and all other columns are the states for that locus. Each column should be named as that stat "e.g., "Rate", "Saturation" etc.
#directory = the directory where the output of deltaGLS function can be found. Specifically, it needs the files named deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv which has the deltaGLS values for each locus across two competing hypotheses.

correlateDeltaGLSAndStats <- function(allStats, directory) {
  
  setwd(directory)
  # Get all files with the specified pattern
  files <- list.files(pattern = "deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv")
  #print(as.character(dir()))
  #file<-files[1]
 
  
  for (file in files) {
    # Read the file
    deltaGLsValues <- read.csv(file)
    title <- strsplit(file, split="\\.")[[1]][c(2:4)]
    print(paste0(title, collapse = " "))
    # Make the locus names match
    colnames(deltaGLsValues)<-c('loci', 'deltalnL')
    for(i in 1:length(deltaGLsValues[[1]])){
      deltaGLsValues[i,]$'loci'<-strsplit( deltaGLsValues[i,]$'loci'  , split="\\.")[[1]][[1]]
    }
    
    #noOutlierData[(noOutlierData < quantile(noOutlierData, 0.01, na.rm = TRUE))] <- median(noOutlierData, na.rm = TRUE)
    
    
    # Merge with allStats (assuming allStats is already defined in the environment)
    newData<-merge(deltaGLsValues, allStats, by.x = "loci", by.y = "LocusName")
    
    
    
    # Iterate over each column after the third column
    for (i in 3:ncol(newData)) {
      #The next line removes the outlier values from the data so that the correlation isn't driven by a few extreme outliers
      newData[,i][ newData[,i] < quantile(newData[,i], 0.01, na.rm = TRUE)]<- median(newData[,i], na.rm = TRUE)
      newData[,i][ newData[,i] > quantile(newData[,i], 0.99, na.rm = TRUE)]<- median(newData[,i], na.rm = TRUE)
      
      plotData <- data.frame(deltalnL = newData$deltalnL, Value = newData[, i])
      print(paste0(as.character(colnames(plotData))))
      # Create the plot
      plot <- ggplot(plotData, aes(x=Value, y=deltalnL)) +
        geom_point(size=.5, color="blue") +
        geom_smooth(method=lm, color="black") +
        scale_color_brewer(palette="Dark2") +
        theme_minimal() +
        geom_hline(yintercept=0) +
        stat_regline_equation(label.x=.1, label.y=50) +
        stat_cor(aes(label=..rr.label..), label.x=.1, label.y=35)
      
      colnames(plotData) <- c("deltalnL", colnames(newData)[i])
      # Save the plot
      filename <- paste("deltalnL", colnames(newData)[i], title[[1]], title[[2]], title[[3]], "png", sep=".")
      #print(filename)
      ggsave(file=filename, plot=plot)
    }
  }
}



#######Quantitative stats for deltaGLS comparisons############################
#directory = the directory where the output of deltaGLS function can be found. Specifically, it needs the files named deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv which has the deltaGLS values for each locus across two competing hypotheses.

comparingOverallDeltaGLSStats <- function(directory) {
  
  setwd(directory)
  # Get all files with the specified pattern
  files <- list.files(pattern = "deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv")
  #print(as.character(dir()))
  #file<-files[1]
  
  results_list <- list()
 for (file in files) {
    # Read the file
    deltaGLsValues <- read.csv(file)
    title <- strsplit(file, split="\\.")[[1]][c(2:4)]
    #print(title)
    # Make the locus names match
    colnames(deltaGLsValues)<-c('loci', 'deltalnL')
    
    for(i in 1:length(deltaGLsValues[[1]])){
      deltaGLsValues[i,]$'loci'<-strsplit( deltaGLsValues[i,]$'loci'  , split="\\.")[[1]][[1]]
    }
    
    ####Statistical test to see if the mean is different than 0...with outliers
    normalityTestData<-ks.test(deltaGLsValues$deltalnL, "pnorm", mean(deltaGLsValues$deltalnL), sd(deltaGLsValues$deltalnL)) #if the value is <0.05 then the data differ from normality
    
    if(normalityTestData$p.value>0.05){
      testResult0<-t.test(deltaGLsValues$deltalnL, mu = 0)
    }else{
      testResult0<-wilcox.test(deltaGLsValues$deltalnL, mu = 0, paired = FALSE)                  
    }
    #str(testResult0)
    
    
    ####Statistical test to see if the mean is different than 0...without outliers
    noOutlierData<-deltaGLsValues$deltalnL
    noOutlierData[(noOutlierData < quantile(noOutlierData, 0.01, na.rm = TRUE))] <- median(noOutlierData, na.rm = TRUE)
    noOutlierData[(noOutlierData > quantile(noOutlierData, 0.99, na.rm = TRUE))] <- median(noOutlierData, na.rm = TRUE)
    
    normalityTestData<-ks.test(noOutlierData, "pnorm", mean(noOutlierData), sd(noOutlierData)) #if the value is <0.05 then the data differ from normality
    
    
    if(normalityTestData$p.value>0.05){
      testResult1<-t.test(noOutlierData, mu = 0)
    }else{
      testResult1<-wilcox.test(noOutlierData, mu = 0, paired = FALSE)                  
    }
   
    
    ####Statistical test to see if the number of genes supporting one hypothesis versus the other is skewed.
    # Calculate the standard deviation of your dataset
    std_dev <- sd(deltaGLsValues$deltalnL)
    # Number of values in your dataset
    n <- length(deltaGLsValues$deltalnL)
    # Number of simulations
    num_simulations <- 1000
    # Store the count of positive values for each simulation
    positive_counts <- numeric(num_simulations)
    # Perform simulations
    set.seed(123)  # for reproducibility
    for(i in 1:num_simulations) {
      simulated_data <- rnorm(n, mean = 0, sd = std_dev)
      positive_counts[i] <- sum(simulated_data > 0)
    }
    
    # Actual number of positive values in your dataset
    actual_positives <- sum(deltaGLsValues$deltalnL > 0)
    # Compare the actual number of positives to the simulated distribution
    quantile<-sum(actual_positives>positive_counts)/length(positive_counts)
    if(quantile>=0.5){
      pvalueNumberOfGenes<-1-quantile
    }else{pvalueNumberOfGenes<-quantile}
    
    
    ####How many outliers are in support of each hypothesis
    outliersSupportingFirst<-sum(deltaGLsValues$deltalnL > quantile(deltaGLsValues$deltalnL, 0.95, na.rm = TRUE))###outliers in favor of first HYP
    outliersSupportingSecond<-sum(deltaGLsValues$deltalnL < quantile(deltaGLsValues$deltalnL, 0.05, na.rm = TRUE))###outliers in favor of second HYP
    
    file_results <- c(
      paste0(title, collapse=" "), # name of comparison
      if(
        round(mean(deltaGLsValues$deltalnL), 4)>0 &&
        round(median(deltaGLsValues$deltalnL), 4)>0 &&
        round(pvalueNumberOfGenes, 4)<0.050001 &&
        (round(testResult1$p.value, 4)<0.050001||round(testResult0$p.value, 4)<0.050001)){"First hypothesis has more support"}else{
          if(round(mean(deltaGLsValues$deltalnL), 4)<0 &&
             round(median(deltaGLsValues$deltalnL), 4)<0 &&
             round(pvalueNumberOfGenes, 4)<0.050001&&
             (round(testResult1$p.value, 4)<0.050001||round(testResult0$p.value, 4)<0.050001)){"Second hypothesis has more support"}else{
               "Both are supported, or it's complicated"
               }
          
          
          
        },#Summary result
      round(mean(deltaGLsValues$deltalnL), 4), # mean value, rounded
      round(median(deltaGLsValues$deltalnL), 4), # median value, rounded
      round(actual_positives/length(deltaGLsValues$deltalnL), 4), # proportion of genes supporting H1, rounded
      round(pvalueNumberOfGenes, 4), # probability that the number of genes supporting each hypothesis is equal, rounded
      strsplit(testResult0$method, " ")[[1]][1], # test used to see if mean of deltaGLS values is different than 0...with outliers
      round(testResult0$p.value, 4), # probability that mean of deltaGLS values is indistinguishable from 0...with outliers, rounded
      strsplit(testResult1$method, " ")[[1]][1], # test used to see if mean of deltaGLS values is different than 0...without outliers
      round(testResult1$p.value, 4), # probability that mean of deltaGLS values is indistinguishable from 0...without outliers, rounded
      outliersSupportingFirst, # number of outlier genes supporting H1
      outliersSupportingSecond # number of outlier genes supporting H2
    )
    
    results_list[[file]] <- file_results
    
    
 }
  
  results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
  colnames(results_df) <- c("Hyps_Compared","TestSummary", "MeanlnL", "MedianlnL", "ProportionFavoringFirstHyp", "PValueNumGenes", "TestlnLDist_W.Outliers", "PValueLnLDist.W.Outliers", "TestlnLDist_WO.out.Outliers", "PValueLnLDist.WO.Outliers", "OutliersFirstH", "OutliersSecondH")
  
  rownames(results_df)<-NULL
  
  # Return the dataframe
  return(results_df)
}


#######ANOVA to compare deltaGLS values by certain stats############################
#allStats = a dataframe where the first column (named "LocusName") are the names of the loci (without file suffix) and all other columns are the states for that locus. Each column should be named as that stat "e.g., "Rate", "Saturation" etc.
#directory = the directory where the output of deltaGLS function can be found. Specifically, it needs the files named deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv which has the deltaGLS values for each locus across two competing hypotheses.
#modelFormula = the linear model formula to pass to the aov() function. The variables are the colnames in your allSTats obj An example is: deltalnL ~ Saturation + ASRV 


anovaDeltaGLSAndStats <- function(allStats, directory, modelFormula) {
  setwd(directory)
  # Get all files with the specified pattern
  files <- list.files(pattern = "deltaGLS\\.Hyp.*\\.vs\\.Hyp.*\\.csv")
  
  anovaResults <- list()
  
  for (file in files) {
    # Read the file
    deltaGLsValues <- read.csv(file)
    title <- strsplit(file, split="\\.")[[1]][c(2:4)]
    label <- paste0(title, collapse = " ")
    
    # Make the locus names match
    colnames(deltaGLsValues) <- c('loci', 'deltalnL')
    deltaGLsValues$loci <- sapply(deltaGLsValues$loci, function(x) strsplit(x, split="\\.")[[1]][1])
    
    newData <- merge(deltaGLsValues, allStats, by.x = "loci", by.y = "LocusName")
    
    # Process newData (median replacement, dealing with outliers, etc.)
    # Your existing code for processing newData goes here
    
    # Perform ANOVA using the provided model formula
    multiFactor <- aov(modelFormula, data = newData)
    
    # Store the summary of the ANOVA results with the label
    anovaResults[[label]] <- summary(multiFactor)
  }
  
  return(anovaResults)
}
