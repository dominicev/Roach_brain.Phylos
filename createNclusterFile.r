####Function that takes the cluster list and turns it into a cluster file
# GPT-4 generated code as of April 2023
listToClusterFile<-function(listOfClusters, fileName){
  output_lines<-sapply(names(listOfClusters), function(name) {
                    elements <- paste(listOfClusters[[name]], collapse = ", ")
                    return(paste(name, " = ", elements, ";", sep = ""))
                    }
                  )

  # Write to a text file
  writeLines(output_lines, paste0(fileName, ".txt"))
}
  
#listToClusterFile(listOfClusters,"test")


###make cluster file
createClusterFile<-function(cladeDataFrame, clusterNames, outputFileName = "unnamed.clusterFile"){
  
  #Initialize list
  listOfClusters<-list()
  
  #create defined clusters
  for (clusterName in clusterNames) {
    # Append the unlisted data from cladeDataFrame for the current clusterName
    listOfClusters[[length(listOfClusters) + 1]] <- unlist(cladeDataFrame[clusterName])
    names(listOfClusters[[length(listOfClusters)]]) <- NULL
    object_name_as_string <- ifelse(length(clusterName) > 1, 
                                    paste0("cluster", as.character(length(listOfClusters))), 
                                    clusterName)
    names(listOfClusters)[[length(listOfClusters)]] <- object_name_as_string
  }
  
  #create IGNORE cluster if any
  if (length(unlist(listOfClusters)) == length(tree$tip.label)) {
    # If true then there is nothing being ignored
    print("All taxa accounted for. Not ignoring anything.")
  } else {
    # If False, turn into ignore line
    taxa2Ignore <- setdiff(tree$tip.label, unlist(listOfClusters))
    listOfClusters[[length(listOfClusters) + 1]] <- taxa2Ignore
    names(listOfClusters[[length(listOfClusters)]]) <- NULL
    names(listOfClusters)[length(listOfClusters)] <- "IGNORE"
  }
  
  
  listToClusterFile(listOfClusters,outputFileName)
  
  
}

########EXAMPLE USAGE#############
#clusterNames<-list("Odonata", Polyneoptera, Paraneoptera, Holometabola)

#createClusterFile(cladeDataFrame, clusterNames, "test.1")




