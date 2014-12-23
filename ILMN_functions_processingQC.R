#=============================================================================#
# FUNCTIONS                                                                   #
# ArrayAnalysis - ILMN_functions_processingOC                                 #
# a tool for quality control and pre-processing of Illumina array data        #                                                                        
#=============================================================================#

####################
## checkUserInput ##
####################

## Check if given species, array type and array annotation combi is ok       
## --------------------------------------------------------------------

checkUserInput <-function(arrayTypes, arrayAnno, species) {
  check <- species %in% names(arrayTypes) && arrayType %in% arrayTypes[[species]] && annoType %in% arrayAnno[[arrayType]]
  if (!check) { 
    stop ('\n' , "The combination of species, array type and array annotation file is not correct:", '\n' ,
          "- Species: ", species, '\n' ,
          "- Array type: ",arrayType, '\n' ,
          "- Annotation file: ", annoType, sep=" ")
  } else {
    cat("..::..::..\n", 
        "CHECK: ALL THE SPECIFIED USER INPUT ARRAY VARAIBLES OK.\n",sep="")
  
  }
  #return(check)
}


####################
## import.rawData ##
####################

## Loads illumina rawdata beadarray data
## -------------------------------------

import.rawData <-function (expData, detectionTh ,convertNuID, checkDupId, 
                                lib.mapping, dec, parseColumnName, rawDataQC) {
  rawData = lumiR(expData,
                  detectionTh = detectionTh,
                  convertNuID = convertNuID,  
                  checkDupId = checkDupId, 
                  lib.mapping = lib.mapping,
                  dec = dec,
                  parseColumnName = parseColumnName,
                  QC = rawDataQC
                  )
  return(rawData)
}


###################
## lumi.normData ##
###################

## Pre-processing and normalize data with lumiExpresso function
## returns a normalized lumi batch object
## ------------------------------------------------------------

lumi.normData <-function(x.lumi, 
                         bg.correct, bgCorrect.m,
                         variance.stabilize, variance.m,
                         normalize, normalization.m, 
                         normDataQC) {
  
  normData = lumiExpresso(x.lumi, 
                          bg.correct = bg.correct, 
                          bgcorrect.param = list(method=bgcorrect.m), 
                          variance.stabilize = variance.stabilize, 
                          varianceStabilize.param = list(method=variance.m), 
                          normalize = normalize, 
                          normalize.param = list(method=normalization.m), 
                          QC.evaluation = normDataQC
                          )
  return(normData)  
}


##################
##neqc.normData ##
##################

## pre-processing and normalize data with neqc function
## returns a normalized lumi batch object
## ----------------------------------------------------

neqc.normData <- function (x.lumi, controlData) {
  print("Normalizing data using the 'limma::neqc function' using log2 and quantile normalization")
  #Select only rows with negative probes
  controlData.neg <- controlData[grep("negative", tolower(controlData$controlType)),]
  #remove first 3 columns of controlData.neg
  controlData.neg <- controlData.neg[,3:length(colnames(controlData.neg))]
  #Get matrix of gene expression values:
  expData <- exprs(x.lumi)
  #Combine data, status vector (stating for each row if it corresponds to gene or control)
  totalData <- rbind(expData, controlData.neg)
  status <- c(rep("regular", nrow(expData)), rep("negative", nrow(controlData.neg)))
  #Normalize the data with limma neqc function
  normData.neqc <- limma::neqc(totalData, status)
  #Create lumiBatch object with normalized data
  normData <- x.lumi
  #Overwrite lumiBatch object with neqc normalized gene expression values from normData
  exprs(normData) <- normData.neqc[1:nrow(normData.neqc), ]
  return(normData)
}


#################  
##createSummary##
#################

## Create summary file from lumi.batch QC object 
## ---------------------------------------------

createSummary <- function(x.lumi, fn) {
  #create table which contains QC summary of the raw data
  summary<-x.lumi@QC$sampleSummary
  summary<-t(summary)
  write.table(summary, file=fn, sep ='\t', col.names=T, row.names=T)
  return(summary)
}


########################
## Filtering function ##
########################

## To speed up the processing and reduce false positives, remove the unexpressed probes.
## Beads detection pvalue (<0.01) and at >0 (detectprobes) samples should have this pvalue.
## ---------------------------------------------------------------------------------------- 

filterFun <-function (x.lumiRaw, x.lumiNorm, filter.Th, filter.dp) {
  
  #create datamatrix of pre-processed data.
  dataMatrix = exprs(x.lumiNorm)
  # Estimate which probes are detectable with an detection p.value = filter.Th
  presentCount <- detectionCall(x.lumiRaw, Th=filter.Th) 
  # create subset DataMatrix with only probes which have detectable 
  # probes (filter.dp) in more than 1 sample 
  filtered.normData <- dataMatrix[presentCount > filter.dp,]
  rmProbes = length(featureNames(normData)) - nrow(filtered.normData)
  
  cat ("..::..::..\n",
       "CREATE FILTERED DATA TABLE TO SPEED UP THE PROCESSING AND REDUCE FALSE POSITIVES; REMOVING THE UNEXPRESSED PROBES/GENES:\n",
       "SETTINGS: Rows which do not have p-value < ", filter.Th, " in at least ", filter.dp+1, " column will be excluded.\n",
       " - Normalized data table contains: ", length(featureNames(normData)), " probe rows;\n",
       " - Filtered table contains: ", nrow(filtered.normData), " probe rows;\n", 
       " - Removed ", rmProbes, " unexpressed probe rows from original normalized data file.\n",sep="")
  return(filtered.normData)
}


###################
## createAnnoFun ##
###################

## Create annotation file from a data.frame or 
## matrix. df must contain a first row with nuIDs
## ----------------------------------------------

createAnnoFun <- function(df, ns, lib.mapping, lib.All.mapping) {
  #list of nuIDs 
  nuIDs = rownames(df)
  annoFile = as.data.frame(nuIDs)
  if (require(lib.All.mapping, character.only = TRUE) & require(annotate)) {
    # get anno
    annoFile$ILMN_GENE      <- nuID2IlluminaID(nuIDs, lib.mapping, idType=c("Symbol"))
    annoFile$geneSymbol     <- getSYMBOL(nuIDs, lib.All.mapping)
    annoFile$PROBE_ID       <- nuID2IlluminaID(nuIDs, lib.mapping, idType=c("Probe"))
    annoFile$ENTREZ_GENE_ID <- unlist(lookUp(nuIDs, lib.All.mapping, "ENTREZID"))
    annoFile$GENE_NAME      <- sapply(lookUp(nuIDs, lib.All.mapping, 'GENENAME'), function(x) x[1])
    annoFile$ACCESSION      <- unlist(lookUp(nuIDs, lib.All.mapping, "ACCNUM"))
    annoFile$PROBE_SEQUENCE <- unlist(id2seq(nuIDs)) 
  }
  return(annoFile)
}


####################
## colorsByFactor ##
####################

# create colors for the plots and the legends
# -------------------------------------------

colorsByFactor <- function(experimentFactor) {
  
  #check whether a factor has been provided
  if(class(experimentFactor)!="factor") stop("Parameter 'experimentFactor' must be of class 'factor'")

  if(length(levels(experimentFactor))==1) {
    #if there is only one group (or no groups are provided) take equally spread colors over the rainbow palette
    plotColors <- rainbow(length(experimentFactor),s=.8,v=.7)
  #set group legend color to white, as there is not a specific group color
  legendColors <- "white"
  } else {
    #compute the number of colors needed for each class
    tab.tmp <- table(experimentFactor)

    #set the two extreme colors for each class
    colors.light <- rainbow(length(levels(experimentFactor)),s=1-sapply(tab.tmp,min,5)*.1)
    colors.dark <- rainbow(length(levels(experimentFactor)),v=1-sapply(tab.tmp,min,5)*.14)

    #create the colors to plot, and colors for the legend (average one per experimental group)
    plotColors <- NULL
    legendColors <- NULL
    for(l in 1:length(levels(experimentFactor))) {
      colorFun <- colorRampPalette(c(colors.light[l],colors.dark[l]))
      tmpColors <- colorFun(tab.tmp[l])
      plotColors[experimentFactor==levels(experimentFactor)[l]] <- tmpColors
      legendColors[l] <- tmpColors[ceiling(length(tmpColors)/2)]
    }
  }
  return(list(plotColors=plotColors,legendColors=legendColors))
}
