#=============================================================================#
# ArrayAnalysis - ILMNAnalysisQCandPreProcessing                              #
# a tool for quality control and pre-processing of Illumina array data        #                                       
#=============================================================================#

version_nb <- "1.0.0"
cat("..::..::..\n",
    "Script run using R version ",R.Version()$major,".",R.Version()$minor,
    " and Illumina pre-processing pipeline version ",version_nb,"\n",sep="")

#set memory to maximum on Windows 32bit machines
if(length(grep("w32",R.Version()$os,fixed=TRUE))>0) memory.size(4095)

###############################################################################
## Load ILMN_AnalysisQC functions and R libraries                            ##  				   
###############################################################################

reload <- function() {
  source(paste(SCRIPT.DIR,"install_libraries.R",sep=""))
  source(paste(SCRIPT.DIR,"ILMN_functions_processingQC.R",sep=""))
  source(paste(SCRIPT.DIR,"ILMN_functions_imagesQC.R",sep=""))
  cat("..::..::..\n", 
      "FUNCTIONS HAVE BEEN LOADED.\n", sep="")
  }

## reload functions from source files 
## ----------------------------------
reload();

## Create list of mandatory packages 
## ---------------------------------------------------
pkgs <- c( "limma", "ALL","bioDist", "gplots",
           "annotate", "arrayQualityMetrics",
           switch(species,
                  Human = pkgs <- c("org.Hs.eg.db"),
                  Mouse = pkgs <- c("org.Mm.eg.db"),
                  Rat =   pkgs <- c("org.Rn.eg.db") 
                  ), "lumi",
           paste("lumi", species, "IDMapping", sep=""),
           paste("lumi", species, "All.db", sep="") )

## load and install R libraries (if not found)
## -------------------------------------------------------

cat("..::..::..\n", 
    "LOADING OR INSTALLING R PACKAGES IF MISSING PACKAGES FOUND.\n", sep="")
loadPackages(pkgs);

if(standALONE){
  dateTIME = paste(format(Sys.Date(), "%Y%m%d"), sub(":", ".", sub(":", "." ,format(Sys.time(), "%X"))), sep="_" )
  ns <- removeExt(expFile)
  ns = paste(ns, dateTIME, sep="_")
  }

###############################################################################
## Create arrayTypes and arrayAnnotation lists                               ##
###############################################################################

arrayTypes = list(
Human = c( "HumanHT-12", "HumanRef-8", "HumanWG-6"),
Mouse = c( "MouseRef-8", "MouseWG-6"),
Rat   = c( "RatRef-8")
)

arrayAnno = list(
`HumanHT-12` = c("HumanHT-12_V4_0_R2_15002873_B_WGDASL", 
                 "HumanHT-12_V4_0_R2_15002873_B", 
                 "HumanHT-12_V4_0_R1_15002873_B", 
                 "HumanHT-12_V3_0_R2_11283641_A",
                 "HumanHT-12_V3_0_R3_11283641_A"),
`HumanRef-8` = c("HumanRef-8_V3_0_R3_11282963_A", 
                 "HumanRef-8_V3_0_R2_11282963_A", 
                 "HUMANREF-8_V3_0_R1_11282963_A_WGDASL", 
                 "HumanRef-8_V2_0_R4_11223162_A"),
`HumanWG-6`  = c("HumanWG-6_V2_0_R4_11223189_A", 
                 "HumanWG-6_V3_0_R2_11282955_A",
                 "HumanWG-6_V3_0_R3_11282955_A"),
`MouseRef-8` = c("MouseRef-8_V1_1_R4_11234312_A", 
                 "MouseRef-8_V2_0_R2_11278551_A",
                 "MouseRef-8_V2_0_R3_11278551_A"),
`MouseWG-6` = c( "MouseWG-6_V1_1_R4_11234304_A", 
                 "MouseWG-6_V2_0_R2_11278593_A",
                 "MouseWG-6_V2_0_R3_11278593_A"),
`RatRef-12` = c( "RatRef-12_V1_0_R5_11222119_A")
)

checkUserInput(arrayTypes, arrayAnno, species);

###############################################################################
## Load RAW data, perform pre-processing and normalization                   ##
###############################################################################
## Read in Illumina expression data. 
## We assume all the files were saved in a comma or tab separated text file.

# create path to datafile
expData <-  paste(DATA.DIR, expFile, sep="")

# load raw data as lumi batch file
cat("..::..::..\n", 
    "LOADING RAW DATA FILE: \n", sep="")
rawData <- import.rawData(expData, detectionTh ,convertNuID, 
                          checkDupId, lib.mapping, dec, 
                          parseColumnName, rawDataQC);
cat("..::..::..\n", 
    "RAW DATA FILE LOADED.\n", sep="")

#create generic sampleNames with function make.names
sampleNames(rawData)<- make.names(sampleNames(rawData))


###############################################################################
## Check description file                                                    ##
###############################################################################
cat("..::..::..\n", 
    "CHECKING DISCRIPTION FILE.\n", sep="")

#create new column with format sampleNames as read-in with make.names
description$arraySampleNames = make.names(description[,1])

if( sum( length(sampleNames(rawData)) - length(description[,1]))  > 0)
  stop("error: Number of array names in raw data file and number of array 
       names in description file is not of the same size")

#Match sampleNames from datafile with first column from description file
file_order <- match(description[,4],sampleNames(rawData))

#Check on NA values in file_order; if na in file_order stop
if(sum(is.na(file_order)) > 0) 
  stop("error: Assigned array names in raw data file and file names 
       in description file do not match")

#Check if every sampleName is unique in description file
if(length(description[,2]) != length(unique(description[,2])) ) 
  stop("error: Assigned sampleNames are not unique")

cat("..::..::..\n", 
    "DISCRIPTION FILE OK!\n", sep="")

#Change order of rawData in order of file_order
rawData <- rawData[,file_order]

###############################################################################
## reorder rawData lumibatch file on Group and sampleNames                   ##
###############################################################################
if(perGroup) {
  cat("..::..::..\n", 
      "RE-ORDERING 'RAW' DATA FILE PER GROUP\n", sep="")
  description2 =  description[order(description[,3], description[,2]), ]
  #Match sampleNames from datafile with first column from description file
  file_order2 <- match(description2[,4],sampleNames(rawData))
  if(sum(is.na(file_order2)) > 0) 
    stop("error: file names in data directory and file names in description file do not match")
  rawData <- rawData[,file_order2]
  # change sampleNames into reordered description file
  sampleNames(rawData)<- as.character(description2[,2]) 
  } else {
    # change sampleNames into loaded description file
    sampleNames(rawData)<- as.character(description[,2]) 
    }

# if bgSub = TRUE
if(bgSub) { 
  
  # normalize lumi batch raw data object
  cat("..::..::..\n", 
      paste("NORMALIZING 'RAW' DATA FILE:", expFile), "\n", sep="")
  normData <- lumi.normData(rawData, 
                            bg.correct=FALSE, bgcorrect.m,
                            variance.stabilize, variance.m,
                            normalize, normalization.m, 
                            normDataQC);
  # if bgSub = FALSE                          
  } else {
    cat("..::..::..\n", 
        paste("LOADING 'CONTROL' DATA FILE ", bgFile), "\n", sep="")
    controlData <- paste(DATA.DIR, bgFile, sep="")
    # checks headers of controlData file with rawData object
    # add control data to the rawData lumi batch object file
    cat("..::..::..\n", 
        "COMBINING 'CONTROL' DATA WITH 'RAW' DATA FILE.\n", sep="")
    rawData.ctrl <- addControlData2lumi(controlData, rawData)
    # get control data in a data.frame
    controlData <- as.data.frame(getControlData(rawData.ctrl), row.names = NULL )
    # normalize lumi batch raw data object using 'lumi' or 'neqc' function
    cat("..::..::..\n", 
        paste("NORMALIZING 'RAW' DATA FILE:", bgFile), ".\n", sep="")
    switch (normType,
            lumi = normData <- lumi.normData(rawData.ctrl, 
                                             bg.correct=TRUE, bgCorrect.m,
                                             variance.stabilize, variance.m,
                                             normalize, normalization.m,
                                             normDataQC),
            neqc = normData <- neqc.normData(rawData.ctrl, controlData)
            )
  }

# create rawData exprs eset table
#---------------------------------
eset.rawData <- exprs(rawData)

# create normData exprs eset table
#---------------------------------
eset.normData <- exprs(normData)

###############################################################################
## Creating summary tables                                                   ##
###############################################################################

##raw data                             
if (rawSummary){
  cat("..::..::..\n", 
      "CREATE 'RAW' SUMMARY FILE. \n", sep="")
  rawSum.table = createSummary(rawData, fn=paste(ns,"summary.rawData.txt",sep='_') );
}  
                        
##normalized data
if(normSummary && QC.evaluation){
  cat("..::..::..\n", 
      "CREATE 'NORMALIZED' SUMMARY FILE. \n", sep="")
  normSum.table = createSummary(normData, fn=paste(ns,"summary.normData.txt",sep='_') );
  }                        

###############################################################################
## save lumiBatch file as a R object                                         ##
###############################################################################

# raw data
if(save.rawData) {
  cat("..::..::..\n", 
      "CREATE R OBJECT FILE of 'RAW' DATA.\n", sep="")
  save(rawData, file = paste(ns, 'rawData.Rdata', sep='_') )
  }

# normalized data
if(save.normData) {
  cat("..::..::..\n", 
      "CREATE R OBJECT FILE of 'NORMALIZED' DATA.\n", sep="")
  save(normData, file = paste(ns, 'normData.Rdata', sep='_') )
  }
  
###############################################################################
# Create array groups, array names and plot variables                         #
###############################################################################
cat("..::..::..\n", 
    "CREATING PLOT COLORSET FOR EACH ARRAY GROUP.\n", sep="")
# Create colorset for the array groups
#-------------------------------------
if(perGroup){
  # use reordered  description file ordered per group
  experimentFactor <- as.factor(description2[,3])
  colList          <- colorsByFactor(experimentFactor)
  plotColors       <- colList$plotColors
  legendColors     <- colList$legendColors
  rm(colList)
  } else {
    # use originaly loaded description file
    experimentFactor <- as.factor(description[,3])
    colList          <- colorsByFactor(experimentFactor)
    plotColors       <- colList$plotColors
    legendColors     <- colList$legendColors
    rm(colList)
    }

# Create symbolset for the array groups
#--------------------------------------
plotSymbols <- 18-as.numeric(experimentFactor)
legendSymbols <- sort(unique(plotSymbols), decreasing=TRUE)

###############################################################################
# Define display parameters for the images  		                              #
###############################################################################

WIDTH <- 1000
HEIGHT <- 1414
POINTSIZE <- 24
if(!exists("MAXARRAY")) MAXARRAY <- 41

###############################################################################
## Raw data Quality Control graphs                                           ##
###############################################################################
cat("..::..::..\n", 
    "COMPUTING PLOTS FOR RAW DATA:\n", sep="")

if(raw.boxplot) {
  print ("Plot boxplot for raw intensities") 
  box.plot(rawData, paste(ns, "RAW", sep="_"), col=plotColors, maxArray=50)
  }

if(raw.density) {
  print ("Plot density histogram for raw intensities")
  density.plot(rawData, paste(ns, "RAW", sep="_"), col=plotColors, maxArray=16)
  }

if(raw.cv) {
  print ("Plot density for coefficient of varience for raw intensities")
  cv.plot(rawData, paste(ns, "RAW", sep="_"), col=plotColors, maxArray=16)
  }

if(raw.sampleRelation) {
  print ("Hierarchical clustering of raw data") 
  clusterFun(rawData, normalized=FALSE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }

if(raw.pca) {  
  print("PCA graph for raw data")
  groupsInLegend =  !( length(unique(levels(experimentFactor))) ) >=10
  
  pcaFun(rawData, normalized=FALSE ,experimentFactor=experimentFactor, 
         plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
         legendSymbols=legendSymbols, groupsInLegend=groupsInLegend ,namesInPlot=((max(nchar(sampleNames(rawData)))<=10)&&
           (length(sampleNames(rawData))<=(MAXARRAY/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
         POINTSIZE=POINTSIZE)
  }

if(raw.correl){  
  print ("Correlation plot of raw data")
  correlFun(rawData, normalized=FALSE, experimentFactor=experimentFactor, 
            clusterOption1=clusterOption1, clusterOption2=clusterOption2,
            legendColors=legendColors,
            WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY)  
  }

###############################################################################
## Evaluation of the pre-processing (After normalization graphs)             ##
###############################################################################
cat("..::..::..\n", 
    "COMPUTING PLOTS FOR NORMALIZED DATA:\n", sep="")

if(norm.boxplot) {
  print ("Plot boxplot for normalized intensities") 
  box.plot(normData, paste(ns, "NORM", sep="_"), col=plotColors, maxArray=50)
  }

if(norm.density) {
  print ("Plot density for normalized intensities")
  density.plot(normData, paste(ns, "NORM", sep="_"), col=plotColors, maxArray=16)
  }

if(norm.cv) {
  print ("Plot density for coefficient of varience for normalized intensities")
  cv.plot(normData, paste(ns, "NORM", sep="_"), col=plotColors, maxArray=16)
  }

if(norm.sampleRelation) {
  print ("Hierarchical clustering of normalized data")
  clusterFun(normData, normalized=TRUE,  experimentFactor=experimentFactor,
             clusterOption1=clusterOption1, clusterOption2=clusterOption2,
             plotColors=plotColors, legendColors=legendColors,
             plotSymbols=plotSymbols, legendSymbols=legendSymbols,
             WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY) 
  }


if(norm.pca) {
  print("PCA graph for normalized data")
  groupsInLegend =  !( length(unique(levels(experimentFactor))) ) >=10
  
  pcaFun(normData, normalized =TRUE, experimentFactor=experimentFactor, 
  plotColors=plotColors, legendColors=legendColors, plotSymbols=plotSymbols,
         legendSymbols=legendSymbols, groupsInLegend=groupsInLegend, namesInPlot=((max(nchar(sampleNames(normData)))<=10)&&
           (length(sampleNames(normData))<=(MAXARRAY/2))),WIDTH=WIDTH,HEIGHT=HEIGHT,
         POINTSIZE=POINTSIZE)
  }

if(norm.correl){
  print ("Correlation plot of normalized data") 
  correlFun(normData, normalized=TRUE, experimentFactor=experimentFactor, 
            clusterOption1=clusterOption1, clusterOption2=clusterOption2,
            legendColors=legendColors,
            WIDTH=WIDTH,HEIGHT=HEIGHT,POINTSIZE=POINTSIZE,MAXARRAY=MAXARRAY)
  }


###############################################################################
## Filtering                                                                 ##
###############################################################################

if(filtering) {
  # will save a filtered.normData file in the working dir.
  filtered.normData = filterFun(rawData, normData, filter.Th, filter.dp)
  }

###############################################################################
## Create annotation table                                                   ##
###############################################################################

if(createAnno) {
  
  ##raw data with anno
  cat("..::..::..\n", "MERGE ANNOTATION COLUMNS WITH 'RAW' DATA FILE.\n", sep="")

  anno.rawData = createAnnoFun(eset.rawData, ns, lib.mapping, lib.All.mapping);
  eset.anno.rawData  = cbind(anno.rawData, eset.rawData);
  
  eset.rawData.fn= paste(ns, "rawData.txt", sep='_')
  write.table(eset.anno.rawData, file = eset.rawData.fn, quote= FALSE,
              sep='\t', row.names= F, col.names= T )
  
  ##norm data with anno
  cat("MERGE ANNOTATION COLUMNS WITH 'NORMALIZED' DATA FILE.\n", sep="")
  anno.normData = createAnnoFun(eset.normData, ns, lib.mapping, lib.All.mapping);
  eset.anno.normData = cbind(anno.normData, eset.normData);
  
  eset.normData.fn= paste(ns, "normData", normType,  ".txt", sep='_')
  write.table(eset.anno.normData, file = eset.normData.fn, quote= FALSE,
              sep='\t', row.names= F, col.names= T )
}
  
if(createAnno && filtering) {
  ##filtered data with anno
  cat("MERGE ANNOTATION COLUMNS WITH 'FILTERED' DATA FILE.\n", sep="")
  anno.filtData = createAnnoFun(filtered.normData, ns, lib.mapping, lib.All.mapping);
  eset.anno.filtData = cbind(anno.filtData, filtered.normData);
  
  eset.norm.filtData.fn= paste(ns, "normData", "Filtered", normType ,".txt", sep='_')
  write.table(eset.anno.filtData, file = eset.norm.filtData.fn, quote= FALSE,
              sep='\t', row.names= F, col.names= T )
}

if(standALONE){
  savehistory(file = paste(ns, ".Rhistory"))
}

cat("..::..::..\n", 
    "ILLUMINA ANALYSIS RUN COMPLETE!\n",
    "..::..::..", "\n\n", sep="")
