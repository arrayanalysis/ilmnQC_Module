#=============================================================================#
# ArrayAnalysis - ILMNAnalysisQCandPreProcessing                              #
# a tool for quality control and pre-processing of Illumina array data        #                                                                        
#=============================================================================#

#####################
## Set directories ## 
#####################

### Change these Paths if needed ###
SCRIPT.DIR = "/Users/varshna/Documents/TNO/pipeline_TNO_BiGCaT/arrayanalysis_wc/trunk/src/IlmnModule/" #getwd()
# If you want to use the last version, uncomment the following line after 
# updating the milestone number:  
# SCRIPT.DIR <- "http://svn.bigcat.unimaas.nl/arrayanalysis/tags/version_1.0.0/src/"
WORK.DIR = SCRIPT.DIR 
DATA.DIR = SCRIPT.DIR
ANNO.DIR = DATA.DIR 

#amend paths of .DIRs if not started or closed off correctly with /
#-----------------------------------------------------------
correctDIR <- function(d) { 
  lastChar <- substr(d,nchar(d),nchar(d))
  if((lastChar != "/") && (lastChar != "/")) d <- paste(d,"/",sep="")
  return(d)
}

if(exists("DATA.DIR"))   DATA.DIR   <- correctDIR(DATA.DIR)
if(exists("SCRIPT.DIR")) SCRIPT.DIR <- correctDIR(SCRIPT.DIR)
if(exists("WORK.DIR"))   WORK.DIR   <- correctDIR(WORK.DIR)
if(exists("ANNO.DIR"))   ANNO.DIR   <- correctDIR(ANNO.DIR)

#change working directory into WORK.DIR
#--------------------------------------
setwd(WORK.DIR) #getwd()

version_nb <- "1.0.0"
cat("Script run using R version ",R.Version()$major,".",R.Version()$minor,
    " and Illumina pre-processing pipeline version_",version_nb,"\n",sep="")

#set memory to maximum on Windows 32bit machines
if(length(grep("w32",R.Version()$os,fixed=TRUE))>0) memory.size(4095)

###############################################################################
## Load ILMN_AnalysisQC functions and R libraries                            ##    			   
###############################################################################

reload <- function() {
  source(paste(SCRIPT.DIR,"install_libraries.R",sep=""))
  source(paste(SCRIPT.DIR,"ILMN_functions_processingQC.R",sep=""))
  source(paste(SCRIPT.DIR,"ILMN_functions_imagesQC.R",sep=""))
  cat("..::..::..\n", "FUNCTIONS HAVE BEEN LOADED.\n", sep="")
}

## reload functions from source files 
## ----------------------------------
reload();

#
#
##########################
## Set INPUT parameters ##
##########################
#
#

#Data variables:
species   <- "Human" #"Mouse" #"Rat"    

#speciesOptions <- c('Human','Mouse','Rat')
#species <- select.list(specieOptions)

arrayType <- "HumanHT-12" #"MouseRef-8
annoType  <- "HumanHT-12_V4_0_R2_15002873_B" #"MouseRef-8_V2_0_R3_11278551_A" 

#input files
#-----------
infiles = list.files(DATA.DIR)
expFile = select.list(infiles)
bgFile = select.list(infiles) #fn of the control (background) file
descFN = select.list(infiles) #REQUIRED of the desc file
#-----------

standALONE = TRUE

perGroup = TRUE #reorder rawData lumibatch file FIRST on Group and THEN ON sampleNames
 
#load description file # use this in stand-alone version
#create path to datafile
cat("..::..::..\n", 
    "LOADING DISCRIPTION FILE.\n", sep="")
descFile = paste(DATA.DIR, descFN, sep = "")
description <- read.table(descFile,
                          header=T,  
                          stringsAsFactors = F,
                          sep='\t',
                          quote="")

#lumi library packages:
lib.mapping = paste( "lumi", species, "IDMapping", sep="");
lib.All.mapping = paste( "lumi", species, "All.db", sep="");

#Importing raw data with lumi:
bgSub <- FALSE
sep = NULL
detectionTh = 0.01
na.rm = TRUE 
convertNuID = TRUE
dec = '.'
parseColumnName = FALSE
checkDupId = TRUE
rawDataQC = TRUE
rawSummary = TRUE
save.rawData = TRUE

# devault normalization function to use 
normType = "lumi"

# lumi normalization options
bg.correct = FALSE
bgcorrect.m = "bgAdjust"
variance.stabilize = TRUE
variance.m = "log2"
normalize = TRUE
normalization.m = "quantile" 
normDataQC = TRUE

normSummary = TRUE
QC.evaluation = TRUE
save.normData = TRUE

#Filtering:
filtering = TRUE
filter.Th = 0.01 #threshold is default set on <0.01
filter.dp = 0    #detect probes >0

#Annotation:
createAnno = TRUE

#raw plots
raw.boxplot = TRUE
raw.density = TRUE
raw.cv = TRUE
raw.sampleRelation = TRUE
raw.pca = TRUE
raw.correl = TRUE

#norm plots
norm.boxplot = TRUE
norm.density = TRUE
norm.cv = TRUE
norm.sampleRelation = TRUE
norm.pca = TRUE
norm.correl= TRUE

# devault correlation options
clusterOption1 = "Pearson"
clusterOption2 = "complete" #"single"   

source(paste(SCRIPT.DIR,"run_ILMNAnalysisQCandPreProcessing.R",sep=""), local=TRUE)
###########################
## PARAMETER DESCRIPTION ##
###########################

#[b]= boolean 
#[s]= string 
#[n]= numeric

#[b]= If script is run locally check TRUE

#[s] ns = General name of the study (used as prefix for output files)
#[s] species Choices: "Human", "Mouse", "Rat" 
############
#[s] arrayType Choices:
## Human: "HumanWG-6", "HumanRef-8", "HumanHT-12"
## Mouse: "MouseWG-6", "MouseRef-8"
## Rat: 	"RatRef-12"
############
#[s] arrayAnno choices:
# HumanHT-12 = "HumanHT-12_V4_0_R2_15002873_B_WGDASL", "HumanHT-12_V4_0_R2_15002873_B",
#              "HumanHT-12_V4_0_R1_15002873_B", "HumanHT-12_V3_0_R2_11283641_A",
#              "HumanHT-12_V3_0_R3_11283641_A"
# HumanRef-8 = "HumanRef-8_V3_0_R3_11282963_A", "HumanRef-8_V3_0_R2_11282963_A", 
#              "HUMANREF-8_V3_0_R1_11282963_A_WGDASL", "HumanRef-8_V2_0_R4_11223162_A")
# HumanWG-6 =  "HumanWG-6_V2_0_R4_11223189_A", "HumanWG-6_V3_0_R2_11282955_A",
#              "HumanWG-6_V3_0_R3_11282955_A",
# MouseRef-8 = "MouseRef-8_V1_1_R4_11234312_A", "MouseRef-8_V2_0_R2_11278551_A",
#              "MouseRef-8_V2_0_R3_11278551_A",
# MouseWG-6 =  "MouseWG-6_V1_1_R4_11234304_A", "MouseWG-6_V2_0_R2_11278593_A",
#              "MouseWG-6_V2_0_R3_11278593_A",
# RatRef-12 =  "RatRef-12_V1_0_R5_11222119_A"
############

#[b] perGroup = # ORDER by groupFactors THEN on sampleNames
#[b] bgSub = Is the data background subtracted in genome/bead studio
#[s] sep = The separation character used in the text file. The function can automatically determine the separation character if it is Tab or comma. 
#[n] detectionTh = The p-value threshold of determining detectability of the expression.
#[b] na.rm = Determine whether to remove NA
#[b] convertNuID = Determine whether convert the probe identifier as nuID
#[s] lib.mapping = A Illumina ID mapping package, e.g, lumiHumanIDMapping, used by addNuID2lumi
#[s] dec = The character used in the file for decimal points = "." | ","
#[b] parseColumnName = determine whether to parse the column names and retrieve the sample information (Assume the sample information is separated by "\_".)
#[b] checkDupId = Determine whether to check duplicated TargetIDs or ProbeIds. The duplicated ones will be averaged.
#[b] rawDataQC = Determine whether to do quality control assessment after read in the data; if false no summary can be computed

#[b] rawSummary = create raw summary table in work dir
#[b] save.rawData = save lumi.batch as R object 

#Pre-processing:
#[s] normType = choose a normalization type "lumi" | "neqc"
#[b] bg.correct =  decide whether to do background correction or not
#[s] bgCorrect.m = list of parameters c('none', 'bgAdjust', 'forcePositive', 'bgAdjust.affy')
#[b] variance.stabilize = decide whether to do variance stabilization or not 
#[s] variance.Stab.m = list of parameters c("vst", 'log2', 'cubicRoot')
#[b] normalize = decide whether to do normalization or not
#[s] normalization.m = list of parameters c("quantile", "rsn", "ssn", "loess", "vsn", "rankinvariant")
#[b] normDataQC = decide whether to do QC evaluation after normalization; if false no summary can be computed

#[b] normSummary = create raw summary table in work dir
#[b] save.normData = save lumi.batch as R object

#Filtering:
#[b] filtering = decide whether to do filtering on clean data
#[n] filter.Th = threshold is default set on: <0.01
#[n] filter.dp = Detect probes is default set on less stringent: >0 

#Annotation:
#[b] createAnno = whether to create a annotation table contains the lumiID, probeID, etc.


#Rawdata plots:
#[b]raw.boxplot = TRUE
#[b]raw.density = create density plot
#[b]raw.cv = create density plot
#[b]raw.sampleRelation = TRUE
#[b]raw.pca = TRUE
#[b]raw.correl = TRUE

#Normdata plots:
#[b]norm.boxplot = TRUE
#[b]norm.density = TRUE
#[b]norm.cv = TRUE
#[b]norm.sampleRelation = TRUE
#[b]norm.pca = TRUE
#[b]norm.correl = TRUE

#[s]clusterOption1 = Distance calculation method: Pearson, Spearman, Euclidean
#[s]clusterOption2 = Clustering method: Ward, McQuitty, avarage, median, single, complete, centroid

