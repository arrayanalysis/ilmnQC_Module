#=============================================================================#
# ArrayAnalysis - ILMN_functions_imagesQC                                     #
# a tool for quality control and pre-processing of Illumina array data        #                                                                        
#=============================================================================#


##############
## box.plot ##
##############

box.plot <- function(x.lumi, ns, col=NULL, maxArray=NULL) {
  
  #get nr of pages
  #divide number of sampleNames with max number of objects per page
  nrPages = ceiling(length(sampleNames(x.lumi))/maxArray)
  #divide the number of objects per page over the nr of pages
  plotsPerPage = ceiling(length(sampleNames(x.lumi))/nrPages)
  
  #pdf file name
  pdf.file = paste (ns, "boxplot.pdf", sep = "_")
  #open pdf
  pdf( paper = "a4r", file = pdf.file, onefile=T, width = 18, height = 12)
  
  for (p in 1:nrPages) {
    #create range of object for each page
    start = (p-1)*plotsPerPage+1
    end   = min(p*plotsPerPage+(p>1), length(sampleNames(x.lumi)))
    #plot function
    boxplot(x.lumi[, start:end], col=col[start:end])
  }
  #close pdf file
  dev.off()
}


####################
## density.plot   ##
####################

density.plot <- function(x.lumi, ns, col=NULL, maxArray=NULL) {
  
  #get nr of pages
  #divide number of sampleNames with max number of objects per page
  nrPages = ceiling(length(sampleNames(x.lumi))/maxArray)
  #divide the number of objects per page over the nr of pages
  plotsPerPage = ceiling(length(sampleNames(x.lumi))/nrPages)
  
  #pdf file name
  pdf.file = paste (ns, "density.pdf", sep = "_")
  #open pdf
  pdf( paper = "a4r", file = pdf.file, onefile=T, width = 18, height = 12)
  
  for (p in 1:nrPages) {
    #create range of object for each page
    start = (p-1)*plotsPerPage+1
    end   = min(p*plotsPerPage+(p>1), length(sampleNames(x.lumi)))
    #plot function
    plot(x.lumi[, start:end], col=col[start:end], what="density")
  }
  #close pdf file
  dev.off()
}


#############
## cv.plot ##
#############

cv.plot <- function(x.lumi, ns, col=NULL, maxArray=NULL) {
  
  #get nr of pages
  #divide number of sampleNames with max number of objects per page
  nrPages = ceiling(length(sampleNames(x.lumi))/maxArray)
  #divide the number of objects per page over the nr of pages
  plotsPerPage = ceiling(length(sampleNames(x.lumi))/nrPages)
  
  #pdf file name
  pdf.file = paste (ns, "cv_plot.pdf", sep = "_")
  #open pdf
  pdf( paper = "a4r", file = pdf.file, onefile=T, width = 18, height = 12)
  
  for (p in 1:nrPages) {
    #create range of object for each page
    start = (p-1)*plotsPerPage+1
    end   = min(p*plotsPerPage+(p>1), length(sampleNames(x.lumi)))
    #plot function
    plot(x.lumi[, start:end], col=col[start:end], what="cv")
  }
  #close pdf file
  dev.off()
}


#########################
## sampleRelation.plot ##
#########################

clusterFun <- function(x.lumi, normalized=FALSE ,experimentFactor=NULL, 
                       clusterOption1=clusterOption1, clusterOption2=clusterOption2, 
                       plotColors=NULL, legendColors=NULL, plotSymbols=NULL, legendSymbols=NULL, 
                       WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41) {

  if(is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified")
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")
  if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if(is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required")
  
  if(normalized){
    #normalized data
    Type="NORM"  
    main <- paste("Cluster dendrogram of",normalization.m,"normalized data")
  } else{
    #raw data
    if(normalized==FALSE) 
      Type <- "RAW"
      main <- "Cluster dendrogram of raw data"
  }
  
  if(length(sampleNames(x.lumi))<3) {
    warning("Only ",length(sampleNames(x.lumi))," sample(s) in dataset, no clustering plot made")
  } else {  
    switch(tolower(clusterOption1), 
           "pearson" = {
             correl <- cor.dist(t(exprs(x.lumi)),abs=FALSE)
           },
           "spearman" = {
             correl <- spearman.dist(t(exprs(x.lumi)),abs=FALSE)
           },
           "euclidean" = {
             correl <- euc(t(exprs(x.lumi)))
           }
           )
    clust <- hclust(correl, method = tolower(clusterOption2))
    png(file = paste(ns,"_",Type,"_","dataCluster_",clusterOption1,"_",clusterOption2,".png",sep=""),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    if(length(sampleNames(x.lumi))<MAXARRAY) {
      cexval1 <- 0.75
      cexval2 <- 1.23
      cexval3 <- 0.55
    } else {
      cexval1 <- 0.55
      cexval2 <- 1.6
      cexval3 <- 0.41
    }   
    par(cex=cexval1,oma=c(14,1,0,0))	
    par(cex.axis=cexval2,cex.lab=cexval2,cex.main=cexval2)	       
    plot(clust, hang=-1, main=main, xlab=paste("distance:",clusterOption1), sub=paste(" cluster method:",clusterOption2))
    points(1:length(clust$order),rep(0,length(clust$order)),pch=15,col="white",cex=1.5)
    points(1:length(clust$order),rep(0,length(clust$order)),pch=plotSymbols[clust$order],col=plotColors[clust$order])
    if(length(levels(experimentFactor))>1) { 
      legend("topright",levels(experimentFactor),
             pch=legendSymbols,col=legendColors)
    }
    par(cex=cexval3)    
    dev.off()
  }
}


############
## pcaFun ##
############

pcaFun <- function(x.lumi, normalized=FALSE, experimentFactor=NULL, 
                   scaled_pca=TRUE, plotColors=NULL, legendColors=NULL, 
                   plotSymbols=NULL, legendSymbols=NULL, namesInPlot=FALSE, groupsInLegend=NULL, 
                   WIDTH=1000, HEIGHT=1414, POINTSIZE=24){
  # Scaled PCA by default
  if(groupsInLegend & is.null(experimentFactor)) stop("The 'experimentFactor' parameter must be specified when groupsInLegend is true")
  if(is.null(plotColors)) stop("the 'plotColors' parameter is required")
  if(groupsInLegend & is.null(legendColors)) stop("the 'legendColors' parameter is required when groupsInLegend is true")
  if(is.null(plotSymbols)) stop("the 'plotSymbols' parameter is required")
  if(groupsInLegend & is.null(legendSymbols)) stop("the 'legendSymbols' parameter is required when groupsInLegend is true")
  
  if(groupsInLegend & !is.null(experimentFactor)) {
    if(length(levels(experimentFactor))<=1) {
      warning("groupsInLegend set to true, but no groups indicated; groups not added to legend")
    }
  }
  
  if(length(sampleNames(x.lumi))<3) {
    warning("Only",length(sampleNames(x.lumi)),"sample(s) in dataset, no PCA plot made")
  } else { 
        
    if(normalized){
      #normalized data
      Type = "NORM"
      tmain <- paste("PCA analysis of", normalization.m, "normalized data", sep=" ")
      } else {
        #raw data
        Type <- "RAW"
        tmain <- "PCA analysis of raw data"
    }    
    
    pca1 <- NULL  
    try(pca1 <- prcomp(t(exprs(x.lumi)[apply(exprs(x.lumi),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=scaled_pca),TRUE)
    if(is.null(pca1) & scaled_pca) {
      try(pca1 <- prcomp(t(exprs(Data)[apply(exprs(x.lumi),1,function(r) {sum(is.na(r))==0}),]), retx=T, center=T, scale=FALSE),TRUE)
      if(!is.null(pca1)) warning("pca with scaling unsuccessful, successfully retried without scaling")
    }
    if(!is.null(pca1)) {
      perc_expl1 <- round(((pca1$sdev[1:3]^2)/sum(pca1$sdev^2))*100,2)
      
      cex.circle <- 1.5
      cex.text <- 0.7
      cex.legend <- 0.75
      tcol <- "#444444"
      
      png(file = paste(ns,Type,"dataPCA_analysis.png",sep="_"), width=WIDTH+200*(!namesInPlot), height=HEIGHT+283*(!namesInPlot),
          pointsize=POINTSIZE)
      
      if(!namesInPlot) {
        layout(rbind(c(1,1,2,2,5),c(3,3,4,4,5)))
      } else {
        layout(rbind(c(1,1,2,2),c(1,1,2,2),c(3,3,4,4),c(3,3,4,4)))
      }
      par(oma=c(20,0,5,0))
      plot(pca1$x[,1],pca1$x[,2],cex=cex.circle,pch=plotSymbols,
           col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
           ylab=paste("PC2 (",perc_expl1[2],"%)",sep=""))
      if(namesInPlot) text(pca1$x[,1],pca1$x[,2], sampleNames(x.lumi),pos=4,cex=cex.text,col=tcol) 
      plot(pca1$x[,1],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
           col=plotColors,xlab=paste("PC1 (",perc_expl1[1],"%)",sep=""),
           ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
      if(namesInPlot) text(pca1$x[,1],pca1$x[,3], sampleNames(x.lumi),pos=4,cex=cex.text,col=tcol)
      plot(pca1$x[,2],pca1$x[,3],cex=cex.circle,pch=plotSymbols,
           col=plotColors,xlab=paste("PC2 (",perc_expl1[2],"%)",sep=""),
           ylab=paste("PC3 (",perc_expl1[3],"%)",sep=""))
      if(namesInPlot) text(pca1$x[,2],pca1$x[,3], sampleNames(x.lumi),pos=4,cex=cex.text,col=tcol)
      barplot((100*pca1$sdev^2)/sum(pca1$sdev^2),xlab="components",ylab="% of total variance explained")
      
      #determine whether groups should be added to the legend
      groupsAdded <- FALSE
      if(groupsInLegend) {
        if(length(levels(experimentFactor))>1) {
          groupsAdded <- TRUE
        }
      }
      
      if(namesInPlot) {
        if(groupsAdded) { 
          legend("topright",levels(experimentFactor),
                 pch=legendSymbols,col=legendColors,cex=cex.legend)
        }
      } else {
        par(mar=c(0,0,0,0))	
        plot(1,type="n",xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
        if(groupsAdded) {
          legend("topleft",c(levels(experimentFactor),"",sampleNames(x.lumi)),
                 #             pch=c(rep(20,length(unique(experimentFactor))+1),plotSymbols,
                 pch=c(legendSymbols,20,plotSymbols),
                 col=c(legendColors,"white",plotColors),cex=(cex.legend+0.1)
                 #             ,fill=c(legendColors,rep("white",length(experimentFactor)+1)),
                 #             border=c(legendColors,rep("white",length(experimentFactor)+1))
                 )
        } else {
          legend("topleft",sampleNames(x.lumi),pch=plotSymbols,
                 col=plotColors,cex=0.7, bty = "n")
        }
      }
      
      mtext(tmain, side = 3, outer = TRUE, font = 2, cex = 1.2)
      dev.off()
    } else {
      warning("PCA on the",Type,"data set unsuccessful, image skipped")
    }
  }
}


#################
## correl.plot ##
#################

correlFun <- function(x.lumi, normalized=FALSE, 
                      clusterOption1=clusterOption1, clusterOption2=clusterOption2,
                      experimentFactor=NULL, legendColors=NULL, 
                      WIDTH=1000, HEIGHT=1414, POINTSIZE=24,MAXARRAY=41){  
  
  if(is.null(experimentFactor)) stop("the 'exerimentFactor' parameter is required")	
  if(is.null(legendColors)) stop("the 'legendColors' parameter is required")	
    
  if(normalized){
    #normalized data
    Type="NORM"  
    text1 <- paste("Array correlation plot\nafter",normalization.m,"normalization")
  } else{
    #raw data
    Type <- "RAW"
    text1 <- "Raw data correlation plot"
  }
  
  if(length(sampleNames(x.lumi))<2) {
    warning("Only one array in dataset, no correlation plot made")
  } else {
    #print(Type)
    #print(normalized)
    png(file = paste(ns, Type,"dataArrayCorrelation.png",sep="_"),width=WIDTH,height=HEIGHT,pointsize=POINTSIZE)
    if(length(sampleNames(x.lumi))<MAXARRAY) {
      par(oma=c(17,0,0,0),cex.axis=0.7,cex.main=0.8)
      #subval <- 10
    } else {
      par(oma=c(17,0,0,0),srt=90,las=2,cex.axis=0.5,cex.main=0.8)
      #subval <- 16
    }        
    
    #note: for computing array correlation, euclidean would not make sense
    #only use euclidean distance to compute the similarity of the correlation vectors for the arrays
    COpt1 <- "pearson"
    if (tolower(clusterOption1) == "spearman") COpt1 <- "spearman"
    crp <- cor(exprs(x.lumi), use="complete.obs", method=COpt1)

    text1 <- paste(text1,"\ncorrelation method:",COpt1,"\ncluster method:",clusterOption2)
    
    switch(tolower(clusterOption1), 
           "pearson" = {
             my.dist <- function(x) cor.dist(x, abs=FALSE)
           },
           "spearman" = {
             my.dist <- function(x) spearman.dist(x, abs=FALSE)
           },
           "euclidean" = {
             my.dist <- function(x) euc(x)
           }
           )
    
    my.hclust <- function(d) hclust(d, method=clusterOption2)
    
    #in order to create some space to put colored symbols as well
    #sampleNames(x.lumi) <- paste(sampleNames(x.lumi)," ")
    
    sideColors <- legendColors[as.numeric(experimentFactor)]
    
    #cols = colorRampPalette(c("red", "yellow", "white"))(20)
    heatmap.2(crp, distfun=my.dist, hclustfun=my.hclust, trace="none", symm=TRUE, density.info="density",
              main=text1, dendrogram="row", ColSideColors=sideColors)
    
    #correlationPlot(x.lumi)    
    #axis(1,side=3,at=seq(from=0.5, to=(length(sampleNames(Data)))-0.5,by=1),
    #    labels=substr(as.character(sampleNames(x.lumi)),1,subval),las=2)
    #par(srt=0) 
    #plot(c(0,2), type = 'n', ann = FALSE, axes = FALSE, 
    #    frame.plot = FALSE, xlim = c(0, 2), ylim = c(0,2))
    #text(1,1,text1,cex=1)  
    dev.off()
  }
}
