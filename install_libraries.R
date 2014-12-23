#=============================================================================#
# Created by Varshna Goelela                                                  #
# Updated: 2011 July 6                                                        #
#=============================================================================#

###############################################################################
# calls function checkInstallPackages                                         #
# loads all the mandatory library packages                                    #
# Needs a character list which contains names of R packages                   #
###############################################################################

loadPackages <- function(man) {
#function to check if all mandatory packages are installed
  checkInstalledPackages(man)
  cat("..::..::..\n", 
      "Loading library packages:\n", sep="")
# for each factor (package) in fMan load package in quiet mode
  for (i in 1:length(man) ) {
    library(man[i] , character.only = TRUE, quietly=TRUE )
    print( paste (" '",man[i] ,"'" , " has been loaded", sep= "") )
  }
}

###############################################################################
# check installation of mandatory packages                                    #
# perform installation of mandatory packages                                  #
###############################################################################

checkInstalledPackages <- function(man) {
  # perform installation of mandatory packages

  # get a list of all the installed packages
  inst<-installed.packages()
    
  # check if mandatory packages are listed in installed list
  notInst<-man[!(man %in% inst[,1])]
  
  if ( length( notInst)!=0 ) {
    print("Some of the mandatory packages are missing, we are now going to install them", sep="")
    # install packages not installed
    source("http://www.bioconductor.org/biocLite.R")
    biocLite( notInst)
    # get a new list of all the installed packages
    inst<-installed.packages()
    # check for second time if mandatory packages are listed in new installed list
    notInst2 <- man[!(man %in% inst[,1])]
  
      if ( length( notInst2)==0) {
        print( paste("Ok, required package is installed and loaded: ", man, sep=" ") )
      } else {
        
        res <- ( paste("!!!Check if followwing package:", notInst2, 
                       "is available in the repository.",
                       "And try a manual install of the package", sep= " ") )
      }
  } else {
  	res<-( paste("Ok, required package installed:" , man, sep=" ") )
	}
return(res)
}
